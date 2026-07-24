from __future__ import annotations

from pathlib import Path
import gzip
import logging
import re
import subprocess
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from .classes import LibraryFiles, BasePaths, TenXPaths, HashtagsPaths, ParsePaths, AnalysisConfig, RunSettings
from . import io


def get_reference(
    genome_file: Path,
    gtf_file: Path,
    genome_url: str,
    gtf_url: str,
    logger: logging.Logger,
) -> None:
    '''Download genome reference given the urls to the fasta and gtf files'''
    logger.info("Downloading genome reference files")
    io.download_file(genome_url, genome_file, logger)
    io.download_file(gtf_url, gtf_file, logger)
    
    logger.debug("Unzipping genome GTF and FASTA files for STAR")
    io.run_command(
        ["gunzip", "-k", str(genome_file)],
        logger=logger
    )
    io.run_command(
        ["gunzip", "-k", str(gtf_file)],
        logger=logger
    )

def prefetch_one_sra(srr: str, paths: BasePaths, logger: logging.Logger) -> None:
    """Download one SRA file with prefetch."""
    logger.info("Prefetching %s to %s", srr, paths.sra_dir)
    io.run_command(
        [
            "prefetch",
            srr,
            "--max-size",
            "u",
            "-O",
            str(paths.sra_dir),
        ],
        logger,
    )


def prefetch_sra(
    srrs: list[str],
    paths: BasePaths,
    logger: logging.Logger,
    max_workers: int = 2,
) -> None:
    """Download multiple SRA files concurrently with prefetch."""
    logger.info("Starting prefetch for %d SRR accession(s)", len(srrs))

    failures: list[tuple[str, Exception]] = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_srr = {
            executor.submit(prefetch_one_sra, srr, paths, logger): srr
            for srr in srrs
        }

        for future in as_completed(future_to_srr):
            srr = future_to_srr[future]
            try:
                future.result()
                logger.info("Finished prefetch for %s", srr)
            except Exception as e:
                logger.exception("Prefetch failed for %s", srr)
                failures.append((srr, e))

    if failures:
        failed_srrs = ", ".join(srr for srr, _ in failures)
        raise RuntimeError(f"Prefetch failed for: {failed_srrs}")


def dump_sra(
    srr: str,
    library: LibraryFiles,
    paths: BasePaths,
    threads: int,
    logger: logging.Logger,
) -> None:
    '''Download SRA files from the cloud using srr accession number, then convert to FASTA.'''

    logger.info("Dumping %s to %s", srr, paths.dumped_dir / library.name)

    ## Dump FASTA files from SRA files
    logger.debug("Running fasterq-dump for %s", srr)
    io.run_command(
        [
            "fasterq-dump",
            "--outdir",
            str(paths.dumped_dir),
            "--temp",
            str(paths.tmp_dir),
            "--outfile",
            f"{library.name}.fasta",
            "--split-files",
            "--skip-technical",
            "-f",
            str(paths.sra_dir / srr / f"{srr}.sra"),
            "--threads",
            str(threads),
            "--fasta-unsorted",
        ],
        logger,
    )

    logger.debug("Zipping %s read 1 and read 2", library.name)
    io.run_command(["pigz", "-f", str(library.read1_fasta), "-p", str(threads)], logger)
    io.run_command(["pigz", "-f", str(library.read2_fasta), "-p", str(threads)], logger)


def multiplex_fastqs(
    multiplexed_files: list[Path],
    batch_file: Path,
    libraries: list[LibraryFiles],
    threads: int,
    logger: logging.Logger,
    bclen: int | None = None,
) -> None:
    '''Combine all downloaded FASTA files into one FASTQ file with splitcode

    With `bclen` set, splitcode prepends a synthetic per-sublibrary barcode of that length,
    which is what lets Parse sublibraries share one FASTQ without their cells colliding.
    splitcode embeds the barcode in the first output file, so the file order is swapped to
    put the barcode read first; multiplexed_files itself keeps its cDNA-first order.
    Libraries whose batch name matches share a barcode.
    '''

    # splitcode assigns barcodes in order of first appearance in the batch file, so the
    # barcode a sublibrary receives depends on the order libraries are listed here.
    barcode_read_first = bclen is not None

    logger.debug("Writing batch file for splitcode multiplexing")
    with open(batch_file, "w") as batch:
        for library in libraries:
            r1, r2 = library.gz_files
            first, second = (r2, r1) if barcode_read_first else (r1, r2)
            batch.write(f"{library.batch_name}\t{first}\t{second}\n")

    outputs = list(reversed(multiplexed_files)) if barcode_read_first else multiplexed_files

    command = [
        "splitcode",
        "--remultiplex",
        "--nFastqs=2",
        "--gzip",
        "-o",
        f"{str(outputs[0])},{str(outputs[1])}",
    ]
    if bclen is None:
        # Without a sublibrary barcode there is nothing to embed or write out
        command.append("--no-outb")
    else:
        command.append(f"--bclen={bclen}")
    command += [str(batch_file), "-t", str(threads)]

    logger.info("Multiplexing %d libraries -> %s (splitcode)", len(libraries), multiplexed_files[0].parent)
    io.run_command(command, logger)


def filter_parse_fastqs(
    paths: ParsePaths,
    threads: int,
    logger: logging.Logger,
    bc1_location: str = "1,78,86",
) -> None:
    '''Filter out reads that do not have the expected barcodes with splitcode'''

    logger.info("Filtering %s by Parse barcodes (splitcode)", paths.multiplexed_files[0].name)

    with open(paths.parse_config, "w") as config_file:
        config_file.write("tags\tdistances\tids\tgroups\tminFindsG\tlocations\n")
        config_file.write(str(paths.randO_barcodes) + f"\t1\tr1_R\tround1\t1\t{bc1_location}\n")
        config_file.write(str(paths.polyT_barcodes) + f"\t1\tr1_T\tround1\t1\t{bc1_location}\n")

    with open(paths.parse_keep_file, "w") as keep_file:
        keep_file.write(f"round1 {str(paths.filtered_files[0]).split('_0')[0]}")

    io.run_command(
        [
            "splitcode",
            "-c", str(paths.parse_config),
            "--keep-grp", str(paths.parse_keep_file),
            "--nFastqs=2",
            "--gzip",
            "--no-output",
            "--no-outb",
            str(paths.multiplexed_files[0]),
            str(paths.multiplexed_files[1]),
            "-t", str(threads),
        ],
        logger,
    )


def extract_rando_polyt(
    paths: ParsePaths,
    threads: int,
    logger: logging.Logger,
) -> None:
    '''Generate FASTQ file of randO reads with splitcode'''

    logger.info("Splitting %s into polyT and randO subsets", paths.filtered_files[0].name)

    with open(paths.randOpolyT_keep_file, "w") as keep_file:
        keep_file.write(f"r1_R {str(paths.randO_files[0]).split('_0')[0]}\n")
        keep_file.write(f"r1_T {str(paths.polyT_files[0]).split('_0')[0]}")

    io.run_command([
        "splitcode",
        "-c", str(paths.parse_config),
        "--keep", str(paths.randOpolyT_keep_file),
        "--nFastqs=2",
        "--gzip",
        "--no-output",
        "--no-outb",
        str(paths.filtered_files[0]),
        str(paths.filtered_files[1]),
        "-t", str(threads)
    ], logger)


def pseudoalign_10x(
    paths: TenXPaths,
    fastq_files: list[Path],
    kb_out_dir: Path,
    tech: str,
    threads: int,
    logger: logging.Logger,
) -> None:
    '''Psuedoalign multiplexed files to reference. Build index if needed'''

    logger.info("Building kallisto index at %s", paths.index_dir)
    io.run_command(
        [
            "kb",
            "ref",
            "--workflow", "nac",
            "-i", str(paths.index_file),
            "-g", str(paths.t2g_file),
            "-c1", str(paths.cdna_file),
            "-c2", str(paths.nascent_file),
            "-f1", str(paths.cdna_fasta_file),
            "-f2", str(paths.nascent_fasta_file),
            str(paths.genome_file),
            str(paths.gtf_file),
        ],
        logger,
    )

    logger.info("Pseudoaligning 10x reads -> %s", kb_out_dir)
    io.run_command(
        [
            "kb",
            "count",
            "--overwrite",
            "--h5ad",
            "--workflow=nac",
            "-t", str(threads),
            "-i", str(paths.index_file),
            "-g", str(paths.t2g_file),
            "-c1", str(paths.cdna_file),
            "-c2", str(paths.nascent_file),
            "-x", tech,
            "-o", str(kb_out_dir),
            str(fastq_files[0]),
            str(fastq_files[1]),
        ],
        logger,
    )


def _trim_r2(fastq_r2: Path, trim_length: int, threads: int, logger: logging.Logger) -> Path:
    '''Trim R2 FASTQ to trim_length bases with seqtk, writing compressed output alongside the input.'''
    trimmed = fastq_r2.parent / fastq_r2.name.replace(".fastq.gz", "_trimmed.fastq.gz")
    logger.info("Trimming R2 to %d bp -> %s", trim_length, trimmed)

    with open(trimmed, "wb") as out_f:
        seqtk = subprocess.Popen(
            ["seqtk", "trimfq", "-L", str(trim_length), str(fastq_r2)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        pigz = subprocess.Popen(
            ["pigz", "-f", "-p", str(threads)],
            stdin=seqtk.stdout,
            stdout=out_f,
            stderr=subprocess.PIPE,
        )
        assert seqtk.stdout is not None
        seqtk.stdout.close()
        pigz_stderr = pigz.stderr.read().decode() if pigz.stderr else ""
        pigz_rc = pigz.wait()
        seqtk_stderr = seqtk.stderr.read().decode() if seqtk.stderr else ""
        seqtk_rc = seqtk.wait()

    if seqtk_stderr:
        logger.error("seqtk error trimming %s: %s", fastq_r2, seqtk_stderr)
    if pigz_stderr:
        logger.error("pigz error compressing %s: %s", trimmed, pigz_stderr)
    if seqtk_rc != 0:
        raise subprocess.CalledProcessError(seqtk_rc, seqtk.args)
    if pigz_rc != 0:
        raise subprocess.CalledProcessError(pigz_rc, pigz.args)

    return trimmed


def pseudoalign_10x_hashtags(
    paths: HashtagsPaths,
    fastq_files: list[Path],
    kb_out_dir: Path,
    tech: str,
    threads: int,
    logger: logging.Logger,
    trim_length: int | None = None,
) -> None:
    '''Psuedoalign multiplexed files to reference. Build index if needed'''

    logger.info("Building kallisto index for 10X Hashtags at %s", paths.index_dir)
    io.run_command(
        [
            "kb",
            "ref",
            "--workflow", "kite",
            "--overwrite",
            "-i", str(paths.index_file),
            "-g", str(paths.t2g_file),
            "-f1", str(paths.cdna_file),
            str(paths.genome_file),
        ],
        logger,
    )

    r2 = fastq_files[1]
    if trim_length is not None:
        r2 = _trim_r2(fastq_files[1], trim_length, threads, logger)

    logger.info("Pseudoaligning 10x Hashtag reads -> %s", kb_out_dir)
    io.run_command(
        [
            "kb",
            "count",
            "--workflow", "kite",
            "--overwrite",
            "--h5ad",
            "-t", str(threads),
            "-i", str(paths.index_file),
            "-g", str(paths.t2g_file),
            "-x", tech,
            "-o", str(kb_out_dir),
            str(fastq_files[0]),
            str(r2),
        ],
        logger,
    )


def pseudoalign_parse(
    paths: ParsePaths,
    fastq_files: list[Path],
    kb_out_dir: Path,
    tech: str,
    threads: int,
    logger: logging.Logger,
) -> None:
    '''Pseudoalign all parse-specific files (FASTQ with all parse reads, polyT FASTQ, randO FASTQ)'''

    logger.info("Building kallisto index at %s", paths.index_dir)
    io.run_command(
        [
            "kb",
            "ref",
            "--workflow", "nac",
            "-i", str(paths.index_file),
            "-g", str(paths.t2g_file),
            "-c1", str(paths.cdna_file),
            "-c2", str(paths.nascent_file),
            "-f1", str(paths.cdna_fasta_file),
            "-f2", str(paths.nascent_fasta_file),
            str(paths.genome_file),
            str(paths.gtf_file),
        ],
        logger,
    )

    logger.info("Pseudoaligning Parse reads [%s] -> %s", fastq_files[0].name, kb_out_dir)
    io.run_command(
        [
            "kb",
            "count",
            "--overwrite",
            "--h5ad",
            "--workflow=nac",
            "--strand=forward",
            "--parity=single",
            "-w", str(paths.kb_onlist),
            "-t", str(threads),
            "-r", str(paths.kb_replace_config),
            "-i", str(paths.index_file),
            "-g", str(paths.t2g_file),
            "-c1", str(paths.cdna_file),
            "-c2", str(paths.nascent_file),
            "-x", tech,
            "-o", str(kb_out_dir),
            str(fastq_files[0]),
            str(fastq_files[1]),
        ],
        logger,
    )


def subsample_fastqs(
    fastq_files: list[Path],
    output_files: list[Path],
    num_reads: int,
    threads: int,
    logger: logging.Logger
) -> None:
    # Subsample FASTQ files to a specified number of reads with seqtk

    logger.info("Subsampling %s -> %s (%d reads)", fastq_files[0].name, output_files[0].name, num_reads)

    def subsample_one_fastq(
        input_file: Path,
        output_file: Path,
    ) -> None:

        with open(output_file, "wb") as out_f:
            seqtk = subprocess.Popen(
                ["seqtk", "sample", "-s", "42", str(input_file), str(num_reads)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            pigz = subprocess.Popen(
                ["pigz", "-f", "-p", str(threads)],
                stdin=seqtk.stdout,
                stdout=out_f,
                stderr=subprocess.PIPE,
            )

            # Let seqtk receive SIGPIPE properly if pigz exits early
            assert seqtk.stdout is not None
            seqtk.stdout.close()

            seqtk_stderr = seqtk.stderr.read().decode() if seqtk.stderr else ""
            pigz_stderr = pigz.stderr.read().decode() if pigz.stderr else ""

            seqtk_rc = seqtk.wait()
            pigz_rc = pigz.wait()

        if seqtk_stderr:
            logger.error("seqtk error sampling %s: %s", input_file, seqtk_stderr)
        if pigz_stderr:
            logger.error("pigz error compressing %s: %s", output_file, pigz_stderr)

        if seqtk_rc != 0:
            raise subprocess.CalledProcessError(seqtk_rc, seqtk.args)
        if pigz_rc != 0:
            raise subprocess.CalledProcessError(pigz_rc, pigz.args)

    subsample_one_fastq(fastq_files[0], output_files[0])
    subsample_one_fastq(fastq_files[1], output_files[1])


def era_ftp_url(err: str, read_num: int) -> str:
    '''Compute the ENA FTP URL for a given run accession and read number.'''
    prefix = err[:6]
    n = len(err)
    if n <= 9:
        path = f"{prefix}/{err}"
    elif n == 10:
        path = f"{prefix}/00{err[-1]}/{err}"
    elif n == 11:
        path = f"{prefix}/0{err[-2:]}/{err}"
    else:
        path = f"{prefix}/{err[-3:]}/{err}"
    return f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{path}/{err}_{read_num}.fastq.gz"


def download_era_fastq(
    err: str,
    library: LibraryFiles,
    config: AnalysisConfig,
    logger: logging.Logger,
) -> None:
    '''Download one ERA run's paired FASTQ files from ENA FTP.'''
    logger.info("Downloading ERA run %s", err)
    r1_num, r2_num = config.read_nums_for(library.sublibrary)
    io.download_ftp_file(era_ftp_url(err, r1_num), library.read1_fasta, logger)
    io.download_ftp_file(era_ftp_url(err, r2_num), library.read2_fasta, logger)


def make_barnyard_reference(
    paths: TenXPaths,
    index_config,
    logger: logging.Logger,
) -> None:
    '''Build a barnyard (dual-species) genome reference.

    Downloads both species' FASTA and GTF files, prefixes chromosome names and
    gene identifiers with the species label (e.g. HUMAN_, MOUSE_), then
    concatenates into a single index-ready pair at paths.genome_file / paths.gtf_file.
    The species field in index_config must be formatted as "{sp1}_{sp2}" (e.g. "human_mouse").
    '''
    if not index_config.fasta_url_2:
        raise ValueError("make_barnyard_reference requires fasta_url_2 and gtf_url_2 in IndexConfig")

    sp1, sp2 = index_config.species.split("_", 1)
    prefix_1, prefix_2 = sp1.upper(), sp2.upper()

    tmp = paths.tmp_dir
    fa_1  = tmp / f"{sp1}.fa.gz"
    fa_2  = tmp / f"{sp2}.fa.gz"
    gtf_1 = tmp / f"{sp1}.gtf.gz"
    gtf_2 = tmp / f"{sp2}.gtf.gz"

    logger.info("Downloading %s reference", sp1)
    io.download_file(index_config.fasta_url, fa_1, logger)
    io.download_file(index_config.gtf_url, gtf_1, logger)

    logger.info("Downloading %s reference", sp2)
    io.download_file(index_config.fasta_url_2, fa_2, logger)
    io.download_file(index_config.gtf_url_2, gtf_2, logger)

    logger.info("Building combined barnyard FASTA")
    with gzip.open(paths.genome_file, "wt") as out:
        for fa_in, prefix in [(fa_1, prefix_1), (fa_2, prefix_2)]:
            with gzip.open(fa_in, "rt") as f:
                for line in f:
                    if line.startswith(">"):
                        seq_id, _, rest = line[1:].partition(" ")
                        out.write(f">{prefix}_{seq_id} {rest}" if rest else f">{prefix}_{seq_id}\n")
                    else:
                        out.write(line)

    logger.info("Building combined barnyard GTF")
    with gzip.open(paths.gtf_file, "wt") as out:
        for gtf_in, prefix in [(gtf_1, prefix_1), (gtf_2, prefix_2)]:
            with gzip.open(gtf_in, "rt") as f:
                for line in f:
                    if line.startswith("#"):
                        out.write(line)
                    else:
                        cols = line.split("\t", 1)
                        line = f"{prefix}_{cols[0]}\t{cols[1]}"
                        line = re.sub(r'(gene_id ")([^"]+)(")', rf'\1{prefix}_\2\3', line)
                        line = re.sub(r'(gene_name ")([^"]+)(")', rf'\1{prefix}_\2\3', line)
                        out.write(line)


def _multiplex_into_fastq(
    settings: RunSettings,
    paths: BasePaths,
    libraries: list[LibraryFiles],
    logger: logging.Logger,
    bclen: int | None = None,
) -> None:
    '''Multiplex pre-downloaded library files into a single paired FASTQ with splitcode.'''
    processed_exist = all(p.is_file() for p in paths.multiplexed_files)
    if not processed_exist or settings.overwrite:
        multiplex_fastqs(
            multiplexed_files=paths.multiplexed_files,
            batch_file=paths.batch_file,
            libraries=libraries,
            threads=settings.threads,
            logger=logger,
            bclen=bclen,
        )
    else:
        logger.info(
            "%s and %s already exist. Skipping FASTA file multiplexing with splitcode.",
            paths.multiplexed_files[0],
            paths.multiplexed_files[1],
        )


def core_pipeline(
    settings: RunSettings,
    paths: BasePaths,
    config: AnalysisConfig,
    assay: str,
    logger: logging.Logger,
    bclen: int | None = None,
) -> None:
    '''Download SRA reads and multiplex into a single paired FASTQ. Reference download is handled by the caller.'''

    logger.info("[%s/%s] Starting SRA pipeline", config.name, assay)

    libraries = io.build_libraries(config, paths)

    prefetch_sra(config.sra, paths, logger, settings.max_workers)

    for srr, library in zip(config.sra, libraries):
        dumped_exist = all(p.is_file() for p in library.gz_files)
        if not dumped_exist or settings.overwrite:
            dump_sra(
                srr=srr,
                library=library,
                paths=paths,
                threads=settings.threads,
                logger=logger,
            )
        else:
            logger.info(
                "Files for %s have already been dumped and zipped. Skipping.",
                library.name,
            )

    _multiplex_into_fastq(settings, paths, libraries, logger, bclen)


def era_core_pipeline(
    settings: RunSettings,
    paths: BasePaths,
    config: AnalysisConfig,
    assay: str,
    logger: logging.Logger,
    bclen: int | None = None,
) -> None:
    '''Download ERA reads from ENA FTP and multiplex into a single paired FASTQ.'''

    logger.info("[%s/%s] Starting ERA pipeline", config.name, assay)

    libraries = io.build_era_libraries(config, paths)

    for err, library in zip(config.era, libraries):
        downloaded_exist = all(p.is_file() for p in library.gz_files)
        if not downloaded_exist or settings.overwrite:
            download_era_fastq(err, library, config, logger)
        else:
            logger.info("Files for %s already downloaded. Skipping.", err)

    _multiplex_into_fastq(settings, paths, libraries, logger, bclen)


def local_pipeline(
    settings: RunSettings,
    paths: BasePaths,
    config: AnalysisConfig,
    assay: str,
    logger: logging.Logger,
    bclen: int | None = None,
) -> None:
    '''Use pre-existing FASTQ files from dumped_dir when no SRA/ERA accessions are provided.'''
    logger.info("No accessions in config for %s — checking for local files in %s", assay, paths.dumped_dir)
    libraries = io.build_local_libraries(config, paths)
    if not libraries:
        raise FileNotFoundError(
            f"No SRA/ERA accessions in config and no local files found in {paths.dumped_dir}. "
            f"Expected files named Lib0_{config.r1_num}.fastq.gz / Lib0_{config.r2_num}.fastq.gz, etc."
        )
    missing = [p for lib in libraries for p in lib.gz_files if not p.is_file()]
    if missing:
        raise FileNotFoundError(f"Local library files missing: {missing}")
    _multiplex_into_fastq(settings, paths, libraries, logger, bclen)

def count_reads(fastq_gz: Path) -> int:
    '''Count reads in a gzipped FASTQ by dividing line count by 4.'''
    zcat = subprocess.Popen(["zcat", str(fastq_gz)], stdout=subprocess.PIPE)
    wc = subprocess.Popen(["wc", "-l"], stdin=zcat.stdout, stdout=subprocess.PIPE)
    assert zcat.stdout is not None
    zcat.stdout.close()
    out, _ = wc.communicate()
    zcat.wait()
    return int(out.strip()) // 4

def run_star_10x(
    paths: TenXPaths,
    settings: RunSettings,
    assay: str,
    logger: logging.Logger,
    tag: str = None,
    overwrite: bool = False,
) -> str:
    '''Run STARsolo on 10X data to generate a BAM file for downstream 
    gene body coverage analysis. Returns the outfile prefix.
    
    Optional tag prefix for analyses with multiple dataset combinations.'''

    if tag:
        sampled_dir = paths.fasta_dir / f"Sampled_{tag}"
        sampled_files = [sampled_dir / f"{assay}_{i}.fastq.gz" for i in range(2)]
        outfile_prefix = str(paths.star_dir / tag) + "/10x_"
    else:
        sampled_files = paths.sampled_files
        outfile_prefix = str(paths.star_dir) + "/10x_"

    if not (paths.star_index_dir / "genomeParameters.txt").is_file():
        logger.info("Building STAR index at %s", paths.star_index_dir)
        io.run_command(
            [
                "STAR",
                "--runThreadN", str(settings.threads),
                "--runMode", "genomeGenerate",
                "--genomeDir", str(paths.star_index_dir),
                "--genomeFastaFiles", str(paths.genome_file).removesuffix(".gz"),
                "--sjdbGTFfile", str(paths.gtf_file).removesuffix(".gz"),
            ],
            logger
        )
    else:
        logger.info("STAR index already exists at %s. Skipping build.", paths.star_index_dir)

    output = Path(outfile_prefix + "Aligned.sortedByCoord.out.bam")
    if not output.is_file() or not output.stat().st_size > 0 or overwrite:
        logger.info("Running STARsolo for %s", assay)
        io.run_command(
            [
                "STAR",
                "--soloType", "CB_UMI_Simple",
                "--soloCBwhitelist", str(paths.kb_onlist),
                "--soloBarcodeReadLength", "0",
                "--runThreadN", str(settings.threads),
                "--genomeDir", str(paths.star_index_dir),
                "--outFileNamePrefix", outfile_prefix,
                "--readFilesIn", str(sampled_files[1]), str(sampled_files[0]),
                "--readFilesCommand", "zcat",
                "--outSAMtype", "BAM", "SortedByCoordinate",
            ],
            logger
        )
    else:
        logger.info("Output BAM already for %s exists. Skipping alignment with STAR.", assay)

    return outfile_prefix

def run_star_parse(
    paths: ParsePaths,
    config: AnalysisConfig,
    settings: RunSettings,
    assay: str,
    fastq_files: list[Path],
    tag: str,
    x_string: str,
    logger: logging.Logger,
    overwrite: bool = False,
) -> str:
    '''Run STARsolo (CB_UMI_Complex) on Parse data to generate a BAM for downstream
    gene body coverage analysis.

    Call once per read subset, passing the appropriate fastq_files and a tag that becomes
    the filename prefix (e.g. "all", "polyT", "randO"). `x_string` must be the same
    (sublibrary-shifted) technology string used for pseudoalignment. Returns the outfile prefix
    '''
    from . import parse_config as pc

    if not (paths.star_index_dir / "genomeParameters.txt").is_file():
        logger.info("Building STAR index at %s", paths.star_index_dir)
        io.run_command(
            [
                "STAR",
                "--runThreadN", str(settings.threads),
                "--runMode", "genomeGenerate",
                "--genomeDir", str(paths.star_index_dir),
                "--genomeFastaFiles", str(paths.genome_file).removesuffix(".gz"),
                "--sjdbGTFfile", str(paths.gtf_file).removesuffix(".gz"),
            ],
            logger
        )
    else:
        logger.info("STAR index already exists at %s. Skipping build.", paths.star_index_dir)

    cb_positions, umi_position = pc.x_string_to_star_params(x_string)
    configs_dir = paths.config_dir / config.name / assay

    # One whitelist per --soloCBposition, in x_string order. The sublibrary barcode comes
    # first when the reads were remultiplexed with one.
    whitelists = [configs_dir / f"star_bc{i}.txt" for i in (3, 2, 1)]
    if len(cb_positions) == len(whitelists) + 1:
        whitelists.insert(0, configs_dir / "lib_bc.txt")

    outfile_prefix = str(paths.star_dir) + f"/{tag}_"

    output = Path(outfile_prefix + "Aligned.sortedByCoord.out.bam")
    if not output.is_file() or not output.stat().st_size > 0 or overwrite:
        logger.info("Running STARsolo for %s [%s]", assay, tag)
        io.run_command(
            [
                "STAR",
                "--soloType", "CB_UMI_Complex",
                "--soloCBwhitelist", *[str(w) for w in whitelists],
                "--soloCBposition", *cb_positions,
                "--soloUMIposition", umi_position,
                "--soloBarcodeReadLength", "0",
                "--runThreadN", str(settings.threads),
                "--genomeDir", str(paths.star_index_dir),
                "--outFileNamePrefix", outfile_prefix,
                "--readFilesIn", str(fastq_files[0]), str(fastq_files[1]),
                "--readFilesCommand", "zcat",
                "--soloCBmatchWLtype", "1MM",
                "--outSAMtype", "BAM", "SortedByCoordinate",
            ],
            logger
        )
    else:
        logger.info("Output BAM for %s [%s] already exists. Skipping alignment with STAR.", assay, tag)

    return outfile_prefix

def download_hk_genes(
    hk_genes_file: Path,
    hrt_atlas_url: str,
    logger: logging.Logger,
    hrt_atlas_url_2: str = None,
) -> None:
    """Download HRT Atlas housekeeping transcript IDs and cache to hk_genes_file.

    Accepts one URL (single-species) or two (barnyard). The CSV files use either
    ',' or ';' as delimiter; transcript IDs are always in column 0. Version suffixes
    are stripped so IDs match gffread BED output.
    """
    if hk_genes_file.is_file():
        return

    ids: set[str] = set()
    for url in filter(None, [hrt_atlas_url, hrt_atlas_url_2]):
        logger.info("Downloading HRT Atlas housekeeping gene list from %s", url)
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        text = response.text
        sep = ";" if ";" in text.splitlines()[0] else ","
        for line in text.splitlines()[1:]:
            if line:
                tid = line.split(sep)[0].strip().split(".")[0]
                if tid:
                    ids.add(tid)

    with open(hk_genes_file, "w") as f:
        for tid in sorted(ids):
            f.write(tid + "\n")

    logger.info("Cached %d housekeeping transcript IDs to %s", len(ids), hk_genes_file)


def filter_bed_to_housekeeping(
    bed_file: Path,
    hk_genes_file: Path,
    hk_bed_file: Path,
    logger: logging.Logger,
) -> None:
    """Subset a gffread BED12 file to transcripts in the HRT Atlas housekeeping list.

    Matches on column 4 (transcript ID), stripping any HUMAN_/MOUSE_ barnyard prefix
    and version suffix before lookup.
    """
    hk_ids: set[str] = set()
    with open(hk_genes_file) as f:
        for line in f:
            hk_ids.add(line.strip())

    total = kept = 0
    with open(bed_file) as f_in, open(hk_bed_file, "w") as f_out:
        for line in f_in:
            total += 1
            cols = line.split("\t")
            if len(cols) > 3:
                tid = cols[3]
                for prefix in ("HUMAN_", "MOUSE_"):
                    if tid.startswith(prefix):
                        tid = tid[len(prefix):]
                        break
                tid = tid.split(".")[0]
                if tid in hk_ids:
                    f_out.write(line)
                    kept += 1

    logger.info("Filtered BED to %d/%d housekeeping transcripts", kept, total)


def generate_genebody_plot(
    settings: RunSettings,
    gtf_file: Path,
    bed_file: Path,
    bam_files: list[Path],
    gene_body_dir: Path,
    logger: logging.Logger,
    hk_genes_file: Path = None,
    hrt_atlas_url: str = None,
    hrt_atlas_url_2: str = None,
    hk_bed_file: Path = None,
    tag: str = None,
) -> None:
    """From the STARsolo generated BAM files generate a genebody coverage plot
    with geneBody_coverage.py from RSeQC.

    When hk_bed_file is provided along with hrt_atlas_url, the coverage plot is
    restricted to HRT Atlas housekeeping transcripts, which is much faster than
    running on the full genome BED.

    The tag field is used for analyses with multiple combinations of datasets to
    specify the particular set of datasets compared.
    """
    logger.info("Generating BED file from GTF")
    if not bed_file.is_file() or settings.overwrite:
        with open(bed_file, "w") as bed_out:
            subprocess.run(
                ["gffread", str(gtf_file).removesuffix(".gz"), "--bed"],
                stdout=bed_out,
                check=True,
            )

    if hk_bed_file is not None and hrt_atlas_url:
        download_hk_genes(hk_genes_file, hrt_atlas_url, logger, hrt_atlas_url_2)
        if not hk_bed_file.is_file() or settings.overwrite:
            filter_bed_to_housekeeping(bed_file, hk_genes_file, hk_bed_file, logger)
        bed_file = hk_bed_file

    for bam in bam_files:
        # Re-index whenever the BAM is newer than its index: STAR rewrites BAMs on every
        # run_kb pass, so an existing .bai is not evidence that it matches the current BAM.
        bai = Path(str(bam) + ".bai")
        stale = bai.is_file() and bai.stat().st_mtime < Path(bam).stat().st_mtime
        if not bai.is_file() or stale or settings.overwrite:
            io.run_command(["samtools", "index", str(bam)], logger=logger)

    if not tag:
        outdir = str(gene_body_dir) + "/"
    else:
        outdir = str(gene_body_dir) + f"/{tag}"

    io.run_command(
        [
            "geneBody_coverage.py",
            "-r", str(bed_file),
            "-i", f"{bam_files[0]},{bam_files[1]},{bam_files[2]},{bam_files[3]}",
            "-o", outdir
        ],
        logger=logger,
    )
