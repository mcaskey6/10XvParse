from __future__ import annotations

from pathlib import Path
import gzip
import logging
import re
import subprocess
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
) -> None:
    '''Combine all downloaded FASTA files into one FASTQ file with splitcode'''

    logger.debug("Writing batch file for splitcode multiplexing")
    with open(batch_file, "w") as batch:
        for library in libraries:
            r1, r2 = library.gz_files
            batch.write(f"{library.name}\t{r1}\t{r2}\n")

    logger.info("Multiplexing FASTA files with splitcode")
    io.run_command(
        [
            "splitcode",
            "--remultiplex",
            "--nFastqs=2",
            "--gzip",
            "-o",
            f"{str(multiplexed_files[0])},{str(multiplexed_files[1])}",
            "--no-outb",
            str(batch_file),
            "-t",
            str(threads),
        ],
        logger,
    )


def filter_parse_fastqs(
    paths: ParsePaths,
    threads: int,
    logger: logging.Logger,
) -> None:
    '''Filter out reads that do not have the expected barcodes with splitcode'''

    logger.info("Filtering multiplexed FASTQ files with splitcode")

    with open(paths.parse_config, "w") as config_file:
        config_file.write("tags\tdistances\tids\tgroups\tminFindsG\tlocations\n")
        config_file.write(str(paths.randO_barcodes) + "\t1\tr1_R\tround1\t1\t1,78,86\n")
        config_file.write(str(paths.polyT_barcodes) + "\t1\tr1_T\tround1\t1\t1,78,86\n")

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

    logger.info("Extracting RandO reads")

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

    logger.info("Building kallisto index")
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

    logger.info("Pseudoaligning 10X multiplexed reads to genome index")
    io.run_command(
        [
            "kb",
            "count",
            "--overwrite",
            "--h5ad",
            "-t", str(threads),
            "-i", str(paths.index_file),
            "-g", str(paths.t2g_file),
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
            ["pigz", "-p", str(threads)],
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
        logger.error(seqtk_stderr)
    if pigz_stderr:
        logger.error(pigz_stderr)
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

    logger.info("Building kallisto index for 10X Hashtags")
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

    logger.info("Pseudoaligning 10X Hashtag multiplexed reads to genome index")
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

    logger.info("Building kallisto index")
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

    logger.info("Pseudoaligning Parse multiplexed reads to genome index")
    io.run_command(
        [
            "kb",
            "count",
            "--overwrite",
            "--h5ad",
            "--strand=forward",
            "--parity=single",
            "-w", str(paths.kb_onlist),
            "-t", str(threads),
            "-r", str(paths.kb_replace_config),
            "-i", str(paths.index_file),
            "-g", str(paths.t2g_file),
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

    logger.info("Subsampling FASTQ files to %d reads with seqtk", num_reads)

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
                ["pigz", "-p", str(threads)],
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
            logger.error(seqtk_stderr)
        if pigz_stderr:
            logger.error(pigz_stderr)

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
    io.download_ftp_file(era_ftp_url(err, config.r1_num), library.read1_fasta, logger)
    io.download_ftp_file(era_ftp_url(err, config.r2_num), library.read2_fasta, logger)


def make_barnyard_reference(
    paths: TenXPaths,
    config: AnalysisConfig,
    logger: logging.Logger,
) -> None:
    '''Build a barnyard (dual-species) genome reference.

    Downloads both species' FASTA and GTF files, prefixes chromosome names and
    gene identifiers with the species label (e.g. HUMAN_, MOUSE_), then
    concatenates into a single index-ready pair at paths.genome_file / paths.gtf_file.
    The species field in config must be formatted as "{sp1}_{sp2}" (e.g. "human_mouse").
    '''
    if not config.genome_url_2:
        raise ValueError("make_barnyard_reference requires genome_url_2 and gtf_url_2 in AnalysisConfig")

    sp1, sp2 = config.species.split("_", 1)
    prefix_1, prefix_2 = sp1.upper(), sp2.upper()

    tmp = paths.tmp_dir
    fa_1  = tmp / f"{sp1}.fa.gz"
    fa_2  = tmp / f"{sp2}.fa.gz"
    gtf_1 = tmp / f"{sp1}.gtf.gz"
    gtf_2 = tmp / f"{sp2}.gtf.gz"

    logger.info("Downloading %s reference", sp1)
    io.download_file(config.genome_url, fa_1, logger)
    io.download_file(config.gtf_url, gtf_1, logger)

    logger.info("Downloading %s reference", sp2)
    io.download_file(config.genome_url_2, fa_2, logger)
    io.download_file(config.gtf_url_2, gtf_2, logger)

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
        )
    else:
        logger.info(
            "%s and %s already exist. Skipping FASTA file multiplexing with splitcode.",
            paths.multiplexed_files[0],
            paths.multiplexed_files[1],
        )


def core_pipeline(settings: RunSettings, paths: BasePaths, config: AnalysisConfig, assay: str, logger: logging.Logger) -> None:
    '''Download SRA reads and multiplex into a single paired FASTQ. Reference download is handled by the caller.'''

    logger.info(f"Starting {assay} pipeline")

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

    _multiplex_into_fastq(settings, paths, libraries, logger)


def era_core_pipeline(
    settings: RunSettings,
    paths: BasePaths,
    config: AnalysisConfig,
    assay: str,
    logger: logging.Logger,
) -> None:
    '''Download ERA reads from ENA FTP and multiplex into a single paired FASTQ.'''

    logger.info("Starting %s ERA pipeline", assay)

    libraries = io.build_era_libraries(config, paths)

    for err, library in zip(config.era, libraries):
        downloaded_exist = all(p.is_file() for p in library.gz_files)
        if not downloaded_exist or settings.overwrite:
            download_era_fastq(err, library, config, logger)
        else:
            logger.info("Files for %s already downloaded. Skipping.", err)

    _multiplex_into_fastq(settings, paths, libraries, logger)


def local_pipeline(
    settings: RunSettings,
    paths: BasePaths,
    config: AnalysisConfig,
    assay: str,
    logger: logging.Logger,
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
    _multiplex_into_fastq(settings, paths, libraries, logger)


def _count_reads(fastq_gz: Path) -> int:
    '''Count reads in a gzipped FASTQ by dividing line count by 4.'''
    zcat = subprocess.Popen(["zcat", str(fastq_gz)], stdout=subprocess.PIPE)
    wc = subprocess.Popen(["wc", "-l"], stdin=zcat.stdout, stdout=subprocess.PIPE)
    assert zcat.stdout is not None
    zcat.stdout.close()
    out, _ = wc.communicate()
    zcat.wait()
    return int(out.strip()) // 4
