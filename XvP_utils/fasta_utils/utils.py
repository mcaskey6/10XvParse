from __future__ import annotations

from pathlib import Path
import logging
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from .classes import LibraryFiles, PipelinePaths, AnalysisConfig, RunSettings
from . import helpers

def get_reference(
    genome_file: Path,
    gtf_file: Path,
    genome_url: str,
    gtf_url: str,
    logger: logging.Logger,
) -> None:
    '''Download genome reference given the urls to the fasta and gtf files'''   

    logger.info("Downloading genome reference files")
    helpers.download_file(genome_url, genome_file, logger)
    helpers.download_file(gtf_url, gtf_file, logger)


def prefetch_one_sra(srr: str, paths: PipelinePaths, logger: logging.Logger) -> None:
    """Download one SRA file with prefetch."""
    logger.info("Prefetching %s to %s", srr, paths.sra_dir)
    helpers.run_command(
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
    paths: PipelinePaths,
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
    paths: PipelinePaths,
    threads: int,
    logger: logging.Logger,
) -> None:
    '''Download SRA files from the cloud using srr accession number, then convert to FASTA.'''

    logger.info("Dumping %s to %s", srr, paths.dumped_dir / library.name)

    ## Prefetch SRA files (fasterq-dump did not work for me otherwise)
    logger.debug("Running prefetch for %s", srr)
    helpers.run_command(
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

    ## Dump FASTA files from SRA files
    logger.debug("Running fasterq-dump for %s", srr)
    helpers.run_command(
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
    helpers.run_command(["pigz", str(library.read1_fasta), "-p", str(threads)], logger)
    helpers.run_command(["pigz", str(library.read2_fasta), "-p", str(threads)], logger)


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
            batch.write(
                f"{library.name}\t{library.read1_fasta}.gz\t{library.read2_fasta}.gz\n"
            )

    logger.info("Multiplexing FASTA files with splitcode")
    helpers.run_command(
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

    helpers.run_command(
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
    paths: PipelinePaths,
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

    helpers.run_command(
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

def extractRandOPolyT(
    paths:PipelinePaths, 
    threads: int, 
    logger:logging.Logger
) -> None:
    '''Generate FASTQ file of randO reads with splitcode'''

    logger.info("Extracting RandO reads")

    with open(paths.randOpolyT_keep_file, "w") as keep_file:
        keep_file.write(f"r1_R {str(paths.randO_files[0]).split('_0')[0]}\n")
        keep_file.write(f"r1_T {str(paths.polyT_files[0]).split('_0')[0]}")

    helpers.run_command([
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

def pseudoalign10X(
    paths: PipelinePaths,
    fastq_files: list[Path],
    kb_out_dir: Path,
    tech: str,
    threads: int,
    logger: logging.Logger,
) -> None:
    '''Psuedoalign multiplexed files to reference. Build index if needed'''
    
    logger.info("Building kallisto index")
    helpers.run_command(
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
    helpers.run_command(
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

def pseudoalign10XHashtags(
    paths: PipelinePaths,
    fastq_files: list[Path],
    kb_out_dir: Path,
    tech: str,
    threads: int,
    logger: logging.Logger,
) -> None:
    '''Psuedoalign multiplexed files to reference. Build index if needed'''
    
    logger.info("Building kallisto index for 10X Hashtags")
    helpers.run_command(
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

    logger.info("Pseudoaligning 10X Hashtag multiplexed reads to genome index")
    helpers.run_command(
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
            str(fastq_files[1]),
        ],
        logger,
    )

def pseudoalignParse(
    paths: PipelinePaths,
    fastq_files: list[Path],
    kb_out_dir: Path,
    tech: str,
    threads: int,
    logger: logging.Logger,
) -> None:
    '''Pseudoalign all parse-specific files (FASTQ with all parse reads, polyT FASTQ, randO FASTQ)'''

    logger.info("Building kallisto index")
    helpers.run_command(
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
    helpers.run_command(
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

def core_pipeline(settings: RunSettings, paths: PipelinePaths, config:AnalysisConfig, assay:str, logger) -> None:
    '''This is the core pipeline for both loading in FASTQ files from GEO'''
    
    logger.info(f"Starting {assay} pipeline")
    
    # Generate names for the dumped files for each SRR number
    libraries = helpers.build_libraries(config, paths)

    # Download reference
    ref_files_exist = (paths.genome_file.is_file() or paths.genome_file == None) and (paths.gtf_file == None or paths.gtf_file.is_file())
    if not ref_files_exist or settings.overwrite:
        get_reference(
            genome_file=paths.genome_file,
            gtf_file=paths.gtf_file,
            genome_url=config.genome_url,
            gtf_url=config.gtf_url,
            logger=logger,
        )
    else:
        logger.info("Reference genome files already exist. Skipping reference download")

    # Dump SRA files from cloud to zipped FASTA
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

    # Multiplex all dumped FASTA file into a single FASTQ
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