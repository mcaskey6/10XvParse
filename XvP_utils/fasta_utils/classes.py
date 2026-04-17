from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
import logging
import yaml

def make_dir(directory: Path, logger: logging.Logger) -> None:
    '''Wrapper for Path.mkdir with logging'''
    try:
        directory.mkdir(parents=True, exist_ok=True)
        logger.debug("Directory '%s' is ready.", directory)
    except Exception as e:
        logger.error("Failed to create directory '%s': %s", directory, e)
        raise

@dataclass(frozen=True)
class RunSettings:
    '''Parameters and variables specifying the pipeline to run'''
    root_dir: Path              # The project directory
    config_name: str            # Filename with the config file for the analysis to be run
    overwrite: bool = False     # Whether or not to re-download FASTA files and rerun splitcode
    subsample_num: int = 10000  # Number of reads to subsample for comparison
    threads: int = 8            # Number of threads to pass to commands (splitcode, kb_python, etc.)
    max_workers: int = 2        # Number of workers to use for parallel steps in the pipeline.

@dataclass(frozen=True)
class AnalysisConfig:
    '''Object to store information from the analysis YAML file'''
    name: str
    sra: list[str]
    r1_num: int
    r2_num: int
    genome_url: str
    gtf_url: str
    species: str
    technology: str

    @classmethod
    def from_yaml(cls, config_file: Path, assay: str) -> "AnalysisConfig":
        '''Reads yaml file to initialize object'''
        with open(config_file, "r") as f:
            raw = yaml.safe_load(f)

        return cls(
            name=raw["name"],
            sra=raw["SRA"][assay],
            r1_num=raw["read_num"][assay]["R1"],
            r2_num=raw["read_num"][assay]["R2"],
            genome_url=raw["reference"]["fasta"],
            gtf_url=raw["reference"]["gtf"],
            species=raw["reference"]["species"],
            technology=raw["tech"][assay],
        )

@dataclass(frozen=True)
class PipelinePaths:
    # Config Paths
    config_dir: Path
    config_file: Path

    ## Data Paths (FASTA/FASTQ/SRA)
    outdir: Path
    fasta_dir: Path
    dumped_dir: Path
    processed_dir: Path
    sampled_dir: Path
    sra_dir: Path
    tmp_dir: Path

    # Log Paths
    log_dir: Path
    log_file: Path

    # splitcode Paths
    split_dir: Path
    batch_file: Path

    # Parse Configs (Parse only)
    parse_config: Path
    polyT_barcodes: list[str]
    randO_barcodes: list[str]
    parse_keep_file: Path
    randOpolyT_keep_file: Path
    kb_replace_config: Path
    kb_onlist: Path

    # FASTQ Paths
    multiplexed_files: list[Path]
    sampled_files: list[Path]
    ## Parse only
    filtered_files: list[Path]
    polyT_files: list[Path]
    randO_files: list[Path]
    sampled_polyT_files: list[Path]
    sampled_randO_files: list[Path]

    # kb ref Paths
    index_dir: Path
    index_file: Path
    t2g_file: Path
    cdna_file: Path
    nascent_file: Path
    cdna_fasta_file: Path
    nascent_fasta_file: Path
    genome_file: Path
    gtf_file: Path

    # kb count Paths
    kb_dir: Path
    kb_all_dir: Path
    kb_sub_dir: Path
    ## Parse only
    kb_sub_polyT_dir: Path
    kb_sub_randO_dir: Path

    @classmethod
    def build_base(cls, settings: RunSettings, config: AnalysisConfig, assay:str) -> "PipelinePaths":
        '''Initialize paths required for all pipelines'''
        config_dir = settings.root_dir / "Configs"
        config_file = config_dir / settings.config_name

        outdir = settings.root_dir / "Data" / config.name / assay
        fasta_dir = outdir / "FASTA"
        dumped_dir = fasta_dir / "Dumped"
        processed_dir = fasta_dir / "Processed"
        sampled_dir = fasta_dir / "Sampled"
        sra_dir = outdir / "SRA"
        tmp_dir = outdir / "tmp"
        split_dir = outdir / "splitcode_configs"

        log_dir = settings.root_dir / "Logs" / config.name
        log_file = log_dir / f"{assay}_log.txt"
        
        batch_file = split_dir / "batch.txt"

        kb_dir = outdir / "kb_python"
        kb_all_dir = kb_dir / "all_out"
        kb_sub_dir = kb_dir / "sampled_out"

        multiplexed_files = [processed_dir / f"{assay}_{i}.fastq.gz" for i in range(2)]
        sampled_files = [sampled_dir / f"{assay}_{i}.fastq.gz" for i in range(2)]

        return cls(
            config_dir=config_dir,
            config_file=config_file,
            outdir=outdir,
            fasta_dir=fasta_dir,
            dumped_dir=dumped_dir,
            processed_dir=processed_dir,
            sampled_dir=sampled_dir,
            sra_dir=sra_dir,
            tmp_dir=tmp_dir,
            log_dir=log_dir,
            split_dir=split_dir,
            log_file=log_file,
            batch_file=batch_file,
            parse_config=None,
            polyT_barcodes=None,
            randO_barcodes=None,
            parse_keep_file=None,
            randOpolyT_keep_file=None,
            kb_onlist=None,
            kb_replace_config=None,
            multiplexed_files = multiplexed_files,
            filtered_files=None,
            sampled_files=sampled_files,
            polyT_files=None,
            randO_files=None,
            sampled_polyT_files=None,
            sampled_randO_files=None,
            index_dir=None,
            index_file=None,
            t2g_file=None,
            cdna_file=None,
            nascent_file=None,
            nascent_fasta_file=None,
            cdna_fasta_file=None,
            genome_file=None,
            gtf_file=None,
            kb_dir=kb_dir,
            kb_all_dir=kb_all_dir,
            kb_sub_dir=kb_sub_dir,
            kb_sub_polyT_dir=None,
            kb_sub_randO_dir=None
        )
    
    @classmethod
    def build_10x(cls, settings: RunSettings, config: AnalysisConfig, assay: str) -> "PipelinePaths":
        "Initialize paths for 10X pipeline. This is also the base for the parse pipeline."
        base = cls.build_base(settings, config, assay)

        index_dir = settings.root_dir / "Index" / config.species
        index_file = index_dir / "index.idx"
        t2g_file = index_dir / "t2g.txt"
        cdna_file = index_dir / "cdna.txt"
        nascent_file = index_dir / "nascent.txt"
        nascent_fasta_file = index_dir / "nascent.fasta"
        cdna_fasta_file = index_dir / "cdna.fasta"
        genome_file = index_dir / "ref.fa.gz"
        gtf_file = index_dir / "ref.gtf"

        return replace(
            base,
            index_dir=index_dir,
            index_file=index_file,
            t2g_file=t2g_file,
            cdna_file=cdna_file,
            nascent_file=nascent_file,
            cdna_fasta_file=cdna_fasta_file,
            nascent_fasta_file=nascent_fasta_file,
            genome_file=genome_file,
            gtf_file=gtf_file,
        )
    
    @classmethod
    def build_10xhashtags(cls, settings: RunSettings, config: AnalysisConfig, assay: str) -> "PipelinePaths":
        "Initialize paths for 10X Hashtag pipeline."

        base = cls.build_base(settings, config, assay)

        index_dir = settings.root_dir / "Index" / config.name / assay
        index_file = index_dir / "index.idx"
        t2g_file = index_dir / "t2g.txt"
        cdna_file = index_dir / "cdna.fasta"
        
        genome_file = base.config_dir / config.name / assay / "hashtags.tsv"


        return replace(
            base,
            index_dir=index_dir,
            index_file=index_file,
            t2g_file=t2g_file,
            cdna_file=cdna_file,
            genome_file=genome_file
        )


    @classmethod
    def build_parse(cls, settings: RunSettings, config: AnalysisConfig, assay: str) -> "PipelinePaths":
        base = cls.build_10x(settings, config, assay)

        configs_dir = base.config_dir / config.name / assay
        parse_config = configs_dir / "config_RT_parse.txt"
        polyT_barcodes = configs_dir / "r1_T.txt"
        randO_barcodes = configs_dir / "r1_R.txt"
        parse_keep_file = configs_dir / "parse_keep.txt"
        randOpolyT_keep_file = configs_dir / "randOpolyT_keep.txt"
        kb_onlist = configs_dir / "onlist.txt"
        kb_replace_config = configs_dir / "replace.txt"

        multiplexed_files = [base.processed_dir / f"multiplexed_{i}.fastq.gz" for i in range(2)]
        filtered_files = [base.processed_dir / f"{assay}_{i}.fastq.gz" for i in range(2)]
        polyT_files = [base.processed_dir / f"polyT_{i}.fastq.gz" for i in range(2)]
        randO_files = [base.processed_dir / f"randO_{i}.fastq.gz" for i in range(2)]

        sampled_polyT_files = [base.sampled_dir / f"polyT_{i}.fastq.gz" for i in range(2)]
        sampled_randO_files = [base.sampled_dir / f"randO_{i}.fastq.gz" for i in range(2)]

        kb_sub_polyT_dir = base.kb_dir / "sampled_polyT_out"
        kb_sub_randO_dir = base.kb_dir / "sampled_randO_out"

        return replace(
            base,
            parse_config=parse_config,
            polyT_barcodes=polyT_barcodes,
            randO_barcodes=randO_barcodes,
            parse_keep_file=parse_keep_file,
            randOpolyT_keep_file=randOpolyT_keep_file,
            multiplexed_files=multiplexed_files,
            filtered_files=filtered_files,
            sampled_polyT_files=sampled_polyT_files,
            sampled_randO_files=sampled_randO_files,
            kb_onlist=kb_onlist,
            kb_replace_config=kb_replace_config,
            polyT_files=polyT_files,
            randO_files=randO_files,
            kb_sub_polyT_dir=kb_sub_polyT_dir,
            kb_sub_randO_dir=kb_sub_randO_dir
        )
    
    def directories(self):
        '''List of directories that need to be created and are not assumed to exist before the pipeline is run'''
        return[
            self.outdir,
            self.fasta_dir,
            self.dumped_dir,
            self.processed_dir,
            self.sampled_dir,
            self.sra_dir,
            self.tmp_dir,
            self.log_dir,
            self.split_dir,
            self.index_dir,
            self.kb_dir,
        ]
    
    def ensure_dirs(self, logger):
        '''Ensures directories exist'''
        for dir in self.directories():
            if dir:
                make_dir(dir, logger)
    
@dataclass(frozen=True)
class LibraryFiles:
    '''An object to hold paired-end read files'''
    name: str
    read1_fasta: Path
    read2_fasta: Path

    @property
    def output_files(self) -> list[Path]:
        return [self.read1_fasta, self.read2_fasta]

    @property
    def gz_files(self) -> list[Path]:
        return [
            Path(f"{self.read1_fasta}.gz"),
            Path(f"{self.read2_fasta}.gz"),
        ]