from __future__ import annotations

from dataclasses import dataclass, field, fields
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
    run_kb: bool = True         # Whether or not to run kb-python pseudoalignment
    threads: int = 8            # Number of threads to pass to commands (splitcode, kb_python, etc.)
    max_workers: int = 2        # Number of workers to use for parallel steps in the pipeline.


@dataclass(frozen=True)
class AnalysisConfig:
    '''Object to store information from the analysis YAML file'''
    # Required fields
    name: str
    r1_num: int
    r2_num: int
    species: str
    technology: str
    # Accession lists — one of sra or era will be populated
    sra: list[str] = field(default_factory=list)
    era: list[str] = field(default_factory=list)
    # Reference URLs — genome_url_2/gtf_url_2 set for barnyard (dual-species) configs
    genome_url: str = ""
    gtf_url: str = ""
    genome_url_2: str = ""
    gtf_url_2: str = ""
    # Optional R2 trim length (bp) — used for feature barcode libraries where long reads
    # cause downstream k-mers to collide with kite index entries from other barcodes
    r2_trim_length: int | None = None
    # Optional well subset for Parse kits — empty list means use all wells
    wells: list[str] = field(default_factory=list)

    @property
    def is_barnyard(self) -> bool:
        return bool(self.genome_url_2)

    @classmethod
    def from_yaml(cls, config_file: Path, assay: str) -> "AnalysisConfig":
        '''Reads yaml file to initialize object. Supports both SRA and ERA accessions,
        and both single-species and barnyard (dual-species) reference configs.'''
        with open(config_file, "r") as f:
            raw = yaml.safe_load(f)

        ref = raw["reference"]
        is_barnyard = "human_fasta" in ref

        return cls(
            name=raw["name"],
            sra=raw.get("SRA", {}).get(assay, []),
            era=raw.get("ERA", {}).get(assay, []),
            r1_num=raw["read_num"][assay]["R1"],
            r2_num=raw["read_num"][assay]["R2"],
            genome_url=ref.get("human_fasta" if is_barnyard else "fasta", ""),
            gtf_url=ref.get("human_gtf" if is_barnyard else "gtf", ""),
            genome_url_2=ref.get("mouse_fasta", ""),
            gtf_url_2=ref.get("mouse_gtf", ""),
            species=ref["species"],
            technology=raw["tech"][assay],
            r2_trim_length=raw.get("trim", {}).get(assay),
            wells=raw.get("tech", {}).get("wells", {}).get(assay, []),
        )


@dataclass(frozen=True)
class BasePaths:
    '''Paths shared by all pipeline variants.'''
    # Config paths
    config_dir: Path
    config_file: Path
    parse_info_dir: Path

    # Data paths
    outdir: Path
    fasta_dir: Path
    dumped_dir: Path
    processed_dir: Path
    sampled_dir: Path
    sra_dir: Path
    tmp_dir: Path

    # splitcode paths
    split_dir: Path
    batch_file: Path

    # kb paths
    kb_dir: Path
    kb_all_dir: Path
    kb_sub_dir: Path

    # FASTQ paths
    multiplexed_files: list[Path]
    sampled_files: list[Path]

    @classmethod
    def build(cls, settings: RunSettings, config: AnalysisConfig, assay: str) -> "BasePaths":
        config_dir = settings.root_dir / "Configs"
        config_file = config_dir / settings.config_name
        outdir = settings.root_dir / "Data" / config.name / assay
        fasta_dir = outdir / "FASTA"
        kb_dir = outdir / "kb_python"

        return cls(
            config_dir=config_dir,
            config_file=config_file,
            parse_info_dir=config_dir / "parse_info",
            outdir=outdir,
            fasta_dir=fasta_dir,
            dumped_dir=fasta_dir / "Dumped",
            processed_dir=fasta_dir / "Processed",
            sampled_dir=fasta_dir / "Sampled",
            sra_dir=outdir / "SRA",
            tmp_dir=outdir / "tmp",
            split_dir=outdir / "splitcode_configs",
            batch_file=outdir / "splitcode_configs" / "batch.txt",
            kb_dir=kb_dir,
            kb_all_dir=kb_dir / "all_out",
            kb_sub_dir=kb_dir / "sampled_out",
            multiplexed_files=[fasta_dir / "Processed" / f"{assay}_{i}.fastq.gz" for i in range(2)],
            sampled_files=[fasta_dir / "Sampled" / f"{assay}_{i}.fastq.gz" for i in range(2)],
        )

    def directories(self) -> list[Path]:
        '''Directories that must exist before the pipeline runs.'''
        return [
            self.outdir,
            self.fasta_dir,
            self.dumped_dir,
            self.processed_dir,
            self.sampled_dir,
            self.sra_dir,
            self.tmp_dir,
            self.split_dir,
            self.kb_dir,
        ]

    def ensure_dirs(self, logger: logging.Logger) -> None:
        for d in self.directories():
            make_dir(d, logger)


@dataclass(frozen=True)
class TenXPaths(BasePaths):
    '''Paths for the 10x Genomics pipeline. Adds genome reference and kallisto index paths.'''
    index_dir: Path
    index_file: Path
    t2g_file: Path
    cdna_file: Path
    nascent_file: Path
    cdna_fasta_file: Path
    nascent_fasta_file: Path
    genome_file: Path
    gtf_file: Path

    @classmethod
    def build(cls, settings: RunSettings, config: AnalysisConfig, assay: str) -> "TenXPaths":
        base = BasePaths.build(settings, config, assay)
        index_dir = settings.root_dir / "Index" / config.species
        base_fields = {f.name: getattr(base, f.name) for f in fields(base)}
        return cls(
            **base_fields,
            index_dir=index_dir,
            index_file=index_dir / "index.idx",
            t2g_file=index_dir / "t2g.txt",
            cdna_file=index_dir / "cdna.txt",
            nascent_file=index_dir / "nascent.txt",
            cdna_fasta_file=index_dir / "cdna.fasta",
            nascent_fasta_file=index_dir / "nascent.fasta",
            genome_file=index_dir / "ref.fa.gz",
            gtf_file=index_dir / "ref.gtf",
        )

    def directories(self) -> list[Path]:
        return super().directories() + [self.index_dir]


@dataclass(frozen=True)
class HashtagsPaths(BasePaths):
    '''Paths for the 10x Hashtag demultiplexing pipeline. Uses a feature TSV instead of a genome.'''
    index_dir: Path
    index_file: Path
    t2g_file: Path
    cdna_file: Path
    genome_file: Path   # points to the committed hashtags.tsv feature file

    @classmethod
    def build(cls, settings: RunSettings, config: AnalysisConfig, assay: str) -> "HashtagsPaths":
        base = BasePaths.build(settings, config, assay)
        index_dir = settings.root_dir / "Index" / config.name / assay
        hashtags_file = base.config_dir / config.name / assay / "hashtags.tsv"
        base_fields = {f.name: getattr(base, f.name) for f in fields(base)}
        return cls(
            **base_fields,
            index_dir=index_dir,
            index_file=index_dir / "index.idx",
            t2g_file=index_dir / "t2g.txt",
            cdna_file=index_dir / "cdna.fasta",
            genome_file=hashtags_file,
        )

    def directories(self) -> list[Path]:
        return super().directories() + [self.index_dir]


@dataclass(frozen=True)
class ParsePaths(TenXPaths):
    '''Paths for the Parse Biosciences pipeline. Extends TenXPaths with barcode-split FASTQ paths.'''
    # splitcode config paths
    parse_config: Path
    polyT_barcodes: Path
    randO_barcodes: Path
    parse_keep_file: Path
    randOpolyT_keep_file: Path
    kb_onlist: Path
    kb_replace_config: Path

    # barcode-split FASTQ paths
    filtered_files: list[Path]
    polyT_files: list[Path]
    randO_files: list[Path]
    sampled_polyT_files: list[Path]
    sampled_randO_files: list[Path]

    # kb count output dirs for barcode-split reads
    kb_sub_polyT_dir: Path
    kb_sub_randO_dir: Path

    @classmethod
    def build(cls, settings: RunSettings, config: AnalysisConfig, assay: str) -> "ParsePaths":
        ten_x = TenXPaths.build(settings, config, assay)
        configs_dir = ten_x.config_dir / config.name / assay
        processed_dir = ten_x.processed_dir
        sampled_dir = ten_x.sampled_dir

        ten_x_fields = {f.name: getattr(ten_x, f.name) for f in fields(ten_x)}
        # Parse multiplexes under a generic name, not the assay name
        ten_x_fields["multiplexed_files"] = [processed_dir / f"multiplexed_{i}.fastq.gz" for i in range(2)]

        return cls(
            **ten_x_fields,
            parse_config=configs_dir / "config_RT_parse.txt",
            polyT_barcodes=configs_dir / "r1_T.txt",
            randO_barcodes=configs_dir / "r1_R.txt",
            parse_keep_file=configs_dir / "parse_keep.txt",
            randOpolyT_keep_file=configs_dir / "randOpolyT_keep.txt",
            kb_onlist=configs_dir / "onlist.txt",
            kb_replace_config=configs_dir / "replace.txt",
            filtered_files=[processed_dir / f"{assay}_{i}.fastq.gz" for i in range(2)],
            polyT_files=[processed_dir / f"polyT_{i}.fastq.gz" for i in range(2)],
            randO_files=[processed_dir / f"randO_{i}.fastq.gz" for i in range(2)],
            sampled_polyT_files=[sampled_dir / f"polyT_{i}.fastq.gz" for i in range(2)],
            sampled_randO_files=[sampled_dir / f"randO_{i}.fastq.gz" for i in range(2)],
            kb_sub_polyT_dir=ten_x.kb_dir / "sampled_polyT_out",
            kb_sub_randO_dir=ten_x.kb_dir / "sampled_randO_out",
        )


@dataclass(frozen=True)
class LibraryFiles:
    '''An object to hold paired-end read files'''
    name: str
    read1_fasta: Path
    read2_fasta: Path
    is_gzipped: bool = False    # True when files are already gzipped (e.g. ERA downloads)

    @property
    def output_files(self) -> list[Path]:
        return [self.read1_fasta, self.read2_fasta]

    @property
    def gz_files(self) -> list[Path]:
        if self.is_gzipped:
            return self.output_files
        return [
            Path(f"{self.read1_fasta}.gz"),
            Path(f"{self.read2_fasta}.gz"),
        ]
