"""Filesystem path templates for the 10XvParse Snakemake workflow.

The single source of truth for every path the workflow reads or writes — edit a
naming convention here and the whole DAG follows. Each function returns a plain
string built from the ``Data/{analysis}/{assay}/...`` (per-dataset outputs) and
``Index/{species}/...`` (shared references) conventions, so the same call serves as:

  * a concrete path — ``parse_filtered("Analysis_5", "parse")``
  * a Snakemake rule pattern — ``parse_filtered("{analysis}", "{assay}")``

All paths are relative to the repository root, which is the Snakemake working
directory.
"""
from __future__ import annotations

import re

# --------------------------------------------------------------------------- #
# base directories
# --------------------------------------------------------------------------- #

def data_dir(analysis: str, assay: str) -> str:
    return f"Data/{analysis}/{assay}"


def index_dir(species: str) -> str:
    return f"Index/{species}"


def star_index_dir(species: str) -> str:
    return f"Index/{species}_STAR"


def configs_dir(analysis: str, assay: str) -> str:
    """Per-assay config working dir (generated Parse files, hashtags.tsv)."""
    return f"Configs/{analysis}/{assay}"


PARSE_INFO_DIR = "Configs/parse_info"


def tenx_whitelist(tech: str) -> str:
    return f"Configs/10x_info/{tech}_whitelist.txt"


def plots_dir(analysis: str) -> str:
    return f"Data/{analysis}/Plots"


# --------------------------------------------------------------------------- #
# reference genome + kb (kallisto) index, shared per species
# --------------------------------------------------------------------------- #

def genome_file(species: str) -> str:
    return f"{index_dir(species)}/ref.fa.gz"


def gtf_file(species: str) -> str:
    return f"{index_dir(species)}/ref.gtf.gz"


def index_file(species: str) -> str:
    return f"{index_dir(species)}/index.idx"


def t2g_file(species: str) -> str:
    return f"{index_dir(species)}/t2g.txt"


def cdna_file(species: str) -> str:
    return f"{index_dir(species)}/cdna.txt"


def nascent_file(species: str) -> str:
    return f"{index_dir(species)}/nascent.txt"


def cdna_fasta_file(species: str) -> str:
    return f"{index_dir(species)}/cdna.fasta"


def nascent_fasta_file(species: str) -> str:
    return f"{index_dir(species)}/nascent.fasta"


def kb_index_files(species: str) -> list[str]:
    """All six files produced by ``kb ref --workflow nac`` for a species."""
    return [
        index_file(species),
        t2g_file(species),
        cdna_file(species),
        nascent_file(species),
        cdna_fasta_file(species),
        nascent_fasta_file(species),
    ]


# --------------------------------------------------------------------------- #
# hashtags (kite) index — per-analysis/assay, built from a committed hashtags.tsv
# --------------------------------------------------------------------------- #

def hashtags_index_dir(analysis: str, assay: str) -> str:
    return f"Index/{analysis}/{assay}"


def hashtags_index_file(analysis: str, assay: str) -> str:
    return f"{hashtags_index_dir(analysis, assay)}/index.idx"


def hashtags_t2g_file(analysis: str, assay: str) -> str:
    return f"{hashtags_index_dir(analysis, assay)}/t2g.txt"


def hashtags_cdna_fasta_file(analysis: str, assay: str) -> str:
    """The -f1 output of ``kb ref --workflow kite`` (the feature-barcode FASTA)."""
    return f"{hashtags_index_dir(analysis, assay)}/cdna.fasta"


def hashtags_tsv(analysis: str, assay: str) -> str:
    """Committed feature-barcode table (name<TAB>barcode) that seeds the kite index."""
    return f"{configs_dir(analysis, assay)}/hashtags.tsv"


# --------------------------------------------------------------------------- #
# FASTQ paths (each list is exactly [read0, read1])
# --------------------------------------------------------------------------- #

def sra_dir(analysis: str, assay: str) -> str:
    return f"{data_dir(analysis, assay)}/SRA"


def sra_file(analysis: str, assay: str, srr: str) -> str:
    return f"{sra_dir(analysis, assay)}/{srr}/{srr}.sra"


def tmp_dir(analysis: str, assay: str) -> str:
    return f"{data_dir(analysis, assay)}/tmp"


def dumped_dir(analysis: str, assay: str) -> str:
    return f"{data_dir(analysis, assay)}/FASTA/Dumped"


def dumped_fastq(analysis: str, assay: str, srr: str, read: str) -> str:
    """A per-accession dumped FASTQ. ``read`` is the logical read label 'R1'/'R2'
    (which fasterq-dump suffix each maps to is set by the dump rule's params)."""
    return f"{dumped_dir(analysis, assay)}/{srr}_{read}.fasta.gz"


def batch_file(analysis: str, assay: str) -> str:
    return f"{data_dir(analysis, assay)}/splitcode_configs/batch.txt"


def processed_dir(analysis: str, assay: str) -> str:
    return f"{data_dir(analysis, assay)}/FASTA/Processed"


def sampled_dir(analysis: str, assay: str, stag: str = "") -> str:
    """Subsampled-FASTQ dir. ``stag`` is a *tag token* — ``""`` or ``"_<tag>"``
    (leading underscore included) — so ``Sampled{stag}`` is ``Sampled`` untagged
    or ``Sampled_standard`` when an analysis subsamples to more than one depth
    (Analysis 2 standard/mini). Every sampled output of a comparison group carries
    that group's token uniformly, so one wildcard threads the whole chain; the
    token form (vs a raw tag) is what lets that wildcard also match the empty case."""
    return f"{data_dir(analysis, assay)}/FASTA/Sampled{stag}"


def tenx_multiplexed(analysis: str, assay: str) -> list[str]:
    d = processed_dir(analysis, assay)
    return [f"{d}/{assay}_{i}.fastq.gz" for i in range(2)]


def hashtags_trimmed_r2(analysis: str, assay: str) -> str:
    """R2 of a hashtag library trimmed to the feature-barcode length (seqtk trimfq).
    Named ``*_trimmed.fastq.gz`` next to the multiplexed R2 it is derived from."""
    return f"{processed_dir(analysis, assay)}/{assay}_1_trimmed.fastq.gz"


def parse_multiplexed(analysis: str, assay: str) -> list[str]:
    """Parse multiplexes under a generic name, not the assay name."""
    d = processed_dir(analysis, assay)
    return [f"{d}/multiplexed_{i}.fastq.gz" for i in range(2)]


def parse_filtered(analysis: str, assay: str) -> list[str]:
    d = processed_dir(analysis, assay)
    return [f"{d}/{assay}_{i}.fastq.gz" for i in range(2)]


def parse_split(analysis: str, assay: str, kind: str) -> list[str]:
    """Barcode-split Parse FASTQs. ``kind`` is 'polyT' or 'randO'."""
    d = processed_dir(analysis, assay)
    return [f"{d}/{kind}_{i}.fastq.gz" for i in range(2)]


def sampled_fastqs(analysis: str, assay: str, stag: str = "", kind: str = "") -> list[str]:
    """Both subsampled FASTQs. ``kind`` empty -> ``{assay}_{i}``; 'polyT'/'randO'
    -> that prefix. ``stag`` is the depth token (``""`` or ``"_<tag>"``)."""
    return [sampled_fastq(analysis, assay, i, stag, kind) for i in range(2)]


def sampled_fastq(analysis: str, assay: str, read, stag: str = "", kind: str = "") -> str:
    """One subsampled FASTQ (``read`` = 0 or 1). Each read is subsampled as its own
    job, so rules target this single-file form. ``stag`` is the depth token."""
    d = sampled_dir(analysis, assay, stag)
    prefix = kind or assay
    return f"{d}/{prefix}_{read}.fastq.gz"


# --------------------------------------------------------------------------- #
# kb count output directories + h5ad markers
# --------------------------------------------------------------------------- #

def kb_dir(analysis: str, assay: str) -> str:
    return f"{data_dir(analysis, assay)}/kb_python"


def kb_out(analysis: str, assay: str, leaf: str) -> str:
    """A kb count output dir. ``leaf`` is e.g. '10x_out', 'parse_out',
    'sampled_parse_out', 'sampled_polyT_out'."""
    return f"{kb_dir(analysis, assay)}/{leaf}"


# The single grammar for a *sampled* kb-count leaf, so building it (targets) and
# reading it back (the kb-count input functions) stay in one place.
SAMPLED_LEAF_RE = r"sampled_(10x|parse|polyT|randO)(_[A-Za-z0-9]+)?_out"


def sampled_leaf(sub: str, stag: str = "") -> str:
    """kb output-dir leaf for a sampled count. ``sub`` is '10x'/'parse'/'polyT'/'randO';
    ``stag`` is the depth token ('' or '_<tag>'). e.g. ('parse', '_standard') ->
    'sampled_parse_standard_out'."""
    return f"sampled_{sub}{stag}_out"


def parse_sampled_leaf(leaf: str):
    """Inverse of ``sampled_leaf``: 'sampled_polyT_standard_out' -> ('polyT', '_standard').
    Returns ``None`` for a non-sampled (full-depth) leaf like '10x_out'/'parse_out'."""
    m = re.fullmatch(SAMPLED_LEAF_RE, leaf)
    return (m.group(1), m.group(2) or "") if m else None


def kb_h5ad(analysis: str, assay: str, leaf: str) -> str:
    """The adata.h5ad the nac/kite workflow writes inside a kb count dir; used
    as the rule output marker for a pseudoalignment step."""
    return f"{kb_out(analysis, assay, leaf)}/counts_unfiltered/adata.h5ad"


def kb_h5ad_modified(analysis: str, assay: str, leaf: str) -> str:
    """The additional adata.h5ad kb count writes when given -r (barcode replace):
    the randO bc1 barcodes are merged into their polyT counterparts. Produced by
    every Parse kb count (all use -r), never by the 10x runs."""
    return f"{kb_out(analysis, assay, leaf)}/counts_unfiltered_modified/adata.h5ad"


# --------------------------------------------------------------------------- #
# Parse config files generated by parse_config.generate_parse_configs
# --------------------------------------------------------------------------- #

# Files written by generate_parse_configs (before any splitcode run). Parse always
# runs with sublibrary barcodes, so lib_bc.txt / sublibraries.txt are always written.
PARSE_GENERATED = [
    "r1_R.txt",
    "r1_T.txt",
    "bcs_to_wells.txt",
    "onlist.txt",
    "replace.txt",
    "star_bc1.txt",
    "star_bc2.txt",
    "star_bc3.txt",
    "lib_bc.txt",
    "sublibraries.txt",
]


def parse_generated_files(analysis: str, assay: str) -> list[str]:
    d = configs_dir(analysis, assay)
    return [f"{d}/{name}" for name in PARSE_GENERATED]


def parse_onlist(analysis: str, assay: str) -> str:
    return f"{configs_dir(analysis, assay)}/onlist.txt"


def parse_replace(analysis: str, assay: str) -> str:
    return f"{configs_dir(analysis, assay)}/replace.txt"


def parse_x_string(analysis: str, assay: str) -> str:
    """kb-python x_string written by the parse_configs step; read by kb count -x."""
    return f"{configs_dir(analysis, assay)}/x_string.txt"


def parse_r1_R(analysis: str, assay: str) -> str:
    return f"{configs_dir(analysis, assay)}/r1_R.txt"


def parse_r1_T(analysis: str, assay: str) -> str:
    return f"{configs_dir(analysis, assay)}/r1_T.txt"


# Files written at run time by the filter / split splitcode steps.
def parse_splitcode_config(analysis: str, assay: str) -> str:
    return f"{configs_dir(analysis, assay)}/config_RT_parse.txt"


def parse_keep_file(analysis: str, assay: str) -> str:
    return f"{configs_dir(analysis, assay)}/parse_keep.txt"


def parse_randOpolyT_keep_file(analysis: str, assay: str) -> str:
    return f"{configs_dir(analysis, assay)}/randOpolyT_keep.txt"


# --------------------------------------------------------------------------- #
# cross-assay subsample read-count coupling
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# gene-body coverage: STAR index + STARsolo BAMs + reference BED + RSeQC plot
# --------------------------------------------------------------------------- #

def genome_fasta_unzipped(species: str) -> str:
    """Uncompressed reference FASTA (STAR genomeGenerate / gffread need it plain)."""
    return f"{index_dir(species)}/ref.fa"


def gtf_unzipped(species: str) -> str:
    return f"{index_dir(species)}/ref.gtf"


def star_index_marker(species: str) -> str:
    """A file STAR always writes into the genome dir — the STAR index build marker."""
    return f"{star_index_dir(species)}/genomeParameters.txt"


def bed_file(species: str) -> str:
    return f"{index_dir(species)}/ref.bed"


def hk_genes_file(species: str) -> str:
    return f"{index_dir(species)}/hk_genes.txt"


def hk_bed_file(species: str) -> str:
    return f"{index_dir(species)}/hk_ref.bed"


def star_dir(analysis: str, assay: str) -> str:
    return f"{data_dir(analysis, assay)}/STARsolo"


def star_prefix(analysis: str, assay: str, label: str, stag: str = "") -> str:
    """STAR ``--outFileNamePrefix``. ``label`` is '10x'/'parse'/'polyT'/'randO';
    ``stag`` tags a multiply-subsampled 10x (empty otherwise), matching the sampled
    outputs. STAR appends 'Aligned.sortedByCoord.out.bam' etc. to this."""
    return f"{star_dir(analysis, assay)}/{label}{stag}_"


def star_bam(analysis: str, assay: str, label: str, stag: str = "") -> str:
    return f"{star_prefix(analysis, assay, label, stag)}Aligned.sortedByCoord.out.bam"


# STAR CB/UMI positions + whitelist order for a Parse assay, derived from its
# x_string by scripts/parse_star_params.py.
def parse_star_positions(analysis: str, assay: str) -> str:
    return f"{configs_dir(analysis, assay)}/star_cb_positions.txt"


def parse_star_umi(analysis: str, assay: str) -> str:
    return f"{configs_dir(analysis, assay)}/star_umi_position.txt"


def parse_star_whitelists(analysis: str, assay: str) -> str:
    return f"{configs_dir(analysis, assay)}/star_whitelists.txt"


def genebody_prefix(analysis: str, tag: str = "") -> str:
    """geneBody_coverage.py ``-o`` prefix for a comparison group (tag '' -> the
    plots dir itself, giving '.geneBodyCoverage.*')."""
    return f"{plots_dir(analysis)}/{tag}"


def genebody_marker(analysis: str, tag: str = "") -> str:
    return f"{genebody_prefix(analysis, tag)}.geneBodyCoverage.txt"


def read_counts_file(analysis: str) -> str:
    """One file per analysis holding a ``<fastq> <read_count>`` line for each
    first-read processed FASTQ. Subsample jobs compute their group's minimum on
    the fly from the relevant lines. The file is group-agnostic, so changing the
    comparison groups doesn't invalidate it."""
    return f"Data/{analysis}/read_counts.txt"
