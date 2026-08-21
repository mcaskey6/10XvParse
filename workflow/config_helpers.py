"""Config parsing helpers for the 10XvParse Snakemake workflow.

These are the load-bearing bits salvaged from ``AnalysisConfig`` in
``XvP_utils/preprocessing/classes.py``, rewritten to operate on the plain dicts
Snakemake's ``config`` mechanism provides instead of frozen dataclasses.

The behaviour here must match ``AnalysisConfig.from_yaml`` exactly — in
particular ``flatten_accessions`` preserves sublibrary order, which determines
which splitcode remultiplexing barcode each sublibrary receives (reordering
sublibraries in the YAML would silently invalidate previously generated data).

Reference URLs are NOT loaded here: they live in ``config/indexes.yaml`` and are
read as plain ``config["indexes"][species]`` dict access.
"""
from __future__ import annotations

import re

# Length of the synthetic sublibrary "4th round" barcode splitcode prepends when
# remultiplexing Parse libraries. Mirrors parse_config.LIB_BC_LEN; inlined so the
# Snakefile need not import XvP_utils at DAG-build time.
LIB_BC_LEN = 4

# Parse kit tech strings look like "<kit>_v<chem>" (e.g. WT_v2, WT_mini_v3);
# 10XV3/10XV4 do not match. Kept in sync with parse_config._PARSE_PATTERN, but
# inlined here so importing this module (during Snakefile parsing / DAG build)
# does not require the XvP_utils package or its runtime dependencies.
_PARSE_KIT_PATTERN = re.compile(r"^.+_v\d+$")


def is_parse_kit(tech_str: str) -> bool:
    return bool(_PARSE_KIT_PATTERN.match(tech_str))


# --------------------------------------------------------------------------- #
# accessions + sublibraries
# --------------------------------------------------------------------------- #

def flatten_accessions(entry) -> tuple[list[str], list[str]]:
    """Accept a flat accession list or a ``{sublibrary: [accessions]}`` mapping.

    Returns ``(accessions, sublibrary_labels)``, the labels list parallel to the
    accessions (empty for a flat list). Mapping key order is preserved and is
    load-bearing (see module docstring). Verbatim port of
    ``classes._flatten_accessions``.
    """
    if isinstance(entry, dict):
        accessions: list[str] = []
        sublibraries: list[str] = []
        for sublibrary, sublibrary_accessions in entry.items():
            for accession in sublibrary_accessions:
                accessions.append(accession)
                sublibraries.append(sublibrary)
        return accessions, sublibraries
    return list(entry or []), []


def sublibrary_labels(cfg: dict, assay: str) -> list[str]:
    """One splitcode batch-row label per library, in accession order — the labels
    the config nests accessions under, or the accession itself when listed flatly.
    For Parse these are the sublibraries whose barcodes are assigned in this order
    (load-bearing); for 10x they are just row names (no barcode is assigned).
    """
    _src, accs, subs = accessions_for(cfg, assay)
    return subs if subs else list(accs)


def accessions_for(cfg: dict, assay: str) -> tuple[str, list[str], list[str]]:
    """Return ``(source, accessions, sublibraries)`` for an assay.

    ``source`` is ``"SRA"``, ``"ERA"``, or ``"local"`` when neither is present
    (pre-downloaded files). Mirrors how ``AnalysisConfig.from_yaml`` reads the
    ``SRA``/``ERA`` blocks.
    """
    sra, sra_subs = flatten_accessions(cfg.get("SRA", {}).get(assay, []))
    if sra:
        return "SRA", sra, sra_subs
    era, era_subs = flatten_accessions(cfg.get("ERA", {}).get(assay, []))
    if era:
        return "ERA", era, era_subs
    return "local", [], []


# --------------------------------------------------------------------------- #
# read numbers (R1/R2 FASTQ suffixes, with per-sublibrary overrides)
# --------------------------------------------------------------------------- #

def read_num_block(cfg: dict, assay: str) -> dict:
    return cfg["read_num"][assay]


def default_read_nums(cfg: dict, assay: str) -> tuple[int, int]:
    block = read_num_block(cfg, assay)
    return block["R1"], block["R2"]


def read_num_overrides(cfg: dict, assay: str) -> dict[str, tuple[int, int]]:
    """Per-sublibrary ``{sublibrary: (R1, R2)}`` overrides, as in
    ``AnalysisConfig.from_yaml`` (any read_num key other than R1/R2 whose value
    is a dict)."""
    block = read_num_block(cfg, assay)
    return {
        sub: (spec["R1"], spec["R2"])
        for sub, spec in block.items()
        if sub not in ("R1", "R2") and isinstance(spec, dict)
    }


def read_nums_for(cfg: dict, assay: str, sublibrary: str) -> tuple[int, int]:
    """(R1, R2) for a sublibrary, honouring any per-sublibrary override."""
    overrides = read_num_overrides(cfg, assay)
    return overrides.get(sublibrary, default_read_nums(cfg, assay))


# --------------------------------------------------------------------------- #
# assorted per-assay accessors
# --------------------------------------------------------------------------- #

def technology(cfg: dict, assay: str) -> str:
    return cfg["tech"][assay]


def species(cfg: dict) -> str:
    return cfg["species"]


def wells(cfg: dict, assay: str) -> list[str]:
    return cfg.get("tech", {}).get("wells", {}).get(assay, [])


def r2_trim_length(cfg: dict, assay: str):
    return cfg.get("trim", {}).get(assay)


# --------------------------------------------------------------------------- #
# ENA FTP URLs (ERA source) — verbatim port of utils.era_ftp_url
# --------------------------------------------------------------------------- #

def era_ftp_url(err: str, read_num: int) -> str:
    """The ENA FTP URL for one run accession and read number (1 or 2).

    ENA nests runs under a length-dependent subdirectory scheme; mirrors
    ``utils.era_ftp_url`` exactly so the download rule reproduces the old paths.
    """
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


def assay_type(cfg: dict, assay: str) -> str:
    """Classify an assay into a pipeline chain: 'parse', 'hashtags', or 'tenx'.

    Parse kits are recognised by their ``<kit>_v<chem>`` tech string; hashtag
    feature-barcode assays by 'hashtag' in the assay name (they use a kite index
    from a committed hashtags.tsv rather than a genome).
    """
    if is_parse_kit(technology(cfg, assay)):
        return "parse"
    if "hashtag" in assay:
        return "hashtags"
    return "tenx"


def is_barnyard(cfg: dict) -> bool:
    return "_" in species(cfg)


# --------------------------------------------------------------------------- #
# comparison groups (moved out of the old Scripts/analysisN.py)
# --------------------------------------------------------------------------- #

def comparisons(cfg: dict) -> list[dict]:
    """The subsample comparison groups declared in the analysis config.

    Each entry is ``{"tag": str, "tenx": [assays], "parse": [assays]}``. The tag
    ("" for none) selects the subsample depth; every assay listed in a group is
    subsampled to that group's shared minimum read count.
    """
    return cfg.get("comparisons", [])


def group_by_tag(cfg: dict, tag: str) -> dict:
    """The comparison group with this tag ("" for the single untagged group)."""
    for c in comparisons(cfg):
        if c.get("tag", "") == tag:
            return c
    return {}


def stag_token(tag: str) -> str:
    """A comparison tag -> the path token used across its sampled outputs:
    ``""`` for the untagged group, ``"_<tag>"`` otherwise (Analysis 2)."""
    return f"_{tag}" if tag else ""
