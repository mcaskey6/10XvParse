"""Helpers for reading a per-analysis config (Config/<analysis>.yaml).

The Snakefile reads each analysis config as a plain dict; the functions here
interpret it — the read sources and sublibraries, per-assay read numbers, assay
type, wells, trim, and the subsample comparison groups.

IMPORTANT: sublibrary order is load-bearing. ``flatten_accessions`` preserves the
order sublibraries appear in the YAML, and that order determines which splitcode
remultiplexing barcode each sublibrary receives — so reordering sublibraries in a
config would silently change the barcode assignment of already-processed data.

Reference URLs are read separately from ``Config/indexes.yaml`` (as
``config["indexes"][species]``), not here.
"""
from __future__ import annotations

import re

# Length of the synthetic sublibrary "4th round" barcode splitcode prepends when
# remultiplexing Parse libraries (kept equal to parse_config.LIB_BC_LEN). Defined
# here too so this config module has no dependencies beyond the standard library.
LIB_BC_LEN = 4

# A Parse kit tech string looks like "<kit>_v<chem>" (e.g. WT_v2, WT_mini_v3);
# 10x strings (10XV3/10XV4) don't match. Distinguishes Parse from 10x assays.
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
    load-bearing (see module docstring).
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

    ``source`` is ``"SRA"`` or ``"ERA"`` (accession downloads; a flat list or a
    ``{sublibrary: [accessions]}`` mapping), else ``"local"`` (pre-downloaded reads
    listed in the ``local:`` block — read those via ``local_libraries``, not this).
    """
    for source in ("SRA", "ERA"):
        names, subs = flatten_accessions(cfg.get(source, {}).get(assay, []))
        if names:
            return source, names, subs
    return "local", [], []


def local_libraries(cfg: dict, assay: str) -> list[tuple[str, str, str]]:
    """``[(name, r1_file, r2_file)]`` for a local-source assay, from the ``local:``
    block: per assay a list of ``[R1, R2]`` filename pairs (R1/R2 in the assay's read
    order), one per sublibrary in barcode-assignment order, named ``Lib{i}``. The
    filenames are used as-is (relative to the assay's Dumped dir), so arbitrary
    Illumina names work without renaming."""
    entries = cfg.get("local", {}).get(assay, [])
    return [(f"Lib{i}", r1, r2) for i, (r1, r2) in enumerate(entries)]


# --------------------------------------------------------------------------- #
# read numbers (R1/R2 FASTQ suffixes, with per-sublibrary overrides)
# --------------------------------------------------------------------------- #

def read_num_block(cfg: dict, assay: str) -> dict:
    return cfg["read_num"][assay]


def default_read_nums(cfg: dict, assay: str) -> tuple[int, int]:
    block = read_num_block(cfg, assay)
    return block["R1"], block["R2"]


def read_num_overrides(cfg: dict, assay: str) -> dict[str, tuple[int, int]]:
    """Per-sublibrary ``{sublibrary: (R1, R2)}`` read-number overrides — any
    ``read_num`` key other than R1/R2 whose value is itself an {R1, R2} dict."""
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
# ENA FTP URLs (ERA source)
# --------------------------------------------------------------------------- #

def era_ftp_url(err: str, read_num: int) -> str:
    """The ENA FTP URL for one run accession and read number (1 or 2).

    ENA nests a run under a subdirectory scheme that depends on the accession's
    length; this reproduces that layout so the download rule fetches the right file.
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


def assay_bclen(cfg: dict, assay: str):
    """The splitcode ``--bclen`` for an assay: the sublibrary-barcode length for
    Parse (which always remultiplexes with a 4th-round barcode), else ``None`` (10x
    and hashtag libraries carry no sublibrary barcode)."""
    return LIB_BC_LEN if assay_type(cfg, assay) == "parse" else None


# --------------------------------------------------------------------------- #
# comparison groups (which assays are subsampled together, and to how many depths)
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
    """A comparison tag -> a path token: ``""`` for the untagged group,
    ``"_<tag>"`` otherwise. Used on an output path only via ``output_stag``."""
    return f"_{tag}" if tag else ""


def assay_multiplicity(cfg: dict) -> dict:
    """{assay: how many comparison groups it appears in}. >1 only when the same
    assay is subsampled at more than one depth (e.g. Analysis 2's shared 10x)."""
    counts: dict = {}
    for comp in comparisons(cfg):
        for a in comp.get("tenx", []) + comp.get("parse", []):
            counts[a] = counts.get(a, 0) + 1
    return counts


def output_stag(cfg: dict, assay: str, tag: str) -> str:
    """Depth token to put on a subsample OUTPUT path: the group tag only when this
    assay is subsampled in more than one group (else ``""``). So single-subsampled
    assays (all of Analysis 3–7) stay untagged and only a doubly-subsampled assay
    (Analysis 2's 10x) is tagged to keep its two depths from colliding."""
    return stag_token(tag) if assay_multiplicity(cfg).get(assay, 0) > 1 else ""


def subsample_group(cfg: dict, assay: str, out_stag: str) -> dict:
    """The comparison group a subsample job belongs to. A tagged output names its
    group by tag; an untagged output's assay belongs to exactly one group."""
    if out_stag:
        return group_by_tag(cfg, out_stag[1:])
    for comp in comparisons(cfg):
        if assay in comp.get("tenx", []) or assay in comp.get("parse", []):
            return comp
    return {}
