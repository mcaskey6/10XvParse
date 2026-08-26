"""Build the per-species gene-attributes cache the notebooks read.

Runs XvP_utils' ``query_ensembl`` (or ``query_ensembl_combined`` for a barnyard
species) once, with ``overwrite=True``, so the cache has a single producer in the
DAG. The combo and Comparisons notebooks then read it with ``overwrite=False``
instead of each rebuilding it. Run as a Snakemake ``script:`` in the active env
(XvP_utils is importable there); reads ref.gtf / cdna.fasta / t2g.txt from the
species index dir and writes Notebooks/gene_info/<species>/gene_attributes.csv.
"""
from pathlib import Path

from XvP_utils.processing import query_ensembl, query_ensembl_combined

p = snakemake.params  # noqa: F821 (injected by Snakemake)

# The Snakemake working directory is the repo root, so project_dir "." puts the
# cache at Notebooks/gene_info/<species>/gene_attributes.csv (the declared output).
fn = query_ensembl_combined if p.barnyard else query_ensembl
fn(Path("."), index_dir=Path(p.index_dir), species=p.species, overwrite=True)
