"""Write the splitcode remultiplexing batch file for one assay.

One row per library: ``<name>\t<first>\t<second>``. splitcode groups rows sharing
a name into one sublibrary.

Mirrors utils.multiplex_fastqs' batch writing: when a sublibrary barcode is used
(``bclen`` set, the Parse case) the barcode read must come first, so R2/R1 are
swapped; with ``bclen`` None (10x) the order is R1/R2.

``params.rows`` is a list of (name, r1, r2) tuples computed source-aware by the
Snakefile: accession-named dumped FASTQs for SRA/ERA, pre-existing Lib{i} files
for a local assay.

Run as a Snakemake `script:`.
"""
from pathlib import Path

p = snakemake.params        # noqa: F821 (injected by Snakemake)
batch_file = Path(snakemake.output.batch)  # noqa: F821

rows = list(p.rows)  # (name, r1, r2) per library, in batch (barcode-assignment) order
bclen = p.bclen      # None for 10x, LIB_BC_LEN for Parse
barcode_read_first = bclen is not None

batch_file.parent.mkdir(parents=True, exist_ok=True)
with open(batch_file, "w") as batch:
    for name, r1, r2 in rows:
        first, second = (r2, r1) if barcode_read_first else (r1, r2)
        batch.write(f"{name}\t{first}\t{second}\n")
