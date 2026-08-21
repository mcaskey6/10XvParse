"""Generate all Parse config files for one assay.

Reuses XvP_utils.preprocessing.parse_config.generate_parse_configs (the barcode
CSV parsing + onlist/replace/star whitelist logic) and additionally writes the
three splitcode helper files that the old filter_parse_fastqs / extract_rando_polyt
functions used to write inline — so the downstream splitcode rules can be pure
`shell:` steps. Also writes the kb-python x_string to a file for the parse kb
count rule to read.

Run as a Snakemake `script:` — the injected ``snakemake`` object supplies params
and outputs.
"""
from pathlib import Path

from XvP_utils.preprocessing.parse_config import (
    generate_parse_configs,
    sublibrary_barcodes,
    x_string_to_bc1_location,
)

p = snakemake.params  # noqa: F821 (injected by Snakemake)
o = snakemake.output  # noqa: F821

out_dir = Path(p.out_dir)
wells = list(p.wells) if p.wells else None

# Parse always registers the sublibrary "4th round" barcode: this shifts the
# returned x_string, adds the leading onlist column, and writes lib_bc.txt /
# sublibraries.txt.
sublibraries = list(p.sublibraries)
if not sublibraries:
    raise ValueError(
        f"No sublibraries found for Parse assay in {out_dir}: expected input FASTQs "
        "were not present (e.g. a local-source assay with no Lib*_*.fastq.gz in its "
        "Dumped dir, or an empty SRA/ERA accession list). Parse always remultiplexes "
        "with a sublibrary barcode, so at least one library is required."
    )
x_string = generate_parse_configs(
    kit_name=p.kit,
    parse_info_dir=Path(p.parse_info_dir),
    output_dir=out_dir,
    wells=wells,
    logger=None,
    lib_barcodes=sublibrary_barcodes(sublibraries),
    sublibraries=sublibraries,
)

# splitcode filtering config (was written inside filter_parse_fastqs). Paths are
# relative to the repo root, which is the Snakemake working directory.
bc1_location = x_string_to_bc1_location(x_string)
with open(o.splitcode_config, "w") as c:
    c.write("tags\tdistances\tids\tgroups\tminFindsG\tlocations\n")
    c.write(f"{out_dir/'r1_R.txt'}\t1\tr1_R\tround1\t1\t{bc1_location}\n")
    c.write(f"{out_dir/'r1_T.txt'}\t1\tr1_T\tround1\t1\t{bc1_location}\n")

# keep-group file for the filter step (round1 -> filtered output prefix)
Path(o.keep).write_text(f"round1 {p.filtered_prefix}")

# keep file for the randO/polyT split step
Path(o.rando_keep).write_text(
    f"r1_R {p.randO_prefix}\nr1_T {p.polyT_prefix}"
)

# x_string consumed by the parse kb count rules (kb count -x)
Path(o.x_string).write_text(x_string)
