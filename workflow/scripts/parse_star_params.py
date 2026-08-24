"""Derive STARsolo CB_UMI_Complex parameters for one Parse assay from its x_string.

Writes three small files the starsolo_parse rule reads:
  * star_cb_positions.txt  — space-separated --soloCBposition anchors
  * star_umi_position.txt  — the --soloUMIposition anchor
  * star_whitelists.txt    — the --soloCBwhitelist files, in x_string order
                             (star_bc3/2/1, plus lib_bc when the reads carry the
                             4th-round sublibrary barcode).

Reuses parse_config.x_string_to_star_params (the validated kb->STAR position
conversion). Run as a Snakemake `script:`.
"""
from pathlib import Path

from parse_config import x_string_to_star_params  # sibling module in workflow/scripts/

o = snakemake.output          # noqa: F821 (injected by Snakemake)
configs_dir = Path(snakemake.params.configs_dir)  # noqa: F821

x_string = Path(snakemake.input.x_string).read_text().strip()  # noqa: F821
cb_positions, umi_position = x_string_to_star_params(x_string)

# One whitelist per CB position, in x_string order; the sublibrary barcode comes
# first when the reads were remultiplexed with one (an extra CB position).
whitelists = [str(configs_dir / f"star_bc{i}.txt") for i in (3, 2, 1)]
if len(cb_positions) == len(whitelists) + 1:
    whitelists.insert(0, str(configs_dir / "lib_bc.txt"))

Path(o.positions).write_text(" ".join(cb_positions))
Path(o.umi).write_text(umi_position)
Path(o.whitelists).write_text(" ".join(whitelists))
