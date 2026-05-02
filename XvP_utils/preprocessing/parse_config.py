from __future__ import annotations

import csv
import logging
import re
from pathlib import Path

_PARSE_PATTERN = re.compile(r"^.+_v\d+$")


def is_parse_kit(tech_str: str) -> bool:
    return bool(_PARSE_PATTERN.match(tech_str))


def _parse_kit_name(kit_name: str) -> tuple[str, int]:
    kit, chem_str = kit_name.rsplit("_v", 1)
    return kit, int(chem_str)


def _load_kit_info(parse_info_dir: Path, kit: str, chem: int) -> dict[str, str]:
    kits_file = parse_info_dir / "kits_info.txt"
    with open(kits_file) as f:
        lines = f.readlines()

    header = lines[0].split()
    for line in lines[1:]:
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        if len(parts) < len(header):
            continue
        row = dict(zip(header, parts))
        if row["kit"] == kit and int(row["chem"]) == chem:
            return row

    raise ValueError(f"Kit '{kit}' chem {chem} not found in {kits_file}")


def _load_bc_rows(parse_info_dir: Path, bc_file: str) -> list[dict[str, str]]:
    csv_path = parse_info_dir / "barcodes" / f"bc_data_{bc_file}.csv"
    with open(csv_path) as f:
        return list(csv.DictReader(f))


def generate_parse_configs(
    kit_name: str,
    parse_info_dir: Path,
    output_dir: Path,
    wells: list[str] | None = None,
    logger: logging.Logger | None = None,
) -> str:
    """Generate Parse config files for a given kit and optional well filter.

    Returns the kb-python x_string for the kit.
    """
    kit, chem = _parse_kit_name(kit_name)
    kit_info = _load_kit_info(parse_info_dir, kit, chem)

    bc1_rows = _load_bc_rows(parse_info_dir, kit_info["bc1"])
    bc2_rows = _load_bc_rows(parse_info_dir, kit_info["bc2"])
    bc3_rows = _load_bc_rows(parse_info_dir, kit_info["bc3"])

    # bc1: T-type (polyT) and R-type (randO), in CSV order
    bc1_T = [(r["sequence"], r["well"]) for r in bc1_rows if r["stype"] == "T"]
    bc1_R = [(r["sequence"], r["well"]) for r in bc1_rows if r["stype"] == "R"]
    # bc2 / bc3: L-type library barcodes
    bc2_seqs = [r["sequence"] for r in bc2_rows if r["stype"] == "L"]
    bc3_seqs = [r["sequence"] for r in bc3_rows if r["stype"] == "L"]

    wells_set = set(wells) if wells else None

    output_dir.mkdir(parents=True, exist_ok=True)

    # r1_R.txt — R-type bc1 sequences for the selected wells, CSV order
    r1_R_seqs = [s for s, w in bc1_R if wells_set is None or w in wells_set]
    (output_dir / "r1_R.txt").write_text("\n".join(r1_R_seqs) + "\n")

    # r1_T.txt — T-type bc1 sequences for the selected wells, CSV order
    r1_T_seqs = [s for s, w in bc1_T if wells_set is None or w in wells_set]
    (output_dir / "r1_T.txt").write_text("\n".join(r1_T_seqs) + "\n")

    # bcs_to_wells.txt — all bc1 T-type then R-type, tab-separated seq\twell
    bcs_lines = [f"{s}\t{w}" for s, w in bc1_T] + [f"{s}\t{w}" for s, w in bc1_R]
    (output_dir / "bcs_to_wells.txt").write_text("\n".join(bcs_lines) + "\n")

    # onlist.txt — columns: bc3 | bc2 | bc1 (T+R combined, CSV order), padded with "-"
    bc1_all_seqs = [s for s, _ in bc1_T] + [s for s, _ in bc1_R]
    max_rows = max(len(bc3_seqs), len(bc2_seqs), len(bc1_all_seqs))
    onlist_lines = []
    for k in range(max_rows):
        col1 = bc3_seqs[k] if k < len(bc3_seqs) else "-"
        col2 = bc2_seqs[k] if k < len(bc2_seqs) else "-"
        col3 = bc1_all_seqs[k] if k < len(bc1_all_seqs) else "-"
        onlist_lines.append(f"{col1} {col2} {col3}")
    (output_dir / "onlist.txt").write_text("\n".join(onlist_lines) + "\n")

    # replace.txt — bc1 R-type → *T-type for the same well
    well_to_T = {w: s for s, w in bc1_T}
    replace_lines = [
        f"{r_seq}\t*{well_to_T[w]}"
        for r_seq, w in bc1_R
        if w in well_to_T
    ]
    (output_dir / "replace.txt").write_text("\n".join(replace_lines) + "\n")

    if logger:
        logger.info(
            "Generated Parse configs for %s (%d wells) in %s",
            kit_name,
            len(wells) if wells else len(bc1_T),
            output_dir,
        )

    return kit_info["x_string"]
