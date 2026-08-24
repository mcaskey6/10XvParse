"""Parse Biosciences barcode logic for the 10XvParse workflow.

Turns a Parse kit name + the committed barcode tables (Configs/parse_info/) into
the config files the splitcode / kb count / STARsolo steps consume — the cell
onlist, the randO->polyT replace table, per-round STAR whitelists, and the
kb-python x_string — and registers the synthetic "4th-round" sublibrary barcode
that lets sublibraries share one FASTQ. Self-contained (standard library only);
imported by workflow/scripts/parse_configs.py and parse_star_params.py.
"""
from __future__ import annotations

import csv
import logging
import re
from pathlib import Path

_PARSE_PATTERN = re.compile(r"^.+_v\d+$")

# Length of the synthetic sublibrary barcode splitcode prepends when remultiplexing
# Parse libraries. Parse sublibraries cannot be concatenated without one: the sequencing
# index acts as a fourth cell barcode, so identical bc1/bc2/bc3 combinations from
# different sublibraries are different cells.
LIB_BC_LEN = 4


def is_parse_kit(tech_str: str) -> bool:
    return bool(_PARSE_PATTERN.match(tech_str))


def library_barcode(index: int, bclen: int = LIB_BC_LEN) -> str:
    """Reproduce the remultiplexing barcode splitcode assigns to the index-th sublibrary.

    splitcode encodes the sublibrary index (0-based, in order of first appearance in the
    batch file) as a 2-bit big-endian sequence with A=0, C=1, G=2, T=3, giving
    AAAA, AAAC, AAAG, AAAT, AACA, ... for --bclen=4.

    Computed here rather than read from splitcode's --mapping output, which reports the
    batch-file row index instead of the deduplicated sublibrary index and so is wrong
    whenever a sublibrary spans more than one row.
    """
    if index >= 4 ** bclen:
        raise ValueError(f"Sublibrary index {index} does not fit in a {bclen}-mer barcode")
    return "".join("ACGT"[(index >> (2 * shift)) & 0b11] for shift in reversed(range(bclen)))


def distinct_sublibraries(labels: list[str]) -> list[str]:
    """Sublibrary labels de-duplicated, in order of first appearance.

    Mirrors splitcode's batch file handling, where rows sharing a name are assigned
    the same barcode.
    """
    seen: list[str] = []
    for label in labels:
        if label not in seen:
            seen.append(label)
    return seen


def sublibrary_barcodes(labels: list[str], bclen: int = LIB_BC_LEN) -> list[str]:
    """Barcode for each distinct sublibrary label, in order of first appearance."""
    return [library_barcode(i, bclen) for i in range(len(distinct_sublibraries(labels)))]


def shift_x_string(x_string: str, bclen: int = LIB_BC_LEN, bc_file: int = 1) -> str:
    """Add the sublibrary barcode to a kb-python x_string.

    Prepends a barcode triplet covering the first `bclen` bases of the barcode read and
    shifts every other coordinate on that read right by `bclen`, to account for the
    sequence splitcode prepends when remultiplexing. Triplets on other files (the cDNA
    read) are left alone.
    """
    def shift_part(part: str) -> str:
        nums = [int(n) for n in part.split(",")]
        shifted = []
        for i in range(0, len(nums), 3):
            file_idx, start, end = nums[i:i + 3]
            if file_idx == bc_file:
                start, end = start + bclen, end + bclen
            shifted += [file_idx, start, end]
        return ",".join(str(n) for n in shifted)

    bc_part, umi_part, seq_part = x_string.split(":")
    bc_part = f"{bc_file},0,{bclen}," + shift_part(bc_part)
    return ":".join([bc_part, shift_part(umi_part), seq_part])


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
    lib_barcodes: list[str] | None = None,
    sublibraries: list[str] | None = None,
) -> str:
    """Generate Parse config files for a given kit and optional well filter.

    Pass `lib_barcodes` (from `sublibrary_barcodes`) to register the splitcode
    remultiplexing barcode as a fourth cell barcode: it becomes the leading column of
    onlist.txt, gets its own STARsolo whitelist, and the returned x_string is shifted
    accordingly. `sublibraries` is the matching list of labels, recorded for provenance.

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

    # onlist.txt — one column per barcode, in x_string order, padded with "-":
    # [sublibrary] | bc3 | bc2 | bc1 (T+R combined, CSV order)
    bc1_all_seqs = [s for s, _ in bc1_T] + [s for s, _ in bc1_R]
    columns = [bc3_seqs, bc2_seqs, bc1_all_seqs]
    if lib_barcodes:
        columns.insert(0, lib_barcodes)
    max_rows = max(len(col) for col in columns)
    onlist_lines = [
        " ".join(col[k] if k < len(col) else "-" for col in columns)
        for k in range(max_rows)
    ]
    (output_dir / "onlist.txt").write_text("\n".join(onlist_lines) + "\n")

    # replace.txt — bc1 R-type → *T-type for the same well
    well_to_T = {w: s for s, w in bc1_T}
    replace_lines = [
        f"{r_seq}\t*{well_to_T[w]}"
        for r_seq, w in bc1_R
        if w in well_to_T
    ]
    (output_dir / "replace.txt").write_text("\n".join(replace_lines) + "\n")

    # star_bc*.txt — per-round whitelists for STARsolo CB_UMI_Complex (one barcode per line)
    (output_dir / "star_bc3.txt").write_text("\n".join(bc3_seqs) + "\n")
    (output_dir / "star_bc2.txt").write_text("\n".join(bc2_seqs) + "\n")
    (output_dir / "star_bc1.txt").write_text("\n".join(bc1_all_seqs) + "\n")

    if not lib_barcodes:
        if logger:
            logger.info(
                "Generated Parse configs for %s (%d wells) in %s",
                kit_name,
                len(wells) if wells else len(bc1_T),
                output_dir,
            )
        return kit_info["x_string"]

    # lib_bc.txt — sublibrary barcode whitelist for STARsolo (first --soloCBwhitelist file)
    (output_dir / "lib_bc.txt").write_text("\n".join(lib_barcodes) + "\n")

    # sublibraries.txt — which barcode splitcode assigned to which sublibrary
    labels = distinct_sublibraries(sublibraries or [])
    labels += [f"sublibrary_{i}" for i in range(len(labels), len(lib_barcodes))]
    (output_dir / "sublibraries.txt").write_text(
        "\n".join(f"{bc}\t{label}" for bc, label in zip(lib_barcodes, labels)) + "\n"
    )

    if logger:
        logger.info(
            "Generated Parse configs for %s (%d wells, %d sublibraries) in %s",
            kit_name,
            len(wells) if wells else len(bc1_T),
            len(lib_barcodes),
            output_dir,
        )

    return shift_x_string(kit_info["x_string"], len(lib_barcodes[0]))


def x_string_to_bc1_location(x_string: str) -> str:
    """Return the splitcode location string for the bc1 (RT) barcode.

    Extracts the last barcode's (file_idx, start, end) from the x_string and
    formats it as '{file_idx},{start},{end}' for use in splitcode config files.
    File index is 0-based in both the x_string and splitcode location format.
    """
    bc_part = x_string.split(":")[0]
    nums = [int(n) for n in bc_part.split(",")]
    bc1_file, bc1_start, bc1_end = nums[-3], nums[-2], nums[-1]
    return f"{bc1_file},{bc1_start},{bc1_end}"


def x_string_to_star_params(x_string: str) -> tuple[list[str], str]:
    """Return (cb_positions, umi_position) for STARsolo CB_UMI_Complex.

    Converts kb-python half-open positions to STAR anchor-based positions.
    The kb file index is ignored — barcodes are always in the second file in
    --readFilesIn, and positions are measured from read start (anchor 0).
    Format: 0_start_0_(end-1)
    """
    bc_part, umi_part, *_ = x_string.split(":")
    nums = [int(n) for n in bc_part.split(",")]
    cb_positions = [
        f"0_{nums[i + 1]}_0_{nums[i + 2] - 1}"
        for i in range(0, len(nums), 3)
    ]
    u = [int(n) for n in umi_part.split(",")]
    return cb_positions, f"0_{u[1]}_0_{u[2] - 1}"
