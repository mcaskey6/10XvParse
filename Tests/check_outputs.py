"""Assert the end-to-end test produced the expected count matrices.

Checks every h5ad the test targets exists and is a readable AnnData, and that the
full-depth 10x and Parse matrices are non-empty (barcodes were recovered and reads
pseudoaligned). The sampled matrices are only checked for existence + readability,
since their cell/gene counts depend on the subsample depth.

Run: python Tests/check_outputs.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad

ROOT = Path(__file__).resolve().parent.parent
KB = ROOT / "Data" / "Analysis_test"

# (h5ad path, must_be_non_empty)
EXPECTED = [
    (KB / "10x/kb_python/10x_out/counts_unfiltered/adata.h5ad", True),
    (KB / "10x/kb_python/sampled_10x_out/counts_unfiltered/adata.h5ad", False),
    (KB / "parse/kb_python/parse_out/counts_unfiltered/adata.h5ad", True),
    (KB / "parse/kb_python/sampled_parse_out/counts_unfiltered/adata.h5ad", False),
    (KB / "parse/kb_python/sampled_polyT_out/counts_unfiltered/adata.h5ad", False),
    (KB / "parse/kb_python/sampled_randO_out/counts_unfiltered/adata.h5ad", False),
]


def main() -> int:
    failures = []
    for path, non_empty in EXPECTED:
        if not path.exists():
            failures.append(f"MISSING  {path.relative_to(ROOT)}")
            continue
        try:
            adata = ad.read_h5ad(path)
        except Exception as exc:  # noqa: BLE001
            failures.append(f"UNREADABLE {path.relative_to(ROOT)}: {exc}")
            continue
        n_obs, n_vars = adata.n_obs, adata.n_vars
        total = int(adata.X.sum()) if adata.X is not None else 0
        tag = "non-empty required" if non_empty else "existence only"
        print(f"OK  {path.relative_to(ROOT)}  ({n_obs} obs x {n_vars} vars, sum={total}) [{tag}]")
        if non_empty and (n_obs == 0 or total == 0):
            failures.append(f"EMPTY    {path.relative_to(ROOT)} (obs={n_obs}, sum={total})")

    if failures:
        print("\nFAILED:")
        for f in failures:
            print("  " + f)
        return 1
    print("\nAll expected h5ads present and valid.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
