# Tests

End-to-end test of the Snakemake workflow on a tiny **synthetic** dataset. It runs
the whole pipeline — reference index → remultiplex → Parse barcode filter/split →
subsample → `kb count` → `adata.h5ad` — for both a 10x and a Parse assay, and checks
the count matrices come out non-empty. It needs no network or SRA access: the
reference and reads are generated from code.

## Run it

```bash
conda activate 10XvParse         # or any env with the workflow tools + anndata
bash Tests/run_e2e.sh
```

Knobs (environment variables):
- `CORES` — Snakemake cores (default 4).
- `KEEP_TEST_OUTPUT=1` — leave the generated `Data/Analysis_test/`, `Index/testsp/`,
  and `Generated/Analysis_test/` in place for inspection (otherwise removed on pass).
- `USE_CONDA=1` — run each rule in its own per-rule conda env (`workflow/envs/*.yaml`)
  via `--use-conda`, so those env specs get exercised too (needs conda ≥24.7.1 or mamba).

## What's here

| File | Role |
|------|------|
| `generate_test_data.py` | Synthesizes the mini reference (`Index/testsp/ref.*`) and the 10x + Parse reads (`Data/Analysis_test/<assay>/FASTA/Dumped/`), deterministically. |
| `e2e/analysis_test.yaml` | The test analysis config (a `local:` source, one 10x-vs-Parse comparison group). |
| `e2e/config.yaml` | Overlay merged onto `Config/config.yaml` via `--configfile`, adding the `Analysis_test` analysis and the `testsp` species. |
| `check_outputs.py` | Asserts the six target h5ads exist, are readable AnnData, and that the full-depth matrices are non-empty. |
| `run_e2e.sh` | Stages the data, runs the workflow for the test targets, checks the outputs. |
| `ci-env.yml` | Minimal conda env for CI (Snakemake, kb/kallisto/bustools, splitcode, seqtk, pigz, anndata). |

CI runs `run_e2e.sh` on every push/PR (`.github/workflows/e2e.yml`).

## How the synthetic data is built

- **Reference** — one ~3 kb chromosome with two 2-exon genes; committed as the
  `download_reference` outputs so Snakemake skips the download and `kb ref` builds a
  tiny index.
- **10x reads** — R1 = a real 10XV3 whitelist barcode (so kb's onlist correction
  keeps it) + a random UMI; R2 = a window of a transcript.
- **Parse reads** — two sublibraries of WT_mini_v3. The barcode read is laid out for
  the kit's x_string (UMI + bc3 + bc2 + bc1 with the kit's linker gaps), using
  bc1/bc2/bc3 drawn from `Resources/parse_info`. Both polyT (T-type bc1) and randO
  (R-type bc1) reads are emitted, so the barcode-split branch is exercised. The 4 bp
  sublibrary barcode is *not* in the raw reads — splitcode prepends it during
  remultiplexing, which is exactly what the test verifies.

Because the barcodes and positions are real, a passing run confirms not just that the
DAG is wired correctly but that the barcodes are recovered and reads pseudoaligned
(e.g. the current run yields 12 10x cells and 48 Parse cells across 2 genes).

## Scope

This covers the count-matrix pipeline (Tier 1). Gene-body coverage (STARsolo + RSeQC)
and the SRA/ERA download and hashtag/kite paths are not yet exercised here.
