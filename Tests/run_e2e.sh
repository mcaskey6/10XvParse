#!/usr/bin/env bash
# End-to-end test for the 10XvParse Snakemake workflow on a tiny synthetic dataset.
#
# Stages generated inputs, runs the full pipeline (reference index -> multiplex ->
# barcode filter/split -> subsample -> kb count) for a test analysis, and asserts the
# count matrices are produced. Uses only local inputs, so no network / SRA download.
#
# Requires the workflow tools on PATH (snakemake, kb/kallisto/bustools, splitcode,
# seqtk, pigz) plus python with anndata — see Tests/ci-env.yml. Run from anywhere.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CORES="${CORES:-4}"

# With USE_CONDA=1, run each rule in its own per-rule conda env (workflow/envs/*.yaml)
# so those env specs are exercised too; otherwise use the already-active environment.
CONDA_FLAGS=()
if [ -n "${USE_CONDA:-}" ]; then CONDA_FLAGS=(--use-conda); fi

# Everything the test creates, removed up front so each run is clean and repeatable.
TEST_DIRS=(Data/Analysis_test Index/testsp Generated/Analysis_test)

cleanup() { rm -rf "${TEST_DIRS[@]}"; }
trap 'echo "(leaving test artifacts in place for inspection)"' ERR

echo "== staging synthetic reference + reads =="
cleanup
python Tests/generate_test_data.py

echo "== running workflow =="
TARGETS=(
  Data/Analysis_test/10x/kb_python/10x_out/counts_unfiltered/adata.h5ad
  Data/Analysis_test/10x/kb_python/sampled_10x_out/counts_unfiltered/adata.h5ad
  Data/Analysis_test/parse/kb_python/parse_out/counts_unfiltered/adata.h5ad
  Data/Analysis_test/parse/kb_python/sampled_parse_out/counts_unfiltered/adata.h5ad
  Data/Analysis_test/parse/kb_python/sampled_polyT_out/counts_unfiltered/adata.h5ad
  Data/Analysis_test/parse/kb_python/sampled_randO_out/counts_unfiltered/adata.h5ad
)
# `--` terminates the variadic --rerun-triggers so it can't swallow the targets.
snakemake --configfile Tests/e2e/config.yaml --cores "$CORES" \
  "${CONDA_FLAGS[@]}" --rerun-triggers mtime -- "${TARGETS[@]}"

echo "== checking outputs =="
python Tests/check_outputs.py

echo "== PASS =="
# Keep artifacts unless KEEP_TEST_OUTPUT is set (useful for debugging in CI).
if [ -z "${KEEP_TEST_OUTPUT:-}" ]; then cleanup; fi
