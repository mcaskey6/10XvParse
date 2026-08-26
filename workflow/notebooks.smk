# Notebook execution: render each analysis notebook to an HTML report (convention C)
# while declaring the real data files it reads/writes (convention D), so the DAG is
# truthful and only stale notebooks rerun. Notebooks run in the ACTIVE conda env (they
# import XvP_utils + edgepython/scclr and shell out to the goseq R env), so these rules
# carry no `conda:` directive — run snakemake from the 10XvParse env.
#
# A notebook's h5ad inputs are inferred from the analysis it declares (analysis_h5ads,
# the full+sampled superset). Regular notebooks map to their folder's analysis; the
# cross-analysis Comparisons/Barnyard notebooks list their analyses explicitly below.

import glob as _glob
import os as _os

REPORT_DIR = "Reports"


def _rel(nb):
    return _os.path.relpath(nb, "Notebooks")


def _nb_html(nb):
    return f"{REPORT_DIR}/{_rel(nb)[:-len('.ipynb')]}.html"


def _slug(nb):
    return _rel(nb)[:-len(".ipynb")].replace("/", "_")


def _species_for(analyses):
    return sorted({ch.species(ANALYSES[a]) for a in analyses})


# The three edgepy comparison notebooks write results.csv (via perform_edgepy) into
# their own directory; the Analysis_2 combos read those three files back.
_EDGEPY = [
    "Notebooks/Comparisons/10XvParse/human_edgepy.ipynb",
    "Notebooks/Comparisons/10XvPolyT/human_edgepy.ipynb",
    "Notebooks/Comparisons/PolyTvRandO/human_edgepy.ipynb",
]

# Cross-analysis notebooks: the analyses whose matrices they consume (can't be inferred
# from a folder name) and whether they read the per-species gene_attributes cache. The
# edgepy SAMPLES pull from Analysis_3 + Analysis_5; barnyard is standalone (it downloads
# its own reference data, so it has no pipeline inputs).
_COMPARISONS = {
    "Notebooks/Comparisons/10XvParse/human_edgepy.ipynb":   dict(analyses=["Analysis_3", "Analysis_5"], gene_attr=True),
    "Notebooks/Comparisons/10XvPolyT/human_edgepy.ipynb":   dict(analyses=["Analysis_3", "Analysis_5"], gene_attr=True),
    "Notebooks/Comparisons/PolyTvRandO/human_edgepy.ipynb": dict(analyses=["Analysis_3", "Analysis_5"], gene_attr=True),
    "Notebooks/Comparisons/boxplots.ipynb":                 dict(analyses=["Analysis_2", "Analysis_3", "Analysis_4", "Analysis_5", "Analysis_6"], gene_attr=True),
    "Notebooks/Comparisons/human_pca_comparison.ipynb":     dict(analyses=["Analysis_3", "Analysis_5", "Analysis_6"], gene_attr=False),
    "Notebooks/Barnyard_Comparison/barnyard.ipynb":         dict(analyses=[], gene_attr=False),
}


# Registry: notebook path -> {analyses, gene_attr, data_out, data_in}.
NOTEBOOKS = {}

# Regular per-analysis notebooks, discovered by folder. Each combo writes a
# gene_comparisons CSV; when an analysis has several combos (A2 standard/mini, A3
# H1/H2) they must not share the path, so the file is suffixed by the notebook name's
# tail after "combo" (combo.ipynb -> gene_comparisons.csv, combo_H1.ipynb ->
# gene_comparisons_H1.csv). The notebook sets the same FULL_GENE_INFO_PATH.
#
# The parse tech notebooks also read the generated Parse barcode files
# (Generated/<analysis>/<assay>/{r1_T,r1_R,replace,bcs_to_wells,...}.txt) directly, so
# those are declared inputs — they are gitignored intermediates that would not otherwise
# be rebuilt when the h5ads already exist. parse_configs is a cheap script rule (no read
# inputs), so listing every parse assay's config for the analysis is harmless.
def _parse_configs_of(name):
    cfg = ANALYSES[name]
    return [f for assay in _assays(cfg) if ch.assay_type(cfg, assay) == "parse"
            for f in paths.parse_generated_files(name, assay)]


for _name in ANALYSES:
    _parse_cfgs = _parse_configs_of(_name)
    for _nb in sorted(_glob.glob(f"Notebooks/{_name}/*.ipynb")):
        _base = _os.path.basename(_nb)
        _is_combo = _base.startswith("combo")
        _is_parse = _base.startswith("parse")
        _suffix = _base[len("combo"):-len(".ipynb")] if _is_combo else ""
        NOTEBOOKS[_nb] = dict(
            analyses=[_name],
            gene_attr=_is_combo,   # combos read the per-species gene_attributes cache
            data_out=[f"Notebooks/{_name}/gene_data/gene_comparisons{_suffix}.csv"] if _is_combo else [],
            data_in=_parse_cfgs if _is_parse else [],
        )

# Cross-analysis notebooks.
for _nb, _spec in _COMPARISONS.items():
    NOTEBOOKS[_nb] = dict(
        analyses=_spec["analyses"],
        gene_attr=_spec["gene_attr"],
        data_out=[f"{_os.path.dirname(_nb)}/results.csv"] if _nb in _EDGEPY else [],
        data_in=[],
    )

# The Analysis_2 combos overlay the edgepy results, so they read those results.csv.
for _nb in ("Notebooks/Analysis_2/combo.ipynb", "Notebooks/Analysis_2/combo_mini.ipynb"):
    if _nb in NOTEBOOKS:
        NOTEBOOKS[_nb]["data_in"] += [f"{_os.path.dirname(e)}/results.csv" for e in _EDGEPY]


def _nb_inputs(spec):
    ins = []
    for a in spec["analyses"]:
        ins += analysis_h5ads(a)
    if spec["gene_attr"]:
        ins += [paths.gene_attributes(sp) for sp in _species_for(spec["analyses"])]
    ins += spec["data_in"]
    return sorted(set(ins))


# One render rule per notebook. Inputs/outputs are computed eagerly (plain lists), so
# each rule closes over its own notebook rather than the loop variable.
for _nb, _spec in NOTEBOOKS.items():
    _html = _nb_html(_nb)
    _outs = [_html] + _spec["data_out"]
    _outdirs = sorted({_os.path.dirname(o) for o in _outs})

    rule:
        name:
            f"notebook_{_slug(_nb)}"
        input:
            _nb_inputs(_spec)
        output:
            _outs
        params:
            nb=_nb,
            outdir=_os.path.dirname(_html),
            outname=_os.path.basename(_html),
            mkdirs=" ".join(_outdirs),
        # Notebooks run heavy, all-core steps internally, so each claims every core:
        # this makes Snakemake execute them ONE AT A TIME (never two in parallel),
        # regardless of --cores, without starving them of threads.
        threads:
            workflow.cores
        log:
            f"Logs/notebooks/{_slug(_nb)}.log"
        shell:
            "exec 2> {log}; mkdir -p {params.mkdirs}; "
            "jupyter nbconvert --to html --execute --ExecutePreprocessor.timeout=-1 "
            "--output-dir {params.outdir} --output {params.outname} {params.nb}"


rule notebooks:
    """Render every notebook to reports/<path>.html (and its declared data outputs)."""
    input:
        [_nb_html(nb) for nb in NOTEBOOKS],
