# Import packages
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from upsetty import Upset
import gseapy as gp
import matplotlib
import matplotlib.patches as mpatches
from scipy.stats import gmean, fisher_exact, pearsonr, spearmanr
import statsmodels.formula.api as smf
import re
import edgepython as ep
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from itertools import combinations
from statsmodels.stats.multitest import multipletests
from typing import Tuple, Any
from . import init_processing
from .processing import _warn_unmatched, _load_orthologs, _normalize_gene_name
import scclr
from scipy import sparse

def load_cross_comparison_data(samples: list[tuple[str, str, str, str, str]], comparison: dict[str, str], project_dir:Path)->pd.DataFrame:
    side1, side2 = comparison

    # Output dirs expected by the shared plotting helpers
    for d in [f"{side1}_outliers", f"{side2}_outliers"]:
        Path(d).mkdir(exist_ok=True)

    series1  = {}
    series2 = {}

    # Load h5ad matrices using init_processing from kb_python folder 
    for label, analysis, tenx_assay, parse_assay in samples:
        print(f"Loading {label}...")
        datasets = []
        for side in comparison:
            assay = tenx_assay if side == "10x" else parse_assay
            d = init_processing(
                analysis_name=analysis, project_dir=project_dir,
                data_name=side, assay=assay, data_title=f"{label} {side}",
            )
            datasets.append(d)
        d1, d2 = datasets
        series1[f"{label}_{side1}"]    = pd.Series(np.asarray(d1.X.sum(axis=0)).flatten(),  index=d1.var_names)
        series2[f"{label}_{side2}"] = pd.Series(np.asarray(d2.X.sum(axis=0)).flatten(), index=d2.var_names)

    # Inner join: genes present in all samples
    df1  = pd.DataFrame(series1)
    df2 = pd.DataFrame(series2)
    combined = pd.concat([df1, df2], axis=1).dropna()
    combined = combined.loc[combined.sum(axis=1) > 0]
    print(f"\n{len(combined)} genes retained across all 10 samples")

    return combined

def scclr_pca(combined: pd.DataFrame, labels: list[str], comparing: dict[str,str]):
    side1, side2 = comparing

    comb_T = combined.T
    combined_side1 = comb_T[comb_T.index.str.endswith(side1)]
    combined_side2 = comb_T[comb_T.index.str.endswith(side2)]
    sclr_side1 = scclr.normalize(combined_side1, target="auto").sparse
    sclr_side2 = scclr.normalize(combined_side2, target="auto").sparse
    sclr = sparse.vstack([sclr_side1, sclr_side2])
    res  = scclr.pca(sclr, n_components=2)
    coords = res.scores
    var_explained = res.explained_variance_ratio * 100

    plt.figure(figsize=(10,7))
    plt.scatter(coords[:,0], coords[:,1])

    for i, label in enumerate(labels):
        plt.annotate(
            label, 
            xy=(coords[i,0], coords[i,1]),            
            xytext=(5, 5),               
            textcoords='offset points',  
            ha='left',                   
            va='bottom',                 
        )
    plt.title("scclr PCA Plot")
    plt.xlabel(f"Dimension 1 ({var_explained[0]:.1f}%)")
    plt.ylabel(f"Dimension 2 ({var_explained[1]:.1f}%)")
    plt.ylabel
    plt.show()

def perform_edgepy(combined: pd.DataFrame, comparing: dict[str, str], samples: list[tuple[str, str, str, str, str]]) -> Tuple[Any, Any, Any, Any, Any]:
    side1, side2 = comparing
    
    labels      = [s[0] for s in samples]
    anal_map  = {s[0]: s[1] for s in samples}
    analyses     = [anal_map[l] for l in labels]
    anal_categories=list(set(anal_map.values()))  

    counts_matrix = combined.values         
    gene_names    = combined.index.tolist()
    genes_df      = pd.DataFrame({"gene_name": gene_names})

    # DGEList and TMM normalization
    dge = ep.make_dgelist(counts=counts_matrix, genes=genes_df)
    dge = ep.calc_norm_factors(dge, method="TMM")
    print("TMM norm factors:", np.round(dge["samples"]["norm.factors"].values, 4))

    # Design matrix: ~analysis + tech. The count columns are the side1 block (all samples)
    # followed by the side2 block, so each per-sample list is tiled once per side.
    design_data = pd.DataFrame({
        "analyses": pd.Categorical(analyses * len(comparing), categories=anal_categories),
        "tech":   pd.Categorical([side1] * len(samples) + [side2] * len(samples), categories=[side1, side2]),
    })
    design = ep.model_matrix("~analyses+tech", data=design_data)
    print(f"Design: {design.shape[0]} samples x {design.shape[1]} coefficients  "
        f"(df.residual = {design.shape[0] - design.shape[1]})")
    
    # Filter low-count genes
    keep  = ep.filter_by_expr(dge, design=design)
    dge_f = dge.copy()
    dge_f["counts"] = dge["counts"][keep]
    dge_f["genes"]  = dge["genes"].iloc[keep]
    print(f"Retained {keep.sum()} / {len(keep)} genes after filterByExpr")

    # Estimate Dispersion
    dge_f = ep.estimate_disp(dge_f, design=design)
    print(f"Common dispersion (BCV): {np.sqrt(dge_f['common.dispersion']):.4f}")

    # Fit data to glm
    fit = ep.glm_ql_fit(dge_f, design=design)

    # Test the last design coefficient 
    qlf = ep.glm_ql_ftest(fit, coef=design.shape[1] - 1)

    # Get significant genes
    keep  = ep.filter_by_expr(dge, design=design)
    tt = ep.top_tags(qlf, n=keep.sum(), adjust_method="BH", sort_by="PValue")
    results = tt["table"].copy()

    # Save Results
    results.to_csv("results.csv", index=False)
    
    return dge_f, fit, qlf, design, results

def volcano_plot(results: pd.DataFrame, comparing: dict[str, str], fdr_thresh: float, label_num: int = 8,
                 marker_genes_path: str | Path | None = None,
                 gene_sets: list[str] | None = None,
                 terms: list[str] | None = None,
                 highlight_genes: list[str] | None = None,
                 highlight_label: str = "outside DE",
                 annotate: bool = True):
    side1, side2 = comparing.values()

    sig2 = (results["FDR"] < fdr_thresh) & (results["logFC"] > 0)
    sig1 = (results["FDR"] < fdr_thresh) & (results["logFC"] < 0)
    ns   = ~(sig2 | sig1)

    n_sig2, n_sig1 = int(sig2.sum()), int(sig1.sum())

    nlp = -np.log10(results["FDR"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    ax.scatter(results.loc[ns,   "logFC"], nlp[ns],    color="grey",    s=3, alpha=0.3)
    ax.scatter(results.loc[sig1, "logFC"], nlp[sig1],  color="#1f77b4", s=5, alpha=0.6)
    ax.scatter(results.loc[sig2, "logFC"], nlp[sig2],  color="#d62728", s=5, alpha=0.6)
    ax.axhline(-np.log10(fdr_thresh), color="black", linewidth=0.8, linestyle="--")
    ax.axvline(0, color="black", linewidth=0.5)

    for sub in [results[sig1].nsmallest(label_num, "FDR"),
                results[sig2].nsmallest(label_num, "FDR")]:
        if annotate:
            for _, row in sub.iterrows():
                ax.annotate(row["gene_name"],
                            (row["logFC"], -np.log10(row["FDR"])),
                            fontsize=6, ha="center",
                            xytext=(0, 4), textcoords="offset points")

    legend_handles = [
        mpatches.Patch(color="#d62728", label=f"Higher in {side2} (n={n_sig2})"),
        mpatches.Patch(color="#1f77b4", label=f"Higher in {side1} (n={n_sig1})"),
        mpatches.Patch(color="grey",    label="n.s."),
    ]

    if marker_genes_path is not None:
        orthologs = _load_orthologs()
        _norm = lambda name: _normalize_gene_name(name, orthologs)

        enr_df = pd.read_csv(marker_genes_path, index_col=0)
        if gene_sets is not None:
            enr_df = enr_df[enr_df["Gene_set"].isin(gene_sets)]
        if terms is not None:
            enr_df = enr_df[enr_df["Term"].isin(terms)]

        per_term_genes: dict[str, set[str]] = {}
        for _, row in enr_df.iterrows():
            t = row["Term"]
            if t not in per_term_genes:
                per_term_genes[t] = set()
            per_term_genes[t].update(_norm(g) for g in row["Genes"].split(";"))

        norm_names = results["gene_name"].apply(_norm)
        available = set(norm_names)
        for term, gene_set in per_term_genes.items():
            _warn_unmatched(gene_set, available, "volcano_plot", f'"{term}"', "results")

        def _volcano_fisher(mask: pd.Series) -> float:
            sub  = results[mask]
            rest = results[~mask]
            sa, sb = (sub["logFC"]  > 0).sum(), (sub["logFC"]  < 0).sum()
            ra, rb = (rest["logFC"] > 0).sum(), (rest["logFC"] < 0).sum()
            if sa + sb == 0:
                return 1.0
            _, p = fisher_exact([[sa, sb], [ra, rb]])
            return p

        def _fmt_p(p: float) -> str:
            return f"p={p:.2e}"

        def _frac(n, denom): return f"{n}/{denom} ({100 * n / denom:.1f}%)" if denom else f"{n}/0"

        def _summarize(label, mask, p):
            m = mask.values if hasattr(mask, "values") else mask
            n = int(m.sum())
            print(f"  {label} — {n} annotated genes labeled on plot ({_fmt_p(p)}):")
            print(f"     logFC > 0, higher in {side2}: {_frac(int((m & (results['logFC'] > 0)).sum()), n)}")
            print(f"     logFC < 0, higher in {side1}: {_frac(int((m & (results['logFC'] < 0)).sum()), n)}")
            print(f"     significant (FDR < {fdr_thresh}) higher in {side2}: {_frac(int((m & sig2).sum()), n)}")
            print(f"     significant (FDR < {fdr_thresh}) higher in {side1}: {_frac(int((m & sig1).sum()), n)}")

        print(f"Annotated marker genes (logFC = {side2} / {side1}):")
        union_mask = pd.Series(False, index=results.index)
        for i, (term, gene_set) in enumerate(per_term_genes.items()):
            mask = norm_names.isin(gene_set)
            union_mask |= mask
            p = _volcano_fisher(mask)
            sub = results[mask]
            ax.scatter(sub["logFC"], nlp[mask], color="black", s=18, alpha=0.9, zorder=3)
            for _, row in sub.iterrows():
                if annotate:
                    ax.annotate(row["gene_name"],
                                (row["logFC"], -np.log10(row["FDR"])),
                                fontsize=6, ha="left", va="bottom",
                                xytext=(2, 2), textcoords="offset points",
                                clip_on=True)
            legend_handles.append(mpatches.Patch(color="black", label=f"{term} ({_fmt_p(p)})"))
            _summarize(term, mask, p)

        if len(per_term_genes) > 1:
            _summarize("ALL terms (union)", union_mask, _volcano_fisher(union_mask))

    if highlight_genes is not None:
        orthologs = _load_orthologs()
        _hn = lambda name: _normalize_gene_name(name, orthologs)
        highlight_set = {_hn(g) for g in highlight_genes}
        mask = results["gene_name"].apply(_hn).isin(highlight_set)
        sub = results[mask]
        ax.scatter(sub["logFC"], nlp[mask], facecolors="none", edgecolors="black",
                   s=22, linewidths=0.8, alpha=0.9, zorder=4)
        legend_handles.append(
            mpatches.Patch(facecolor="none", edgecolor="black",
                           label=f"{highlight_label} (n={int(mask.sum())})")
        )

    ax.set_xlabel(f"log\u2082 fold change ({side2} / {side1})")
    ax.set_ylabel("\u2212log\u2081\u2080(FDR)")
    ax.set_title(f"Volcano: {side1} vs {side2} (FDR < {fdr_thresh})")
    ax.legend(handles=legend_handles, loc='lower left', fontsize=8)
    plt.tight_layout()
    plt.show()

def _match_polyA_metric(gene_names: pd.Series, polyA_path: str | Path,
                        value: str = "snr_sum",
                        snr_min: int = 1, snr_max: int | None = None) -> pd.Series:
    """Map a gene-name Series to per-gene values from a CSV keyed by Ensembl ID and symbol.

    The CSV's first column holds unversioned Ensembl gene IDs and a "Name" column holds
    symbols. Current gene names are symbols for annotated genes and versioned Ensembl IDs
    (ENSG...N) for unannotated ones, so ENSG names are matched on the version-stripped ID
    and everything else on the uppercased symbol (duplicate symbols averaged).

    ``value`` selects the per-gene quantity over the SNR<n> columns whose index n falls in
    [snr_min, snr_max] (snr_max=None means no upper bound):
      "snr_sum"     — plain sum of counts, Σ SNRn.
      "snr_density" — run-length-weighted sum divided by gene length, Σ(n·SNRn) / Length,
                      i.e. the fraction of the transcript occupied by internal A-runs.
    Any other value is treated as a literal column name (e.g. "Proportion").

    Returns a float Series aligned to gene_names (NaN where unmatched).
    """
    tab = pd.read_csv(polyA_path)
    ens_col = tab.columns[0]
    if value in ("snr_sum", "snr_density"):
        snr_cols = [c for c in tab.columns
                    if (m := re.fullmatch(r"SNR(\d+)", str(c)))
                    and snr_min <= int(m.group(1)) <= (snr_max if snr_max is not None else np.inf)]
        counts = tab[snr_cols].to_numpy()
        if value == "snr_sum":
            tab = tab.assign(_v=counts.sum(axis=1))
        else:
            run_len = np.array([int(re.fullmatch(r"SNR(\d+)", str(c)).group(1)) for c in snr_cols])
            tab = tab.assign(_v=(counts * run_len).sum(axis=1) / tab["Length"].to_numpy())
    else:
        tab = tab.assign(_v=tab[value])

    ens = pd.Series(tab["_v"].values, index=tab[ens_col].astype(str))
    sym = tab.assign(_s=tab["Name"].astype(str).str.upper()).groupby("_s")["_v"].mean()

    def _lookup(name):
        s = str(name)
        if s.startswith("ENSG"):
            return ens.get(s.split(".")[0], np.nan)
        return sym.get(s.upper(), np.nan)

    return gene_names.map(_lookup).astype(float)


def plot_differential_metrics(results: pd.DataFrame, comparing: dict[str, str], fdr_thresh: float,
                              polyA_path: str | Path | None = None):
    side1, side2 = comparing.values()

    sig = results["FDR"] < fdr_thresh
    pos = sig & (results["logFC"] > 0)   # higher in side 2
    neg = sig & (results["logFC"] < 0)   # higher in side 2

    # Each spec carries the x-values directly (a Series aligned to results), so an
    # externally-matched metric can be added without stashing a column on results.
    specs = [
        (results["gene_length"], "Median Transcript Length (log10 bp)", True),
        (results["gc_content"],  "Median Transcript GC content (%)",     False),
    ]
    if polyA_path is not None:
        snr = _match_polyA_metric(results["gene_name"], polyA_path, value="snr_density", snr_min=8)
        n_matched = int(snr.notna().sum())
        snr = snr.where(snr != 0)   # drop genes with no internal A-runs (score 0)
        print(f"Internal polyA priming: {int(snr.notna().sum())} genes with nonzero density "
              f"(of {n_matched} matched)")
        specs.append((snr, "Internal A-run density (Σn·SNR≥8 / length)", False))

    fig, axes = plt.subplots(1, len(specs), figsize=(6.5 * len(specs), 5))
    axes = np.atleast_1d(axes)
    for ax, (xdata, xlabel, log_x) in zip(axes, specs):
        x = np.log10(xdata) if log_x else xdata
        ax.scatter(x[~sig], results.loc[~sig,  "logFC"], c="grey",    s=3, alpha=0.2, label="n.s.")
        ax.scatter(x[neg],  results.loc[neg,   "logFC"], c="#1f77b4", s=5, alpha=0.5, label=f"Higher in {side1}")
        ax.scatter(x[pos],  results.loc[pos,   "logFC"], c="#d62728", s=5, alpha=0.5, label=f"Higher in {side2}")

        # Univariate line of best fit (logFC ~ x) over all genes with a value
        fit_df = pd.DataFrame({"x": x, "y": results["logFC"]}).replace(
            [np.inf, -np.inf], np.nan).dropna()
        uni_fit = smf.ols("y ~ x", data=fit_df).fit()
        x_range = np.linspace(fit_df["x"].min(), fit_df["x"].max(), 200)
        ax.plot(x_range, uni_fit.predict(pd.DataFrame({"x": x_range})),
                color="black", linewidth=1.5,
                label=f"fit: slope={uni_fit.params['x']:+.3f}, R\u00b2={uni_fit.rsquared:.4f}")

        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f"log\u2082 FC ({side2} / {side1})")
        ax.legend(fontsize=8, markerscale=2)

    plt.suptitle("logFC vs gene properties", fontsize=11)
    plt.tight_layout()
    plt.show()

def plot_differential_biotype(results: pd.DataFrame, comparing: dict[str, str]):
    side1, side2 = comparing.values()

    def assign_biotype(row):
        if row.get("is_mito", False):   return "mito"
        if row.get("is_ribo", False):   return "ribo"
        if row.get("is_pseudo", False): return "pseudo"
        if row.get("is_tf", False):     return "tf"
        if row.get("is_pc", False):     return "pc"
        if row.get("is_lnc", False):    return "lnc"
        return "other"

    results["biotype"] = results.apply(assign_biotype, axis=1)
    order = ["mito", "ribo", "pseudo", "tf", "pc", "lnc", "other"]

    fig, ax = plt.subplots(figsize=(10, 5))
    vdata = [results.loc[results["biotype"] == bio, "logFC"].dropna().values for bio in order]
    ax.violinplot(vdata, positions=range(len(order)), showmedians=True, widths=0.6)

    y_top = max(results["logFC"].max(), 1) * 1.05
    for pos, data in enumerate(vdata):
        ax.text(pos, y_top, f"n={len(data)}", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_ylabel(f"log\u2082 FC ({side2} / {side1})")
    ax.set_title("logFC distribution by biotype")
    plt.tight_layout()
    plt.show()

def _plot_term_enrichment_bar(ern_results: Tuple[pd.DataFrame, pd.DataFrame], comparing: dict[str, str], fdr_thresh: float=0.05, max_terms: int=20,
                              title: str | None=None):
    side1, side2 = comparing.values()

    DB_COLORS = {
        'GO_Biological_Process_2026': '#4e9af1',
        'GO_Molecular_Function_2026': '#7bc67e',
        'GO_Cellular_Component_2026': '#2a9d8f'
    }
    DEFAULT_COLOR = '#999999'

    def shorten(term, maxlen=48):
        term = term.split(' (GO:')[0].strip()
        return term[:maxlen] + '…' if len(term) > maxlen else term

    def top_terms(df, n_top=max_terms):
        sig = df[df['Adjusted P-value'] < fdr_thresh].copy()
        sig['nlp'] = -np.log10(sig['Adjusted P-value'])
        sig['label'] = sig['Term'].apply(shorten)
        return sig.head(n_top)

    side1_top  = top_terms(ern_results[0]).sort_values('nlp')
    side2_top = top_terms(ern_results[1]).sort_values('nlp')

    n_t, n_p = len(side1_top), len(side2_top)
    gap = 1  # blank row between the two halves

    # y positions: second side at bottom, gap, first side at top
    y_side2 = np.arange(n_p)
    y_side1  = np.arange(n_p + gap, n_p + gap + n_t)

    fig, ax = plt.subplots(figsize=(13, (n_t + n_p) * 0.38 + 2))

    ax.barh(y_side1,  side1_top['nlp'],
            color=[DB_COLORS.get(d, DEFAULT_COLOR) for d in side1_top['Gene_set']],  alpha=0.85)
    ax.barh(y_side2, -side2_top['nlp'],
            color=[DB_COLORS.get(d, DEFAULT_COLOR) for d in side2_top['Gene_set']], alpha=0.85)

    # tick labels
    all_y      = np.concatenate([y_side2, y_side1])
    all_labels = np.concatenate([side2_top['label'].values, side1_top['label'].values])
    ax.set_yticks(all_y)
    ax.set_yticklabels(all_labels, fontsize=8.5)

    ax.axvline(0, color='black', linewidth=0.8)
    ax.axhline(n_p + gap/2 - 0.5, color='#888888', linewidth=0.6, linestyle='--')

    # Section labels. A side with zero significant terms leaves an empty nlp series, so
    # nanmax over the concatenation (with a floor) keeps xlim finite and the labels visible.
    all_nlp = np.concatenate([side1_top['nlp'].values, side2_top['nlp'].values])
    xlim = np.nanmax(all_nlp) if len(all_nlp) and not np.all(np.isnan(all_nlp)) else 1.0
    if n_t:
        ax.text( xlim * 0.02, y_side1.mean(),  f"← {side1} enriched",  va='center', fontweight='bold', fontsize=9)
    if n_p:
        ax.text(-xlim * 0.02, y_side2.mean(), f"{side2} enriched →", va='center', ha='right', fontweight='bold', fontsize=9)

    ax.set_xlabel('−log₁₀(adjusted p-value)')
    ax.set_title(title or f"Functional enrichment: genes consistently enriched in {side1} and {side2}",
                 fontsize=10)

    # Legend covers only the databases actually plotted, so a custom `libraries=` set
    # cannot raise a KeyError on an unknown Gene_set.
    db_labels = {'GO_Biological_Process_2026': 'GO Biological Process',
                'GO_Molecular_Function_2026':  'GO Molecular Function',
                'GO_Cellular_Component_2026':  'GO Cellular Component'}
    present = [d for d in dict.fromkeys(
        list(side1_top['Gene_set']) + list(side2_top['Gene_set']))]
    if present:
        ax.legend(handles=[mpatches.Patch(facecolor=DB_COLORS.get(d, DEFAULT_COLOR),
                                          label=db_labels.get(d, str(d)))
                           for d in present],
                  fontsize=8, title='Database', loc='lower right')

    plt.tight_layout()
    plt.show()
    plt.close(fig)

# ---------------------------------------------------------------------------
# GC-bias-corrected GO enrichment (goseq)
# ---------------------------------------------------------------------------

GO_LIBRARIES = ("GO_Biological_Process_2026",
                "GO_Molecular_Function_2026",
                "GO_Cellular_Component_2026")

_GOSEQ_R = Path(__file__).parent / "goseq_enrichment.R"
_RSCRIPT_CACHE: str | None = None
_GENE2CAT_CACHE: dict = {}


def _probe_rscript(candidate: str | Path) -> bool:
    """True if this Rscript exists and can load goseq."""
    try:
        p = subprocess.run(
            [str(candidate), "-e",
             'quit(status = as.integer(!requireNamespace("goseq", quietly=TRUE)))'],
            capture_output=True, text=True, timeout=120,
        )
        return p.returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def _resolve_rscript(rscript: str | Path | None = None) -> str:
    """Find an Rscript that has goseq installed.

    Order: explicit argument -> $XVP_RSCRIPT -> PATH -> sibling conda envs. The last
    step matters because the analysis env may ship its own Rscript *without* goseq, so
    a bare `which` can select the wrong interpreter; goseq typically lives in a separate
    env (e.g. `goseq`, from Envs/goseq.yaml). An explicit argument or $XVP_RSCRIPT that
    fails the probe raises rather than silently falling through. The successful result
    is memoised.
    """
    global _RSCRIPT_CACHE

    for label, cand in (("rscript argument", rscript),
                        ("$XVP_RSCRIPT", os.environ.get("XVP_RSCRIPT"))):
        if cand:
            if _probe_rscript(cand):
                return str(cand)
            raise RuntimeError(
                f"{label} points to {cand!r}, but that R interpreter cannot load goseq.\n"
                'Check with: Rscript -e \'library(goseq)\''
            )

    if _RSCRIPT_CACHE and _probe_rscript(_RSCRIPT_CACHE):
        return _RSCRIPT_CACHE

    tried = []
    candidates = []
    on_path = shutil.which("Rscript")
    if on_path:
        candidates.append(on_path)
    candidates += [str(p) for p in sorted(Path(sys.prefix).parent.glob("*/bin/Rscript"))]

    for cand in candidates:
        tried.append(cand)
        if _probe_rscript(cand):
            _RSCRIPT_CACHE = cand
            return cand

    raise RuntimeError(
        "No Rscript with goseq found. Tried: " + (", ".join(tried) or "(none)") + "\n"
        "Install it, e.g.:\n"
        "  conda create -n goseq -c conda-forge -c bioconda bioconductor-goseq r-biasedurn\n"
        "then point XvP at it:\n"
        "  export XVP_RSCRIPT=$CONDA_PREFIX/../goseq/bin/Rscript\n"
        'or pass rscript="/path/to/Rscript".'
    )


def _load_gene2cat(project_dir: Path, species: str = "human",
                   libraries: tuple[str, ...] = GO_LIBRARIES,
                   organism: str = "Human",
                   min_set_size: int = 0, max_set_size: int = 10**9,
                   overwrite: bool = False, verbose: bool = True) -> pd.DataFrame:
    """Build a gene -> category table from Enrichr gene-set libraries, cached as GMT.

    Caches each library to ``<project_dir>/Notebooks/gene_info/<species>/enrichr_gmt/<library>.gmt``
    (the same convention as ``query_ensembl``/``compute_gene_metrics``), so repeat runs need
    no network.

    NOTE ``max_set_size``: ``gp.get_library`` defaults to ``max_size=2000`` and silently
    DROPS larger sets. Several broad GO terms exceed that (e.g. "Nucleus" ~5300,
    "Regulation of Transcription by RNA Polymerase II" ~2300), so the default here is
    effectively unbounded. Lowering it will silently delete those terms from every result.

    Returns:
        DataFrame with columns gene (uppercase), category ("<library>::<term>"),
        Gene_set, Term.
    """
    key = (str(project_dir), species, libraries, organism, min_set_size, max_set_size)
    if not overwrite and key in _GENE2CAT_CACHE:
        return _GENE2CAT_CACHE[key]

    cache_dir = Path(project_dir) / "Notebooks" / "gene_info" / species / "enrichr_gmt"
    cache_dir.mkdir(parents=True, exist_ok=True)

    frames = []
    for lib in libraries:
        gmt = cache_dir / f"{lib}.gmt"
        if gmt.exists() and not overwrite:
            sets = gp.get_library(str(gmt), min_size=min_set_size, max_size=max_set_size)
        else:
            sets = gp.get_library(lib, organism=organism, min_size=min_set_size,
                                  max_size=max_set_size, save=str(gmt))
        for term, genes in sets.items():
            if not genes:
                continue
            frames.append(pd.DataFrame({
                "gene": pd.Series(genes, dtype=str).str.upper(),
                "Gene_set": lib,
                "Term": term,
            }))

    g2c = pd.concat(frames, ignore_index=True).drop_duplicates(["gene", "Gene_set", "Term"])
    g2c["category"] = g2c["Gene_set"] + "::" + g2c["Term"]
    g2c = g2c[["gene", "category", "Gene_set", "Term"]]

    if verbose:
        print(f"gene2cat: {g2c['category'].nunique()} categories, "
              f"{g2c['gene'].nunique()} genes, {len(g2c)} edges")

    _GENE2CAT_CACHE[key] = g2c
    return g2c


def _run_goseq(genes: pd.DataFrame, gene2cat: pd.DataFrame, rscript: str,
               verbose: bool = False) -> pd.DataFrame:
    """Run the goseq R script on a prepared universe. Returns its raw p-value table."""
    if genes["bias"].nunique() < 2:
        raise ValueError(
            "The bias column has no variation, so goseq's nullp() cannot fit a "
            "probability weighting function (mgcv's spline smoother requires spread). "
            "Check that the GC column was merged correctly."
        )

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        genes_csv, g2c_csv, out_csv = tmp / "genes.csv", tmp / "gene2cat.csv", tmp / "out.csv"
        genes[["gene", "de", "bias"]].to_csv(genes_csv, index=False)
        gene2cat[["gene", "category"]].to_csv(g2c_csv, index=False)

        # Strip R_* so an activated environment cannot drag a foreign library tree into
        # the goseq interpreter.
        env = {k: v for k, v in os.environ.items()
               if k not in ("R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE")}
        proc = subprocess.run(
            [rscript, str(_GOSEQ_R), str(genes_csv), str(g2c_csv), str(out_csv)],
            capture_output=True, text=True, env=env,
        )
        if proc.returncode != 0 or not out_csv.exists():
            tail = "\n".join((proc.stderr or "").strip().splitlines()[-40:])
            raise RuntimeError(
                f"goseq failed (exit {proc.returncode}) running {_GOSEQ_R}:\n{tail}"
            )
        if verbose and proc.stdout:
            print(proc.stdout.strip())
        return pd.read_csv(out_csv)


def goseq_enrichment_analysis(
    results: pd.DataFrame,
    comparing: dict[str, str],
    fdr_thresh: float,
    project_dir: Path,
    species: str = "human",
    bias_col: str = "gc_median_tx",
    exclude_cols: tuple[str, ...] = (),
    enr_adjP_thresh: float = 0.05,
    robust_only: bool = True,
    max_terms: int = 20,
    libraries: tuple[str, ...] = GO_LIBRARIES,
    min_set_size: int = 0,
    max_set_size: int = 10**9,
    out_dir: str | Path = ".",
    out_name: str | None = None,
    rscript: str | Path | None = None,
    plot: bool = True,
    verbose: bool = False,
) -> dict[str, pd.DataFrame]:
    """GC-bias-corrected GO enrichment via goseq (Young et al. 2010, Genome Biol 11:R14).

    Replaces plain Enrichr enrichment. goseq fits a probability weighting function
    P(gene is DE | GC), then tests each GO term with a Wallenius non-central
    hypergeometric test whose odds come from that curve — so a term full of GC-extreme
    genes must clear a correspondingly higher bar. The uncorrected hypergeometric test is
    reported alongside as the "naive" reference, which is what plain Enrichr computes.

    Terms are classified at ``enr_adjP_thresh``:
      - "GC-robust"  significant naively AND after correction -> orthogonal to GC bias
      - "GC-driven"  significant naively but NOT after correction -> a GC artifact
      - "n.s."       otherwise (includes terms that only pass *after* correction, which
                     are typically very broad categories and are correction artifacts)

    The universe is the genes that went through the differential test
    (``results["gene_name"]``), not an unfiltered background.

    ``exclude_cols`` names boolean columns on ``results`` (e.g. ``("is_ribo", "is_oxphos")``
    as defined for ``add_cell_metrics``); matching genes are dropped from the **entire
    universe**, not just the DE list. That is required for correctness: the PWF is fit over
    the universe, so retaining excluded genes with de=0 would misstate the bias curve and
    corrupt every p-value in the run. DE calls are NOT recomputed — they belong to the
    differential test, not the enrichment universe.

    Writes ``{out_dir}/{side}_outliers/{out_name}`` per side, using Enrichr's column names
    where they overlap so the output drops straight into ``_plot_term_enrichment_bar``.

    Args:
        results: from ``perform_edgepy``, with ``gene_name``, ``logFC``, ``FDR``, ``bias_col``
            and any ``exclude_cols`` merged in (logFC = side2 / side1).
        comparing: dict of side key -> display label; keys name the ``{side}_outliers`` dirs.
        fdr_thresh: DE significance threshold applied to ``results["FDR"]``.
        project_dir: project root (for the gene-set cache).
        bias_col: GC column used as the goseq bias variable.
        exclude_cols: boolean columns whose True genes are removed from the universe.
        enr_adjP_thresh: adjusted-p cutoff for term significance and classification.
        robust_only: plot only GC-robust terms.
        out_name: output filename; defaults to ``goseq_enrichment.csv``, or
            ``goseq_enrichment_excl_<...>.csv`` when ``exclude_cols`` is given.
        plot: draw the enrichment bar plot.

    Returns:
        dict of side key -> full per-term DataFrame (all tested terms, unfiltered).
    """
    side_keys = list(comparing)
    side1_key, side2_key = side_keys
    side1_lab, side2_lab = comparing[side1_key], comparing[side2_key]

    missing = [c for c in (bias_col,) + tuple(exclude_cols) if c not in results.columns]
    if missing:
        hint = (f"Run plotting.compute_gene_metrics(...) and merge {bias_col!r} into results"
                if bias_col in missing else
                f"Add {missing!r} to the usecols of the gene_attributes merge")
        raise KeyError(f"results is missing {missing}. {hint}.")

    if out_name is None:
        out_name = ("goseq_enrichment.csv" if not exclude_cols else
                    "goseq_enrichment_excl_"
                    + "_".join(c.removeprefix("is_") for c in exclude_cols) + ".csv")

    # --- Universe: genes that went through the differential test ---
    cols = ["gene_name", "logFC", "FDR", bias_col] + list(exclude_cols)
    df = results[cols].copy()
    df["gene"] = df["gene_name"].astype(str).str.upper()
    df = df.drop_duplicates("gene").dropna(subset=[bias_col])
    n_before = len(df)
    if exclude_cols:
        excl = df[list(exclude_cols)].astype("boolean").fillna(False).astype(bool)
        df = df[~excl.any(axis=1)]
        print(f"Excluded {n_before - len(df)} genes via {list(exclude_cols)}; "
              f"universe {n_before} -> {len(df)}")
    df = df.rename(columns={bias_col: "bias"})
    # Sort by gene for determinism: goseq's nullp() fits the PWF with mgcv's
    # monotonicity-constrained spline (pcls), which is mildly sensitive to row order when
    # the fit sits near its inequality constraints (~1e-5 relative on p-values). Sorting
    # makes results independent of however `results` happened to be ordered upstream.
    df = df.sort_values("gene").reset_index(drop=True)

    gene2cat = _load_gene2cat(project_dir, species=species, libraries=libraries,
                              min_set_size=min_set_size, max_set_size=max_set_size,
                              verbose=verbose)
    g2c = gene2cat[gene2cat["gene"].isin(set(df["gene"]))]

    rscript = _resolve_rscript(rscript)

    out: dict[str, pd.DataFrame] = {}
    tables: list[pd.DataFrame] = []
    for side_key, sign in ((side1_key, -1), (side2_key, +1)):
        de_mask = (df["FDR"] < fdr_thresh) & (
            df["logFC"] < 0 if sign < 0 else df["logFC"] > 0)
        genes = df.assign(de=de_mask.astype(int))
        de_genes = set(genes.loc[de_mask, "gene"])

        raw = _run_goseq(genes, g2c, rscript, verbose=verbose)

        tab = raw.rename(columns={"p_wallenius": "P-value", "p_hypergeom": "P-value (naive)"})
        parts = tab["category"].str.split("::", n=1)
        tab["Gene_set"], tab["Term"] = parts.str[0], parts.str[1]

        # BH within each ontology rather than pooled: the libraries differ enormously in
        # size and term-overlap structure, so a pooled correction is not a homogeneous family.
        for src, dst in (("P-value", "Adjusted P-value"),
                         ("P-value (naive)", "Adjusted P-value (naive)")):
            tab[dst] = np.nan
            for gs, idx in tab.groupby("Gene_set").groups.items():
                tab.loc[idx, dst] = multipletests(tab.loc[idx, src], method="fdr_bh")[1]

        sig_corr = tab["Adjusted P-value"] < enr_adjP_thresh
        sig_naive = tab["Adjusted P-value (naive)"] < enr_adjP_thresh
        tab["classification"] = np.where(sig_naive & sig_corr, "GC-robust",
                                  np.where(sig_naive & ~sig_corr, "GC-driven", "n.s."))

        # Odds ratio from the 2x2 (DE in cat / not) x (in cat / not), matching Enrichr's
        # effect-size column.
        n_de = int(de_mask.sum())
        n_univ = len(genes)
        a = tab["numDEInCat"].astype(float)
        b = n_de - a
        c = tab["numInCat"].astype(float) - a
        d = n_univ - n_de - c
        with np.errstate(divide="ignore", invalid="ignore"):
            tab["Odds Ratio"] = (a / b) / (c / d)

        # Genes = DE genes in the term, same semantics as Enrichr's Genes column.
        hits = (g2c[g2c["gene"].isin(de_genes)]
                .groupby("category")["gene"]
                .apply(lambda s: ";".join(sorted(set(s)))))
        tab["Genes"] = tab["category"].map(hits).fillna("")

        tab = tab[["Gene_set", "Term", "P-value", "Adjusted P-value",
                   "P-value (naive)", "Adjusted P-value (naive)", "Odds Ratio",
                   "numDEInCat", "numInCat", "classification", "Genes"]]
        # Deterministic tiebreak: BH ties are pervasive and the bar plot takes .head().
        tab = tab.sort_values(["Adjusted P-value", "P-value", "numDEInCat"],
                             ascending=[True, True, False]).reset_index(drop=True)

        odir = Path(out_dir) / f"{side_key}_outliers"
        odir.mkdir(parents=True, exist_ok=True)
        keep = ((tab["Adjusted P-value"] < enr_adjP_thresh) |
                (tab["Adjusted P-value (naive)"] < enr_adjP_thresh))
        tab[keep].to_csv(odir / out_name, index=False)

        n_rob = int((tab["classification"] == "GC-robust").sum())
        n_drv = int((tab["classification"] == "GC-driven").sum())
        print(f"  {side_key:<6}: {n_rob} GC-robust, {n_drv} GC-driven "
              f"(of {int(sig_naive.sum())} naively significant); "
              f"{n_de} DE genes, universe {n_univ} -> {odir / out_name}")

        out[side_key] = tab
        tables.append(tab[tab["classification"] == "GC-robust"] if robust_only else tab)

    if plot:
        excl = (" (excl. " + ", ".join(c.removeprefix("is_") for c in exclude_cols) + ")"
                if exclude_cols else "")
        _plot_term_enrichment_bar(
            (tables[0], tables[1]), comparing,
            fdr_thresh=enr_adjP_thresh, max_terms=max_terms,
            title=(f"GC-corrected GO enrichment (goseq): {side1_lab} vs {side2_lab}"
                   f"{excl} — {'GC-robust terms' if robust_only else 'all terms'}"),
        )

    return out


def plot_gc_correction_scatter(term_tables: dict[str, pd.DataFrame],
                               comparing: dict[str, str],
                               adjP_thresh: float = 0.05) -> None:
    """Naive vs GC-corrected term significance, one panel per side.

    Points below the diagonal were weakened by the GC correction. This is the plot that
    visually justifies the GC-robust / GC-driven split produced by
    ``goseq_enrichment_analysis`` (pass its return value straight in).
    """
    colors = {"GC-robust": "#2ca02c", "GC-driven": "#d62728", "n.s.": "#cccccc"}
    keys = list(term_tables)

    fig, axes = plt.subplots(1, len(keys), figsize=(6.0 * len(keys), 5.4), squeeze=False)
    for ax, key in zip(axes[0], keys):
        tab = term_tables[key]
        xn = -np.log10(tab["Adjusted P-value (naive)"].clip(lower=1e-300))
        yw = -np.log10(tab["Adjusted P-value"].clip(lower=1e-300))
        for cls in ("n.s.", "GC-driven", "GC-robust"):
            m = (tab["classification"] == cls).values
            ax.scatter(xn[m], yw[m], s=10, alpha=0.55, color=colors[cls],
                       label=f"{cls} (n={int(m.sum())})")
        lim = float(max(np.nanmax(xn.values, initial=1.0), np.nanmax(yw.values, initial=1.0))) * 1.05
        ax.plot([0, lim], [0, lim], color="black", linewidth=0.6, linestyle=":")
        thr = -np.log10(adjP_thresh)
        ax.axhline(thr, color="grey", linewidth=0.6, linestyle="--")
        ax.axvline(thr, color="grey", linewidth=0.6, linestyle="--")
        ax.set_xlabel("naive −log₁₀(BH p)")
        ax.set_ylabel("GC-corrected (Wallenius) −log₁₀(BH p)")
        ax.set_title(f"Higher in {comparing.get(key, key)}", fontsize=10)
        ax.legend(fontsize=7)

    plt.suptitle("GO enrichment before vs after GC correction "
                 "(below the diagonal = weakened by GC bias)", fontsize=10)
    plt.tight_layout()
    plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# Redundancy reduction: collapse GO terms that ride on the same DE genes
# ---------------------------------------------------------------------------

def _term_similarity(a: set, b: set, metric: str) -> float:
    """Overlap coefficient |A∩B|/min(|A|,|B|), or Jaccard |A∩B|/|A∪B|."""
    if not a or not b:
        return 0.0
    inter = len(a & b)
    if metric == "overlap":
        return inter / min(len(a), len(b))
    if metric == "jaccard":
        return inter / len(a | b)
    raise ValueError(f"metric must be 'overlap' or 'jaccard', got {metric!r}")


def _merge_groups(groups: list[list[int]], sets: list[set], merge_threshold: float,
                  merge_on: str) -> list[list[int]]:
    """Iteratively merge groups that share more than ``merge_threshold`` of their content.

    This is the step that makes DAVID's fuzzy multiple-linkage partitioning usable: without
    it, maximal-clique enumeration on a dense similarity graph returns many near-duplicate
    groups covering the same genes. Merging is agglomerative — the most similar pair above
    the threshold is merged, and similarities are recomputed — until nothing qualifies.

    ``merge_on="genes"`` compares the union of each group's gene sets (targets duplicated
    gene coverage directly); ``merge_on="members"`` compares the sets of member terms, which
    is closer to DAVID's literal definition.
    """
    groups = [list(g) for g in groups]

    def content(g):
        return set.union(*(sets[i] for i in g)) if merge_on == "genes" else set(g)

    while len(groups) > 1:
        conts = [content(g) for g in groups]
        best, best_pair = 0.0, None
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                a, b = conts[i], conts[j]
                s = len(a & b) / max(1, len(a | b))     # Jaccard on group content
                if s > best:
                    best, best_pair = s, (i, j)
        if best_pair is None or best < merge_threshold:
            break
        i, j = best_pair
        groups[i] = sorted(set(groups[i]) | set(groups[j]))
        groups.pop(j)
    return groups


def _collapse_one(term_table: pd.DataFrame, threshold: float, metric: str,
                  gene_col: str, min_size: int, max_nodes: int, max_cliques: int,
                  merge_threshold: float | None, merge_on: str, core_frac: float
                  ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Group one side's terms into overlapping (fuzzy) clusters of shared DE genes.

    Returns (representatives, memberships). See ``collapse_redundant_terms``.
    """
    import networkx as nx

    tab = term_table.copy()
    tab = tab[tab[gene_col].fillna("").astype(str) != ""].reset_index(drop=True)
    if tab.empty:
        empty_rep = term_table.iloc[0:0].copy()
        for c in ("cluster_id", "cluster_size", "cluster_members", "core_genes", "n_core_genes",
                  "common_genes", "n_common_genes", "union_genes", "n_union_genes"):
            empty_rep[c] = pd.Series(dtype=object)
        return empty_rep, pd.DataFrame(columns=["Term", "cluster_id", "n_clusters"])

    if len(tab) > max_nodes:
        raise ValueError(
            f"{len(tab)} terms exceeds max_nodes={max_nodes}. Maximal-clique enumeration is "
            "worst-case exponential; raise `threshold`, tighten `classification`, or raise "
            "`max_nodes` deliberately."
        )

    sets = [set(str(g).split(";")) for g in tab[gene_col]]
    n = len(sets)

    G = nx.Graph()
    G.add_nodes_from(range(n))          # isolated terms survive as singleton cliques
    for i in range(n):
        for j in range(i + 1, n):
            if _term_similarity(sets[i], sets[j], metric) >= threshold:
                G.add_edge(i, j)

    cliques = []
    for cl in nx.find_cliques(G):
        cliques.append(sorted(cl))
        if len(cliques) > max_cliques:
            raise ValueError(
                f"More than max_cliques={max_cliques} maximal cliques; the similarity graph is "
                "too dense. Raise `threshold` to sparsify it."
            )
    # Raw maximal cliques proliferate on a dense graph: many near-duplicate groups carve up
    # the same genes (e.g. four "Mitochondrial Inner Membrane" cliques whose gene unions are
    # 97-100% identical). DAVID avoids this by iteratively merging seed groups that share more
    # than a "multiple linkage threshold" of their content; replicate that here.
    adjp = tab["Adjusted P-value"].to_numpy()
    if merge_threshold is not None:
        cliques = _merge_groups(cliques, sets, merge_threshold, merge_on)

    # A term may still be the most significant member of several groups, which surfaces as the
    # same representative listed twice. Fold those together: the union keeps that member as its
    # minimum, so one pass leaves every group with a distinct representative.
    by_rep: dict[int, set] = {}
    for c in cliques:
        rep = min(c, key=lambda i: adjp[i])
        by_rep.setdefault(rep, set()).update(c)
    cliques = [sorted(v) for v in by_rep.values()]

    cliques = [c for c in cliques if len(c) >= min_size]

    # Most significant member represents each clique; order cliques by that member.
    cliques.sort(key=lambda c: (adjp[c].min(), -len(c)))

    rep_rows, mem_rows = [], []
    for cid, members in enumerate(cliques, start=1):
        rep_i = min(members, key=lambda i: adjp[i])
        core = set.intersection(*(sets[i] for i in members))
        union = set.union(*(sets[i] for i in members))
        # Strict intersection empties out once groups are merged (merged members need not all
        # pairwise overlap), so also report the genes carried by a majority of member terms --
        # this is what actually characterises a merged group.
        counts = Counter(g for i in members for g in sets[i])
        need = max(1, int(np.ceil(core_frac * len(members))))
        common = {g for g, c in counts.items() if c >= need}
        row = tab.iloc[rep_i].to_dict()
        row.update({
            "cluster_id": cid,
            "cluster_size": len(members),
            "cluster_members": ";".join(tab["Term"].iloc[i] for i in members),
            "core_genes": ";".join(sorted(core)),
            "n_core_genes": len(core),
            "common_genes": ";".join(sorted(common)),
            "n_common_genes": len(common),
            "union_genes": ";".join(sorted(union)),
            "n_union_genes": len(union),
        })
        rep_rows.append(row)
        for i in members:
            mem_rows.append({"Term": tab["Term"].iloc[i], "Gene_set": tab["Gene_set"].iloc[i],
                             "cluster_id": cid, "is_representative": i == rep_i,
                             "Adjusted P-value": adjp[i]})

    reps = pd.DataFrame(rep_rows)
    mem = pd.DataFrame(mem_rows)
    if not mem.empty:
        mem["n_clusters"] = mem.groupby("Term")["cluster_id"].transform("nunique")
    return reps, mem


def collapse_redundant_terms(
    term_tables: dict[str, pd.DataFrame],
    comparing: dict[str, str],
    threshold: float = 0.7,
    metric: str = "overlap",
    classification: str | None = "GC-robust",
    gene_col: str = "Genes",
    min_size: int = 1,
    out_dir: str | Path = ".",
    out_name: str | None = None,
    plot: bool = True,
    max_terms: int = 20,
    max_nodes: int = 200,
    max_cliques: int = 100_000,
    merge_threshold: float | None = 0.4,
    merge_on: str = "genes",
    core_frac: float = 0.5,
) -> dict[str, dict[str, pd.DataFrame]]:
    """Collapse GO terms that are significant because of the same differentially expressed genes.

    GO output is heavily redundant: nested parent/child terms, and terms that are semantically
    unrelated but annotated to the same genes, are reported as separate findings. For example
    "Negative Regulation of Myoblast Fusion" can come back highly significant in PBMC data while
    being 100% ribosomal proteins — the translation signal under a different label.

    Terms are grouped by the **overlap coefficient** of their DE gene sets (the ``Genes`` column,
    i.e. DE ∩ term). Overlap is used rather than Jaccard because GO redundancy is dominated by
    *nested* terms and Jaccard penalises the size difference: a 29-gene child fully contained in a
    135-gene parent scores Jaccard 0.21 but overlap 1.00.

    Groups are the **maximal cliques** of the graph whose edges are pairs at or above ``threshold``,
    so every pair within a group clears the threshold and a term may belong to more than one group
    (fuzzy membership, as in DAVID's fuzzy multiple-linkage partitioning). Isolated terms come back
    as singleton cliques, so every input term is covered.

    Note the overlap coefficient is not a metric (no triangle inequality). That does not affect the
    clique guarantee, which is defined directly on the pairwise values, but it is why this uses
    cliques rather than cutting a dendrogram.

    Args:
        term_tables: ``{side_key: term DataFrame}`` as returned by ``goseq_enrichment_analysis``.
        comparing: side key -> display label; keys name the ``{side}_outliers`` dirs.
        threshold: minimum pairwise similarity for two terms to share a group.
        metric: "overlap" (default) or "jaccard".
        classification: keep only rows with this ``classification`` (None keeps all).
        min_size: drop groups smaller than this (1 keeps singletons).
        merge_threshold: after clique enumeration, iteratively merge groups sharing more than
            this fraction (Jaccard) of their content. This is DAVID's multiple-linkage merging
            step; without it maximal cliques proliferate into near-duplicate groups covering the
            same genes. None disables merging (raw maximal cliques). Groups that end up sharing a
            representative term are always folded together, independently of this setting.
        merge_on: "genes" merges on the union of each group's gene sets (targets duplicated gene
            coverage directly); "members" merges on shared member terms, closer to DAVID's
            literal definition.
        core_frac: a gene is reported in ``common_genes`` when it appears in at least this
            fraction of a group's member terms. The strict intersection (``core_genes``) is
            usually empty for merged groups, so this is the useful summary of what a group is
            built on.
        out_name: defaults to ``goseq_enrichment_clusters.csv``.
        plot: draw the bar plot of group representatives.
        max_nodes, max_cliques: guards against pathological clique enumeration.

    Returns:
        ``{side_key: {"representatives": df, "memberships": df}}``. ``representatives`` keeps the
        input schema (so it feeds ``_plot_term_enrichment_bar`` unchanged) and adds ``cluster_id``,
        ``cluster_size``, ``cluster_members``, ``core_genes``/``n_core_genes`` (strict
        intersection), ``common_genes``/``n_common_genes`` (present in >= ``core_frac`` of members
        — the practical summary of what drives a merged group) and ``union_genes``/``n_union_genes``.
    """
    if out_name is None:
        out_name = "goseq_enrichment_clusters.csv"

    out: dict[str, dict[str, pd.DataFrame]] = {}
    reps_by_side = []
    for side_key in comparing:
        tab = term_tables[side_key]
        if classification is not None and "classification" in tab.columns:
            tab = tab[tab["classification"] == classification]
        reps, mem = _collapse_one(tab, threshold=threshold, metric=metric, gene_col=gene_col,
                                  min_size=min_size, max_nodes=max_nodes, max_cliques=max_cliques,
                                  merge_threshold=merge_threshold, merge_on=merge_on,
                                  core_frac=core_frac)

        odir = Path(out_dir) / f"{side_key}_outliers"
        odir.mkdir(parents=True, exist_ok=True)
        reps.to_csv(odir / out_name, index=False)

        n_multi = int((mem["n_clusters"] > 1).sum()) if not mem.empty else 0
        print(f"  {side_key:<6}: {len(tab)} terms -> {len(reps)} groups "
              f"({n_multi} term-memberships in >1 group) -> {odir / out_name}")
        out[side_key] = {"representatives": reps, "memberships": mem}
        reps_by_side.append(reps)

    if plot and len(reps_by_side) == 2:
        side1_lab, side2_lab = comparing.values()
        _plot_term_enrichment_bar(
            (reps_by_side[0], reps_by_side[1]), comparing,
            fdr_thresh=1.0, max_terms=max_terms,
            title=(f"GO enrichment, redundancy-collapsed ({metric} ≥ {threshold}): "
                   f"{side1_lab} vs {side2_lab}"),
        )

    return out


def find_sig_prefixes(results: pd.DataFrame, comparing: dict[str, str], fdr_thresh: float, prefix_threshold: int = 5):
    def extract_prefix(gene_name):
        # Leading letters up to first digit or hyphen: RPS27→RPS, MT-CO1→MT, ATP5F1A→ATP
        m = re.match(r'^([A-Z]+)', gene_name)
        return m.group(1) if m else gene_name

    sig_side1   = results.loc[(results["FDR"] < fdr_thresh) & (results["logFC"] < 0), "gene_name"].tolist()
    sig_side2 = results.loc[(results["FDR"] < fdr_thresh) & (results["logFC"] > 0), "gene_name"].tolist()
    sig_genes = (sig_side1, sig_side2)

    for gene_list, side in zip(sig_genes, comparing):
        gene_list = [g for g in gene_list if not g.startswith("ENSG")]
        prefixes = pd.Series(gene_list).apply(extract_prefix)
        prefix_counts = prefixes.value_counts()

        with open(f"{side}_outliers/enriched_prefixes.txt", "w") as f:
            for prefix in sorted(prefix_counts.index.tolist()):
                f.write(f"{prefix}\n")

        prefix_counts_filtered = prefix_counts[prefix_counts > prefix_threshold].sort_values(ascending=False)

        fig, ax = plt.subplots(figsize=(8, max(4, len(prefix_counts_filtered) * 0.35)))
        ax.barh(prefix_counts_filtered.index[::-1], prefix_counts_filtered.values[::-1])
        ax.set_xlabel(f"Number of genes in {comparing[side]} consistently enriched set")
        ax.set_title(f"Gene name prefixes (n={len(gene_list)} total, prefixes with >{prefix_threshold} genes)")
        plt.tight_layout()
        plt.show()
        plt.close(fig)


# Gene families that are depleted under fixation (cytoplasmic + mitochondrial ribosomal
# protein genes), as label -> regex on the normalized (uppercased) symbol.
RIBO_FAMILIES = {
    "RP (RPL/RPS)":    r"^RP[LS]",
    "MRP (MRPL/MRPS)": r"^MRP[LS]",
}
# OXPHOS subunit families by respiratory complex. ATP synthase (CV) is restricted to
# ATP5* — a bare ^ATP also matches unrelated ion-transport ATPases and washes out the signal.
# NB: the ATP5 symbols differ between HGNC vintages (e.g. ATP5B vs ATP5F1B), so ATP5 genes
# may not match by symbol across datasets even though they're caught within each one.
OXPHOS_FAMILIES = {
    "NDUF (CI)":   r"^NDUF",
    "SDH (CII)":   r"^SDH",
    "UQCR (CIII)": r"^UQCR",
    "COX (CIV)":   r"^COX",
    "ATP5 (CV)":   r"^ATP5",
}
_FAMILY_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd",
                  "#d62728", "#8c564b", "#e377c2", "#17becf"]


def _load_outside_consensus(
    outside_path: str | Path,
    outside_fdr_thresh: float,
    outside_comparisons: tuple[str, ...],
    gene_col: str,
) -> Tuple[pd.DataFrame, int]:
    """Load the outside DE sheet and return (consensus_df, n_total).

    The consensus set is genes significant (FDR < outside_fdr_thresh) in *every*
    comparison in outside_comparisons; ``outside_logFC`` is the mean of the per-comparison
    logFCs and ``_gene`` is the normalized symbol (duplicates collapsed, strongest |logFC|).
    """
    outside = pd.read_excel(outside_path)
    _ortho = _load_orthologs()
    outside["_gene"] = outside[gene_col].astype(str).apply(lambda g: _normalize_gene_name(g, _ortho))

    fdr_cols = [f"FDR.{c}" for c in outside_comparisons]
    lfc_cols = [f"logFC.{c}" for c in outside_comparisons]
    consensus = outside.loc[(outside[fdr_cols] < outside_fdr_thresh).all(axis=1)].copy()
    consensus["outside_logFC"] = consensus[lfc_cols].mean(axis=1)
    consensus = (consensus.reindex(consensus["outside_logFC"].abs()
                                   .sort_values(ascending=False).index)
                          .drop_duplicates("_gene"))
    return consensus, len(outside)


def compare_outside_de(
    results: pd.DataFrame,
    comparing: dict[str, str],
    outside_path: str | Path,
    fdr_thresh: float = 0.01,
    outside_fdr_thresh: float = 0.05,
    outside_comparisons: tuple[str, ...] = ("PRO", "BLA", "COL"),
    gene_col: str = "genes",
    families: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Test whether genes significant in an outside DE analysis correlate with the
    significant genes of this 10x-vs-Parse analysis.

    The outside file (e.g. an FF-vs-FFPE edgeR result) is expected to hold one sheet
    with a gene column and, for each comparison ``C`` in ``outside_comparisons``,
    columns ``logFC.C`` and ``FDR.C``. The "consensus" outside set is the genes
    significant (FDR < ``outside_fdr_thresh``) in *every* listed comparison; their
    consensus logFC is the mean of the per-comparison logFCs.

    Produces three outputs:
      1. A 2x2 overlap contingency (over genes testable in both analyses) + Fisher's
         exact test for over-representation, shown as an annotated heatmap.
      2. A logFC concordance scatter (this analysis vs outside) with Pearson/Spearman r.
      3. Returns the merged per-gene table for further inspection/export.

    Args:
        results: DataFrame from ``perform_edgepy`` with ``gene_name``, ``logFC``, ``FDR``.
        comparing: (side1, side2) naming this analysis's two groups.
        outside_path: Path to the outside DE .xlsx file.
        fdr_thresh: Significance threshold for this analysis.
        outside_fdr_thresh: Significance threshold applied to each outside comparison.
        outside_comparisons: Suffixes of the outside comparison columns to require.
        gene_col: Name of the gene-symbol column in the outside sheet.
        families: Optional label -> regex map; matching genes are recolored on the
            concordance scatter (e.g. ``RIBO_FAMILIES``). None => no family highlight.

    Returns:
        DataFrame of genes shared between the outside consensus set and this analysis's
        testable universe, with columns: gene_name, outside_logFC, logFC, FDR, sig_here.
    """
    side1, side2 = comparing.values()

    # --- Load outside DE results and build the consensus significant set ---
    consensus, n_outside = _load_outside_consensus(
        outside_path, outside_fdr_thresh, outside_comparisons, gene_col)
    consensus_genes = set(consensus["_gene"])

    print(f"Outside DE: {n_outside} genes loaded; "
          f"{len(consensus_genes)} significant in all of {list(outside_comparisons)} "
          f"(FDR < {outside_fdr_thresh}).")

    # --- Restrict to the universe of genes testable in this analysis ---
    res = results.copy()
    _ortho = _load_orthologs()
    res["_gene"] = res["gene_name"].astype(str).apply(lambda g: _normalize_gene_name(g, _ortho))
    res = res.drop_duplicates("_gene")
    universe = set(res["_gene"])
    consensus_in_universe = consensus_genes & universe
    print(f"This analysis tested {len(universe)} genes; "
          f"{len(consensus_in_universe)} of the outside consensus genes are testable here.")

    sig_here = set(res.loc[res["FDR"] < fdr_thresh, "_gene"])

    # --- (1) Overlap + Fisher's exact test ---
    a = len(sig_here & consensus_in_universe)            # sig here & outside
    b = len(sig_here - consensus_in_universe)            # sig here, not outside
    c = len(consensus_in_universe - sig_here)            # outside, not sig here
    d = len(universe) - a - b - c                        # neither
    table = np.array([[a, b], [c, d]])
    odds, pval = fisher_exact(table, alternative="greater")
    expected = len(sig_here) * len(consensus_in_universe) / max(len(universe), 1)
    print(f"\nOverlap: {a} genes are significant here AND in the outside consensus "
          f"(expected by chance ≈ {expected:.1f}).")
    print(f"Fisher's exact (one-sided, enrichment): odds ratio = {odds:.2f}, p = {pval:.3e}")

    fig, ax = plt.subplots(figsize=(4.6, 4))
    im = ax.imshow(table, cmap="Blues")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["outside\nconsensus", "not"])
    ax.set_yticks([0, 1]); ax.set_yticklabels([f"sig here\n(FDR<{fdr_thresh})", "not"])
    for (i, j), v in np.ndenumerate(table):
        ax.text(j, i, f"{v:,}", ha="center", va="center",
                color="white" if v > table.max() / 2 else "black", fontsize=11)
    ax.set_title(f"Gene overlap (universe n={len(universe):,})\n"
                 f"OR={odds:.2f}, Fisher p={pval:.2e}", fontsize=9)
    plt.tight_layout()
    plt.show()
    plt.close(fig)

    # --- (2) logFC concordance scatter ---
    merged = res.merge(consensus[["_gene", "outside_logFC"]], on="_gene", how="inner")
    merged["sig_here"] = merged["FDR"] < fdr_thresh

    x = merged["outside_logFC"].values
    y = merged["logFC"].values
    pear_r, pear_p = pearsonr(x, y)
    spear_r, spear_p = spearmanr(x, y)

    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
    ns_m = ~merged["sig_here"].values
    ax.scatter(x[ns_m], y[ns_m], s=10, alpha=0.4, color="grey", label="n.s. here")
    ax.scatter(x[~ns_m], y[~ns_m], s=14, alpha=0.7, color="#d62728",
               label=f"sig here (FDR<{fdr_thresh})")

    if families:
        for i, (label, pattern) in enumerate(families.items()):
            color = _FAMILY_COLORS[i % len(_FAMILY_COLORS)]
            fam_m = merged["_gene"].str.match(pattern).values
            ax.scatter(x[fam_m], y[fam_m], s=28, color=color, edgecolors="black",
                       linewidths=0.4, alpha=0.9, zorder=3,
                       label=f"{label} (n={int(fam_m.sum())})")

    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlabel("outside log₂ FC (FFPE / FF, mean over "
                  f"{','.join(outside_comparisons)})")
    ax.set_ylabel(f"this analysis log₂ FC ({side2} / {side1})")
    ax.set_title(f"logFC concordance on {len(merged)} shared genes\n"
                 f"Pearson r={pear_r:.2f} (p={pear_p:.1e}), "
                 f"Spearman ρ={spear_r:.2f} (p={spear_p:.1e})", fontsize=9)
    ax.legend(fontsize=8, markerscale=1.5)
    plt.tight_layout()
    plt.show()
    plt.close(fig)

    return merged[["gene_name", "outside_logFC", "logFC", "FDR", "sig_here"]]


def plot_family_logfc(
    results: pd.DataFrame,
    comparing: dict[str, str],
    outside_path: str | Path,
    families: dict[str, str] = RIBO_FAMILIES,
    outside_fdr_thresh: float = 0.05,
    outside_comparisons: tuple[str, ...] = ("PRO", "BLA", "COL"),
    gene_col: str = "genes",
) -> None:
    """Two-panel logFC distribution comparison for selected gene families.

    Left panel: this analysis (logFC = side2 / side1) over all tested genes.
    Right panel: the outside consensus set (logFC = FFPE / FF, mean over comparisons).
    Each panel shows a violin + median for every family plus an "all genes" reference,
    making the family-level shift (e.g. ribosomal depletion under fixation) directly visible.

    Args:
        results: DataFrame from ``perform_edgepy`` with ``gene_name`` and ``logFC``.
        comparing: (side1, side2) naming this analysis's two groups.
        outside_path: Path to the outside DE .xlsx file.
        families: label -> regex map on the normalized symbol (default ``RIBO_FAMILIES``).
        outside_fdr_thresh, outside_comparisons, gene_col: passed to the outside loader.
    """
    side1, side2 = comparing.values()

    res = results.copy()
    _ortho = _load_orthologs()
    res["_gene"] = res["gene_name"].astype(str).apply(lambda g: _normalize_gene_name(g, _ortho))
    consensus, _ = _load_outside_consensus(
        outside_path, outside_fdr_thresh, outside_comparisons, gene_col)

    group_labels = list(families) + ["all genes"]
    panels = [
        (res, "logFC", f"this analysis  (log₂ FC {side2} / {side1})"),
        (consensus, "outside_logFC", "outside consensus  (log₂ FC FFPE / FF)"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, (df, col, title) in zip(axes, panels):
        masks = [df["_gene"].str.match(p) for p in families.values()] + [pd.Series(True, index=df.index)]
        vdata = [df.loc[m, col].dropna().values for m in masks]
        parts = ax.violinplot(vdata, positions=range(len(group_labels)),
                              showmedians=True, widths=0.7)
        body_colors = [_FAMILY_COLORS[i % len(_FAMILY_COLORS)]
                       for i in range(len(families))] + ["#bbbbbb"]
        for pc, color in zip(parts["bodies"], body_colors):
            pc.set_facecolor(color); pc.set_alpha(0.6)

        y_top = max((np.max(v) for v in vdata if len(v)), default=1)
        for pos, (data, m) in enumerate(zip(vdata, masks)):
            med = np.median(data) if len(data) else float("nan")
            ax.text(pos, y_top, f"n={len(data)}\nmed={med:+.2f}",
                    ha="center", va="bottom", fontsize=8)

        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_xticks(range(len(group_labels)))
        ax.set_xticklabels(group_labels, rotation=15, ha="right", fontsize=8)
        ax.set_ylabel("log₂ fold change")
        ax.set_title(title, fontsize=9)

    plt.suptitle("Ribosomal gene families are depleted under fixation in both analyses",
                 fontsize=11)
    plt.tight_layout()
    plt.show()
    plt.close(fig)


def plot_logfc_lm(results: pd.DataFrame, comparing: dict[str, str]) -> None:
    """Fit linear models of logFC on gene length and GC content for three annotation methods.

    Fits four OLS models using statsmodels:
      A (genomic span):    logFC ~ log10(length_genomic)  + gc_genomic
      B (median tx):       logFC ~ log10(length_median_tx) + gc_median_tx
      C (exon union):      logFC ~ log10(length_exon_union) + gc_exon_union
      D (combined):        all six covariates simultaneously

    Requires results to contain columns produced by compute_gene_metrics():
    length_genomic, gc_genomic, length_median_tx, gc_median_tx,
    length_exon_union, gc_exon_union.

    Produces three figures:
      1. Forest plot — coefficient ± 95% CI for each model (one panel per model).
      2. Scatter grid — logFC vs each covariate with the univariate fit overlaid
         (annotated with slope, p-value, and R²).
      3. Covariate correlation grid — pairwise scatters among the length/GC metrics
         (length vs GC per method, all length-vs-length pairs, all GC-vs-GC pairs)
         annotated with Pearson r, to inspect the collinearity behind the VIFs.

    Also prints a coefficient summary table for all four models.

    Args:
        results: DataFrame with logFC and the six gene metric columns.
        comparing: Two-element tuple (side1, side2) naming the comparison groups.
    """
    side1, side2 = comparing.values()

    metric_pairs = [
        ("length_genomic",   "gc_genomic",    "Genomic span"),
        ("length_median_tx", "gc_median_tx",  "Median transcript"),
        ("length_exon_union","gc_exon_union",  "Exon union"),
    ]
    all_cols = ["logFC"] + [c for p in metric_pairs for c in p[:2]]
    df = results[all_cols].dropna().copy()

    for length_col, _, _ in metric_pairs:
        log_col = f"log10_{length_col}"
        df[log_col] = np.log10(df[length_col].clip(lower=1))

    models_spec = [
        ("Genomic span",      f"logFC ~ log10_length_genomic + gc_genomic"),
        ("Median transcript", f"logFC ~ log10_length_median_tx + gc_median_tx"),
        ("Exon union",        f"logFC ~ log10_length_exon_union + gc_exon_union"),
        ("Combined",          (
            "logFC ~ log10_length_genomic + gc_genomic"
            " + log10_length_median_tx + gc_median_tx"
            " + log10_length_exon_union + gc_exon_union"
        )),
    ]

    fits = {}
    for label, formula in models_spec:
        fits[label] = smf.ols(formula, data=df).fit()

    # --- Figure 1: forest plots ---
    fig, axes = plt.subplots(1, 4, figsize=(18, 4), sharey=False)
    for ax, (label, _) in zip(axes, models_spec):
        fit = fits[label]
        coef_df = fit.conf_int()
        coef_df.columns = ["lo", "hi"]
        coef_df["coef"] = fit.params
        coef_df["pval"] = fit.pvalues
        coef_df = coef_df.drop("Intercept")

        y_pos = range(len(coef_df))
        colors = ["#d62728" if p < 0.05 else "#aaaaaa" for p in coef_df["pval"]]
        ax.barh(list(y_pos), coef_df["coef"],
                xerr=[coef_df["coef"] - coef_df["lo"], coef_df["hi"] - coef_df["coef"]],
                color=colors, ecolor="black", capsize=3, height=0.5)
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
        ax.set_yticks(list(y_pos))
        ax.set_yticklabels(coef_df.index, fontsize=8)
        ax.set_title(f"{label}\nN={len(df):,}  R²={fit.rsquared:.4f}", fontsize=9)
        ax.set_xlabel("Coefficient (effect on logFC)", fontsize=8)

    red_patch = mpatches.Patch(color="#d62728", label="p < 0.05")
    grey_patch = mpatches.Patch(color="#aaaaaa", label="p ≥ 0.05")
    fig.legend(handles=[red_patch, grey_patch], loc="upper right", fontsize=8)
    plt.suptitle(f"logFC ({side2} / {side1}) — linear model coefficients", fontsize=11)
    plt.tight_layout()
    plt.show()

    # --- Figure 2: scatter grid with univariate fit lines ---
    # Univariate lines show the actual marginal trend visible in each scatter.
    # These may differ from the bivariate forest-plot coefficients when length
    # and GC are correlated — that discrepancy reveals the confounding structure.
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    scatter_specs = [
        ("log10_length_genomic",   "gc_genomic",   "Genomic span"),
        ("log10_length_median_tx", "gc_median_tx", "Median transcript"),
        ("log10_length_exon_union","gc_exon_union", "Exon union"),
    ]
    for col_idx, (log_len_col, gc_col, title) in enumerate(scatter_specs):
        for row_idx, (xcol, xlabel) in enumerate([
            (log_len_col, f"log₁₀ length ({title})"),
            (gc_col,      f"GC content (%) ({title})"),
        ]):
            ax = axes[row_idx][col_idx]
            ax.scatter(df[xcol], df["logFC"], s=2, alpha=0.15, color="grey")

            uni_fit = smf.ols(f"logFC ~ {xcol}", data=df).fit()
            x_range = np.linspace(df[xcol].min(), df[xcol].max(), 200)
            y_pred = uni_fit.predict(pd.DataFrame({xcol: x_range}))
            slope = uni_fit.params[xcol]
            pval = uni_fit.pvalues[xcol]
            r2 = uni_fit.rsquared
            line_color = "#d62728" if pval < 0.05 else "#1f77b4"
            ax.plot(x_range, y_pred, color=line_color, linewidth=1.5,
                    label=f"slope={slope:+.3f}\np={pval:.2e}  R²={r2:.4f}")
            ax.axhline(0, color="black", linewidth=0.5)
            ax.set_xlabel(xlabel, fontsize=8)
            ax.set_ylabel(f"log₂ FC ({side2} / {side1})", fontsize=8)
            ax.legend(fontsize=7)


    plt.tight_layout()
    plt.show()

    # --- Figure 3: covariate correlation grid ---
    # Pairwise relationships among the model covariates themselves, to inspect the
    # collinearity that drives the VIFs reported in the diagnostics table below.
    # Row 1: length vs GC within each annotation method.
    # Row 2: length vs length across methods (all pairs).
    # Row 3: GC vs GC across methods (all pairs).
    len_cols = ["log10_length_genomic", "log10_length_median_tx", "log10_length_exon_union"]
    gc_cols  = ["gc_genomic", "gc_median_tx", "gc_exon_union"]
    short = {
        "log10_length_genomic":    "log₁₀ len (genomic)",
        "log10_length_median_tx":  "log₁₀ len (median tx)",
        "log10_length_exon_union": "log₁₀ len (exon union)",
        "gc_genomic":              "GC % (genomic)",
        "gc_median_tx":            "GC % (median tx)",
        "gc_exon_union":           "GC % (exon union)",
    }
    row_specs = [
        ("length vs GC (per method)", list(zip(len_cols, gc_cols)),  False),
        ("length vs length",         list(combinations(len_cols, 2)), True),
        ("GC vs GC",                 list(combinations(gc_cols, 2)),  True),
    ]

    fig, axes = plt.subplots(3, 3, figsize=(15, 13))
    for row_idx, (row_label, pairs, like) in enumerate(row_specs):
        for col_idx, (xcol, ycol) in enumerate(pairs):
            ax = axes[row_idx][col_idx]
            ax.scatter(df[xcol], df[ycol], s=2, alpha=0.12, color="grey")

            r, p = pearsonr(df[xcol], df[ycol])
            cov_fit = smf.ols(f"{ycol} ~ {xcol}", data=df).fit()
            xr = np.linspace(df[xcol].min(), df[xcol].max(), 200)
            ax.plot(xr, cov_fit.predict(pd.DataFrame({xcol: xr})),
                    color="#d62728", linewidth=1.5)

            # For like-vs-like metrics (length/length, GC/GC) overlay y = x so the
            # deviation from perfect agreement between annotation methods is visible.
            if like:
                lo = min(df[xcol].min(), df[ycol].min())
                hi = max(df[xcol].max(), df[ycol].max())
                ax.plot([lo, hi], [lo, hi], color="black", linestyle=":",
                        linewidth=1, label="y = x")
                ax.legend(fontsize=7, loc="lower right")

            ax.set_xlabel(short[xcol], fontsize=8)
            ax.set_ylabel(short[ycol], fontsize=8)
            ax.set_title(f"Pearson r={r:.3f}  (p={p:.1e})", fontsize=8)

    plt.suptitle("Covariate correlations — collinearity among length / GC metrics",
                 fontsize=11)
    plt.tight_layout()
    plt.show()

    # --- Confounding diagnostics table ---
    # For each bivariate model (A/B/C), shows how each covariate's coefficient
    # and p-value change from univariate → bivariate, plus the VIF.
    # A coefficient that shrinks toward zero when the partner is added means
    # the univariate association was largely a proxy for the other variable.
    from statsmodels.stats.outliers_influence import variance_inflation_factor

    print(f"\n{'='*90}")
    print(f"  CONFOUNDING DIAGNOSTICS  (logFC ~ length + GC, three annotation methods)")
    print(f"{'='*90}")
    print(f"  N = {len(df):,}")
    print(
        f"\n  {'Method':<18} {'Covariate':<28}"
        f" {'Uni coef':>10} {'Uni p':>10} {'Uni R²':>8}"
        f" {'Biv coef':>10} {'Biv p':>10} {'VIF':>7}"
    )
    print(f"  {'─'*88}")

    for label, log_len_col, gc_col in [
        ("Genomic span",      "log10_length_genomic",   "gc_genomic"),
        ("Median transcript", "log10_length_median_tx", "gc_median_tx"),
        ("Exon union",        "log10_length_exon_union","gc_exon_union"),
    ]:
        biv_fit = fits[label]
        X_biv = df[[log_len_col, gc_col]].assign(const=1.0)
        vif_len = variance_inflation_factor(X_biv.values, 0)
        vif_gc  = variance_inflation_factor(X_biv.values, 1)

        for xcol, vif in [(log_len_col, vif_len), (gc_col, vif_gc)]:
            uni_fit = smf.ols(f"logFC ~ {xcol}", data=df).fit()
            uni_coef = uni_fit.params[xcol]
            uni_p    = uni_fit.pvalues[xcol]
            uni_r2   = uni_fit.rsquared
            biv_coef = biv_fit.params[xcol]
            biv_p    = biv_fit.pvalues[xcol]

            short_name = xcol.replace("log10_length_", "log10_len_")
            print(
                f"  {label:<18} {short_name:<28}"
                f" {uni_coef:>10.4f} {uni_p:>10.2e} {uni_r2:>8.4f}"
                f" {biv_coef:>10.4f} {biv_p:>10.2e} {vif:>7.2f}"
            )


def plot_intron_gc_decomposition(results: pd.DataFrame, gc_col: str = "gc_genomic") -> pd.DataFrame:
    """Show that the GC–length correlation lives in the intronic, not exonic, sequence.

    Genomic span ≈ exon-union length + intron content. This decomposes the negative
    GC vs length correlation by regressing GC against each length component, confirming
    the isochore mechanism: GC-poor (AT-rich) regions carry long introns, while exonic
    (mRNA) length is roughly orthogonal to regional base composition.

    Restricted to multi-exon genes (intron content > 0). Produces a 1x3 scatter row
    (GC vs genomic span / intron content / exon-union length) annotated with Pearson r,
    and prints the correlation of GC against all four length measures.

    Requires columns: length_genomic, length_exon_union, length_median_tx, and ``gc_col``.

    Returns:
        The per-gene frame used (with the derived ``intron_len`` column).
    """
    cols = ["length_genomic", "length_exon_union", "length_median_tx", gc_col]
    df = results[cols].dropna().copy()
    df["intron_len"] = df["length_genomic"] - df["length_exon_union"]
    df = df[df["intron_len"] > 0]   # multi-exon genes only

    specs = [
        ("length_genomic",    "log₁₀ genomic span"),
        ("intron_len",        "log₁₀ intron content (span − exon union)"),
        ("length_exon_union", "log₁₀ exon-union length"),
        ("length_median_tx",  "log₁₀ median-tx length"),
    ]

    print(f"  N multi-exon genes = {len(df):,}    (GC = {gc_col})")
    print(f"  {'GC vs ...':<42}{'Pearson r':>10}{'p':>12}")
    print(f"  {'─'*62}")
    rvals = {}
    for col, label in specs:
        r, p = pearsonr(np.log10(df[col].clip(lower=1)), df[gc_col])
        rvals[col] = r
        print(f"  {label:<42}{r:>10.3f}{p:>12.1e}")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, (col, label) in zip(axes, specs[:3]):
        logx = np.log10(df[col].clip(lower=1))
        ax.scatter(logx, df[gc_col], s=2, alpha=0.1, color="grey")
        fit = smf.ols("y ~ x", data=pd.DataFrame({"x": logx, "y": df[gc_col]})).fit()
        xr = np.linspace(logx.min(), logx.max(), 200)
        ax.plot(xr, fit.predict(pd.DataFrame({"x": xr})), color="#d62728", linewidth=1.5)
        ax.set_xlabel(label, fontsize=8)
        ax.set_ylabel(f"GC % ({gc_col})", fontsize=8)
        ax.set_title(f"Pearson r = {rvals[col]:.3f}", fontsize=9)

    plt.suptitle("GC–length correlation is carried by intron content, not exonic length",
                 fontsize=11)
    plt.tight_layout()
    plt.show()
    plt.close(fig)

    return df