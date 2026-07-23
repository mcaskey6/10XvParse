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
from itertools import combinations
from typing import Tuple, Any
from . import init_processing
from .processing import _warn_unmatched, _load_orthologs, _normalize_gene_name

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
    ax.legend(handles=legend_handles, fontsize=8)
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

def _plot_term_enrichment_bar(ern_results: Tuple[pd.DataFrame, pd.DataFrame], comparing: dict[str, str], fdr_thresh: float=0.05, max_terms: int=20):
    side1, side2 = comparing.values()

    DB_COLORS = {
        'GO_Biological_Process_2026': '#4e9af1',
        'GO_Molecular_Function_2026': '#7bc67e',
        'GO_Cellular_Component_2026': '#2a9d8f'
    }

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
            color=[DB_COLORS[d] for d in side1_top['Gene_set']],  alpha=0.85)
    ax.barh(y_side2, -side2_top['nlp'],
            color=[DB_COLORS[d] for d in side2_top['Gene_set']], alpha=0.85)

    # tick labels
    all_y      = np.concatenate([y_side2, y_side1])
    all_labels = np.concatenate([side2_top['label'].values, side1_top['label'].values])
    ax.set_yticks(all_y)
    ax.set_yticklabels(all_labels, fontsize=8.5)

    ax.axvline(0, color='black', linewidth=0.8)
    ax.axhline(n_p + gap/2 - 0.5, color='#888888', linewidth=0.6, linestyle='--')

    # section labels
    xlim = max(side1_top['nlp'].max(), side2_top['nlp'].max())
    ax.text( xlim * 0.02, y_side1.mean(),  f"← {side1} enriched",  va='center', fontweight='bold', fontsize=9)
    ax.text(-xlim * 0.02, y_side2.mean(), f"{side2} enriched →", va='center', ha='right', fontweight='bold', fontsize=9)

    ax.set_xlabel('−log₁₀(adjusted p-value)')
    ax.set_title(f"Functional enrichment: genes consistently enriched in {side1} and {side2}", fontsize=10)

    db_labels = {'GO_Biological_Process_2026': 'GO Biological Process',
                'GO_Molecular_Function_2026':  'GO Molecular Function',
                'GO_Cellular_Component_2026':  'GO Cellular Component'}
    ax.legend(handles=[mpatches.Patch(facecolor=c, label=db_labels[d]) for d, c in DB_COLORS.items()],
            fontsize=8, title='Database', loc='lower right')

    plt.tight_layout()
    plt.show()
    plt.close(fig)

def enrichment_analysis(combined: pd.DataFrame, results: pd.DataFrame, comparing: dict[str, str], fdr_thresh: float, enr_adjP_thresh: float = 0.05, max_terms: int = 20):
    sig_side1   = results.loc[(results["FDR"] < fdr_thresh) & (results["logFC"] < 0), "gene_name"].tolist()
    sig_side2 = results.loc[(results["FDR"] < fdr_thresh) & (results["logFC"] > 0), "gene_name"].tolist()
    background = combined.index
    sig_genes = (sig_side1, sig_side2)

    enrichment_results = []
    for gene_list, side in zip(sig_genes, comparing):
        # Get genes consistently enriched
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=['GO_Biological_Process_2026', 'GO_Molecular_Function_2026',
                    'GO_Cellular_Component_2026'],
            background=background,
            outdir=None,
            verbose=False,
        )

        enr_results = enr.results[enr.results['Adjusted P-value'] < enr_adjP_thresh].sort_values('Adjusted P-value')
        enr_results.to_csv(f"{side}_outliers/enrichment_results.csv")
        enrichment_results.append(enr_results)

    _plot_term_enrichment_bar(enrichment_results, comparing, fdr_thresh=fdr_thresh, max_terms=max_terms)


def plot_goseq_enrichment_bar(comparison_dir: str | Path, comparing: dict[str, str],
                              pcol: str = "padj_wal", fdr_thresh: float = 0.05, max_terms: int = 20,
                              robust_only: bool = True):
    """Render goseq GC-corrected GO enrichment through the enrichr bar-plot style.

    Reads ``{side}_outliers/goseq_gc_corrected.csv`` (written by goseq.ipynb) for each side
    of the comparison, reshapes goseq's ``category`` (``db::term``) into the Gene_set / Term /
    Adjusted P-value columns ``_plot_term_enrichment_bar`` expects, then draws the same plot.

    For the corrected view (``pcol="padj_wal"``), ``robust_only`` (default) keeps only
    **GC-robust** terms — those also significant *before* correction (padj_naive < fdr_thresh).
    Without it the plot would also show terms that only crossed the threshold *after* GC
    reweighting (typically huge generic categories like "Nucleus"), which are correction
    artifacts, not enrichment that survived correction.

    Args:
        comparison_dir: Folder holding the ``{side}_outliers`` subdirs (e.g. .../PolyTvRandO).
        comparing: dict mapping side key -> display label; keys name the outlier folders.
        pcol: goseq p-value column to plot — "padj_wal" (GC-corrected, default) or "padj_naive".
        fdr_thresh, max_terms: passed through to the bar plot.
        robust_only: for the corrected view, require naive significance too (GC-robust).
    """
    comparison_dir = Path(comparison_dir)
    side1_key, side2_key = list(comparing)

    def _load(side):
        df = pd.read_csv(comparison_dir / f"{side}_outliers" / "goseq_gc_corrected.csv")
        if robust_only and pcol == "padj_wal":
            df = df[df["padj_naive"] < fdr_thresh]   # keep only terms also enriched pre-correction
        parts = df["category"].str.split("::", n=1)
        df["Gene_set"] = parts.str[0]
        df["Term"] = parts.str[1]
        df["Adjusted P-value"] = df[pcol]
        return df.sort_values("Adjusted P-value")

    _plot_term_enrichment_bar((_load(side1_key), _load(side2_key)), comparing,
                              fdr_thresh=fdr_thresh, max_terms=max_terms)

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