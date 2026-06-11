# Import packages
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from upsetty import Upset
import gseapy as gp
import matplotlib
import matplotlib.patches as mpatches
from scipy.stats import gmean
import re
import edgepython as ep
from typing import Tuple, Any
from . import init_processing

def load_cross_comparison_data(samples: list[tuple[str, str, str, str, str]], comparison:Tuple[str,str], project_dir:Path)->pd.DataFrame:
    side1, side2 = comparison

    # Output dirs expected by the shared plotting helpers
    for d in [f"{side1}_outliers", f"{side2}_outliers"]:
        Path(d).mkdir(exist_ok=True)

    series1  = {}
    series2 = {}

    # Load h5ad matrices using init_processing from kb_python folder 
    for label, analysis, tenx_assay, parse_assay, tissue in samples:
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

def perform_edgepy(combined: pd.DataFrame, comparing: Tuple[str,str], samples: list[tuple[str, str, str, str, str]]) -> Tuple[Any, Any, Any, Any, Any]:
    side1, side2 = comparing
    
    labels      = [s[0] for s in samples]
    tissue_map  = {s[0]: s[4] for s in samples}
    tissues     = [tissue_map[l] for l in labels]
    tissue_categories=list(set(tissue_map.values()))  

    counts_matrix = combined.values         
    gene_names    = combined.index.tolist()
    genes_df      = pd.DataFrame({"gene_name": gene_names})

    # DGEList and TMM normalization
    dge = ep.make_dgelist(counts=counts_matrix, genes=genes_df)
    dge = ep.calc_norm_factors(dge, method="TMM")
    print("TMM norm factors:", np.round(dge["samples"]["norm.factors"].values, 4))

    # Design matrix: ~tissue + tech ()
    design_data = pd.DataFrame({
        "tissue": pd.Categorical(tissues * len(tissue_categories), categories=tissue_categories),
        "tech":   pd.Categorical([side1] * len(samples) + [side2] * len(samples), categories=[side1, side2]),
    })
    design = ep.model_matrix("~tissue+tech", data=design_data)
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
    
    return dge_f, fit, qlf, design, results

def volcano_plot(results: pd.DataFrame, comparing: Tuple[str, str], fdr_thresh: float, label_num: int = 8,
                 marker_genes_path: str | Path | None = None,
                 gene_sets: list[str] | None = None,
                 terms: list[str] | None = None):
    side1, side2 = comparing

    sig2 = (results["FDR"] < fdr_thresh) & (results["logFC"] > 0)
    sig1 = (results["FDR"] < fdr_thresh) & (results["logFC"] < 0)
    ns   = ~(sig2 | sig1)

    nlp = -np.log10(results["FDR"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    ax.scatter(results.loc[ns,   "logFC"], nlp[ns],    color="grey",    s=3, alpha=0.3)
    ax.scatter(results.loc[sig1, "logFC"], nlp[sig1],  color="#1f77b4", s=5, alpha=0.6)
    ax.scatter(results.loc[sig2, "logFC"], nlp[sig2],  color="#d62728", s=5, alpha=0.6)
    ax.axhline(-np.log10(fdr_thresh), color="black", linewidth=0.8, linestyle="--")
    ax.axvline(0, color="black", linewidth=0.5)

    for sub in [results[sig1].nsmallest(label_num, "FDR"),
                results[sig2].nsmallest(label_num, "FDR")]:
        for _, row in sub.iterrows():
            ax.annotate(row["gene_name"],
                        (row["logFC"], -np.log10(row["FDR"])),
                        fontsize=6, ha="center",
                        xytext=(0, 4), textcoords="offset points")

    legend_handles = [
        mpatches.Patch(color="#d62728", label=f"Higher in {side2} (n={sig1.sum()})"),
        mpatches.Patch(color="#1f77b4", label=f"Higher in {side1} (n={sig2.sum()})"),
        mpatches.Patch(color="grey",    label="n.s."),
    ]

    if marker_genes_path is not None:
        def _norm(name: str) -> str:
            for prefix in ("HUMAN_", "MOUSE_"):
                if name.startswith(prefix):
                    return name[len(prefix):].upper()
            return name.upper()

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

        for i, (term, gene_set) in enumerate(per_term_genes.items()):
            mask = norm_names.isin(gene_set)
            sub = results[mask]
            ax.scatter(sub["logFC"], nlp[mask], color="black", s=18, alpha=0.9, zorder=3)
            for _, row in sub.iterrows():
                ax.annotate(row["gene_name"],
                            (row["logFC"], -np.log10(row["FDR"])),
                            fontsize=6, ha="left", va="bottom",
                            xytext=(2, 2), textcoords="offset points",
                            clip_on=True)
            legend_handles.append(mpatches.Patch(color="black", label=term))

    ax.set_xlabel(f"log\u2082 fold change ({side2} / {side1})")
    ax.set_ylabel("\u2212log\u2081\u2080(FDR)")
    ax.set_title(f"Volcano: {side1} vs {side2} (FDR < {fdr_thresh})")
    ax.legend(handles=legend_handles, fontsize=8)
    plt.tight_layout()
    plt.show()

def plot_differential_metrics(results: pd.DataFrame, comparing: Tuple[str, str], fdr_thresh: float):
    side1, side2 = comparing

    sig = results["FDR"] < fdr_thresh
    pos = sig & (results["logFC"] > 0)   # higher in side 2
    neg = sig & (results["logFC"] < 0)   # higher in side 2

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    specs = [
        ("gene_length", "Gene length (log10 bp)",  True),
        ("gc_content",  "GC content (%)",           False),
    ]
    for ax, (col, xlabel, log_x) in zip(axes, specs):
        x = np.log10(results[col]) if log_x else results[col]
        ax.scatter(x[~sig], results.loc[~sig,  "logFC"], c="grey",    s=3, alpha=0.2, label="n.s.")
        ax.scatter(x[neg],  results.loc[neg,   "logFC"], c="#1f77b4", s=5, alpha=0.5, label=f"Higher in {side1}")
        ax.scatter(x[pos],  results.loc[pos,   "logFC"], c="#d62728", s=5, alpha=0.5, label=f"Higher in {side2}")
        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f"log\u2082 FC ({side2} / {side1})")
        ax.legend(fontsize=8, markerscale=2)

    plt.suptitle("logFC vs gene properties", fontsize=11)
    plt.tight_layout()
    plt.show()

def plot_differential_biotype(results: pd.DataFrame, comparing: Tuple[str, str]):
    side1, side2 = comparing

    def assign_biotype(row):
        if row.get("is_mito", False): return "mito"
        if row.get("is_ribo", False): return "ribo"
        if row.get("is_pc",   False): return "pc"
        if row.get("is_lnc",  False): return "lnc"
        return "other"

    results["biotype"] = results.apply(assign_biotype, axis=1)
    order = ["mito", "ribo", "pc", "lnc", "other"]

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

def _plot_term_enrichment_bar(ern_results: Tuple[pd.DataFrame, pd.DataFrame], comparing: Tuple[str, str]):
    side1, side2 = comparing

    DB_COLORS = {
        'GO_Biological_Process_2026': '#4e9af1',
        'GO_Molecular_Function_2026': '#7bc67e',
        'KEGG_2021_Human':            '#f4a261',
        'MSigDB_Hallmark_2020':       '#e76f51',
        'GO_Cellular_Component_2026': '#2a9d8f'
    }

    def shorten(term, maxlen=48):
        term = term.split(' (GO:')[0].strip()
        return term[:maxlen] + '…' if len(term) > maxlen else term

    def top_terms(df, n_per_db=4):
        sig = df[df['Adjusted P-value'] < 0.05].copy()
        sig['nlp'] = -np.log10(sig['Adjusted P-value'])
        sig['label'] = sig['Term'].apply(shorten)
        return (sig.groupby('Gene_set', group_keys=False)
                .apply(lambda g: g.nlargest(n_per_db, 'nlp'))
                .reset_index(drop=True))

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
                'GO_Cellular_Component_2026':  'GO Cellular Component',
                'KEGG_2021_Human':             'KEGG',
                'MSigDB_Hallmark_2020':        'MSigDB Hallmark'}
    ax.legend(handles=[mpatches.Patch(facecolor=c, label=db_labels[d]) for d, c in DB_COLORS.items()],
            fontsize=8, title='Database', loc='lower right')

    plt.tight_layout()
    plt.show()
    plt.close(fig)

def enrichment_analysis(combined: pd.DataFrame, results: pd.DataFrame, comparing: Tuple[str,str], fdr_thresh: float, enr_adjP_thresh: float = 0.05):
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
                    'KEGG_2021_Human', 'MSigDB_Hallmark_2020', 'GO_Cellular_Component_2026'],
            background=background,   
            outdir=None,
            verbose=False,
        )

        enr_results = enr.results[enr.results['Adjusted P-value'] < enr_adjP_thresh].sort_values('Adjusted P-value')
        enr_results.to_csv(f"{side}_outliers/enrichment_results.csv")
        enrichment_results.append(enr_results)

    _plot_term_enrichment_bar(enrichment_results, comparing)

def find_sig_prefixes(results: pd.DataFrame, comparing: Tuple[str,str], fdr_thresh: float, prefix_threshold: int = 5):
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
        ax.set_xlabel(f"Number of genes in {side} consistently enriched set")
        ax.set_title(f"Gene name prefixes (n={len(gene_list)} total, prefixes with >{prefix_threshold} genes)")
        plt.tight_layout()
        plt.show()
        plt.close(fig)