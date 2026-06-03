# Import packages
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from upsetty import Upset
import gseapy as gp
import matplotlib.patches as mpatches
from scipy.stats import gmean
import re
from typing import Tuple

def merge_analyses_with_clr(project_dir: Path, analyses_dict: dict[str:str], comparing: Tuple[str, str], pct_threshold: float = 0.5, pct_cutoff: float = 0.001) -> pd.DataFrame:
    type1 = comparing[0]
    type2 = comparing[1]

    type1_outlier_out_dir = Path(type1 + "_outliers")
    type2_outlier_out_dir = Path(type2 + "_outliers")
    type1_outlier_out_dir.mkdir(exist_ok=True)
    type2_outlier_out_dir.mkdir(exist_ok=True)

    # Load data
    data = {}
    for name, path in analyses_dict.items():
        data[name] = pd.read_csv(project_dir / path)
        data[name].drop(columns=["Unnamed: 0"], inplace=True)

    # Remove genes with low percent counts in all datasets
    for name, df in data.items():
        drop_indices = df[(df[f"{type1}_percent_counts"] < pct_cutoff) & (df[f"{type2}_percent_counts"] < pct_cutoff)].index
        df.drop(index=drop_indices, inplace=True)

    # Find outlier genes
    for name, df in data.items():
        CLR = np.log2((df[f"{type1}_percent_counts"] + 1e-6) / (df[f"{type2}_percent_counts"] + 1e-6)) - np.log2(gmean(df[f"{type1}_percent_counts"] + 1e-6) / gmean(df[f"{type2}_percent_counts"] + 1e-6))
        df["CLR"] = CLR
        df[f"{type1}_outliers"] = CLR >= np.log2(1 + pct_threshold)  
        df[f"{type2}_outliers"] = CLR <= np.log2(1/(1+pct_threshold)) 

    # Suffix dataset specific columns
    for name, df in data.items():
        for col_name in df.columns:
            if col_name.endswith("_distance") or col_name.endswith("_percent_counts") or col_name.endswith("_n_cells") or col_name.endswith("_outliers") or col_name in ["CLR"]:
                new_col_name = f"{name}_{col_name}"
                df.rename(columns={col_name: new_col_name}, inplace=True)

    merged_df = pd.DataFrame()
    for name, df in data.items():
        
        if merged_df.empty:
            merged_df = df
        else:
            merged_df = merged_df.merge(df, on=["gene_name", "gene_id", "gene_length", "gc_content",
                                        "is_lnc", "is_pc", "is_mito", "is_ribo"], how="outer")

    for col_name in merged_df.columns:
        if col_name.endswith("_distance") or col_name.endswith("_percent_counts") or col_name.endswith("_n_cells"):
            merged_df[col_name] = merged_df[col_name].fillna(0)
        if col_name.endswith("_outliers"):
            merged_df[col_name] = merged_df[col_name].fillna(False)
    merged_df.set_index("gene_name", inplace=True)

    return merged_df

def plot_merged_outlier_upset(merged_df: pd.DataFrame, comparing: Tuple[str, str]):
    def plot_upset(merged_df: pd.DataFrame, side: str):
        outlier_data = merged_df[[col for col in merged_df.columns if col.endswith(f"_{side}_outliers")]]
        outlier_data.columns = [col.replace(f"_{side}_outliers", "") for col in outlier_data.columns]
        outlier_data = outlier_data.astype(bool)
        outlier_data = outlier_data[outlier_data.any(axis=1)]
        upset = Upset.generate_plot(outlier_data)
        upset.show()

    plot_upset(merged_df, side=comparing[0])
    plot_upset(merged_df, side=comparing[1])

def consistent_outlier_violin_plots(merged_df: pd.DataFrame, comparing: Tuple[str, str], figsize: Tuple[int, int] = (10, 5)):
    background_lengths = merged_df['gene_length'].tolist()
    background_gcs = merged_df['gc_content'].tolist()

    for side in comparing:
        all5 = merged_df[
        merged_df[[c for c in merged_df.columns if c.endswith(f"_{side}_outliers")]].all(axis=1)
        ].index.tolist()

        lengths = merged_df[merged_df.index.isin(all5)]['gene_length'].tolist()
        gcs = merged_df[merged_df.index.isin(all5)]['gc_content'].tolist()

        fig, ax = plt.subplots(1, 2, figsize=figsize)
        ax[0].violinplot([np.log10(lengths), np.log10(background_lengths)], showextrema=False, showmedians=True)
        ax[1].violinplot([np.log10(gcs), np.log10(background_gcs)], showextrema=False, showmedians=True)
        ax[0].set_ylabel("Gene Length (log10)", fontsize=10)
        ax[1].set_ylabel("Percent GC content", fontsize=10)
        ax[0].set_xticks([1,2],labels=["Outliers", "Background"])
        ax[1].set_xticks([1,2],labels=["Outliers", "Background"])
        plt.tight_layout()
        plt.show()
        plt.close(fig)

def consistent_outlier_bar_plots(merged_df: pd.DataFrame, comparing: Tuple[str, str], figsize: Tuple[int, int] = (10, 5)):
    cols = ['is_pc', 'is_lnc', 'is_mito', 'is_ribo']
    col_names = ['protein-coding count', 'lncRNA count', 'mtRNA count', 'rRNA count']
    color = ['yellow', 'blue', 'red', 'green']

    for side in comparing:
        all5 = merged_df[
        merged_df[[c for c in merged_df.columns if c.endswith(f"_{side}_outliers")]].all(axis=1)
        ].index.tolist()

        fig, ax = plt.subplots(1,2, figsize=figsize)
        sums = []
        pcts = []
        for col in cols:
            sum = merged_df[col][merged_df.index.isin(all5)].sum()
            sums.append(sum)
            pcts.append(sum/merged_df[col].sum())

        ax[0].bar(col_names, sums, color=color)
        ax[1].bar(col_names, pcts, color=color)

        ax[0].set_ylabel('Count of Genes in Type')
        ax[1].set_ylabel('Percent of Genes in Type')

        plt.tight_layout()
        plt.show()
        plt.close(fig)

def consistent_outlier_enrichment(merged_df: pd.DataFrame, comparing: Tuple[str, str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # Background = all genes in your merged table
    background = merged_df.index.tolist()

    enrichment_results = []
    for side in comparing:
        # Get genes consistently enriched in all 5 datasets for this side
        all5 = merged_df[
            merged_df[[c for c in merged_df.columns if c.endswith(f'_{side}_outliers')]].all(axis=1)
        ].index.tolist()

        enr = gp.enrichr(
            gene_list=all5,
            gene_sets=['GO_Biological_Process_2026', 'GO_Molecular_Function_2026', 'GO_Cellular_Component_2026',
                    'KEGG_2021_Human', 'MSigDB_Hallmark_2020'],
            background=background,   
            outdir=None,
            verbose=False,
        )

        enr_results = enr.results[enr.results['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value')
        enr_results.to_csv(f"{side}_outliers/enrichment_results.csv")
        enrichment_results.append(enr_results)
    return enrichment_results[0], enrichment_results[1]

def consistent_outlier_prefixes(merged_df: pd.DataFrame, comparing: Tuple[str, str], prefix_threshold: int = 5):
    def extract_prefix(gene_name):
            # Leading letters up to first digit or hyphen: RPS27→RPS, MT-CO1→MT, ATP5F1A→ATP
            m = re.match(r'^([A-Z]+)', gene_name)
            return m.group(1) if m else gene_name
    for side in comparing:
        all5 = merged_df[
            merged_df[[c for c in merged_df.columns if c.endswith(f"_{side}_outliers")]].all(axis=1)
        ].index.tolist()

        prefixes = pd.Series(all5).apply(extract_prefix)
        prefix_counts = prefixes.value_counts()

        with open(f"{side}_outliers/enriched_prefixes.txt", "w") as f:
            for prefix in sorted(prefix_counts.index.tolist()):
                f.write(f"{prefix}\n")

        prefix_counts_filtered = prefix_counts[prefix_counts > prefix_threshold].sort_values(ascending=False)

        fig, ax = plt.subplots(figsize=(8, max(4, len(prefix_counts_filtered) * 0.35)))
        ax.barh(prefix_counts_filtered.index[::-1], prefix_counts_filtered.values[::-1])
        ax.set_xlabel(f"Number of genes in {side} consistently enriched set")
        ax.set_title(f"Gene name prefixes (n={len(all5)} total, prefixes with >{prefix_threshold} genes)")
        plt.tight_layout()
        plt.show()
        plt.close(fig)

def plot_term_enrichment_bar(ern_results: Tuple[pd.DataFrame, pd.DataFrame], comparing: Tuple[str, str]):
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
