from pathlib import Path
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.axes
import matplotlib.collections
import numpy as np
import scanpy as sc
import pandas as pd
import warnings
import urllib.request
import urllib.parse
import gseapy as gp
from scipy.stats import spearmanr, pearsonr, fisher_exact
from matplotlib.colors import LogNorm, Normalize
from typing import Tuple
from . import processing
from .processing import _warn_unmatched, _load_orthologs, _normalize_gene_name

_GO_LIBRARY_CACHE: dict[str, dict[str, list[str]]] = {}

_DEFAULT_GO_LIBRARIES = [
    'GO_Biological_Process_2026',
    'GO_Molecular_Function_2026',
    'GO_Cellular_Component_2026',
]

def _fetch_enrichr_library(library_name: str) -> dict[str, list[str]]:
    url = (f"https://maayanlab.cloud/Enrichr/geneSetLibrary"
           f"?mode=text&libraryName={urllib.parse.quote(library_name)}")
    with urllib.request.urlopen(url, timeout=60) as resp:
        raw = resp.read().decode("utf-8")
    lib: dict[str, list[str]] = {}
    for line in raw.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 3:
            lib[parts[0]] = [g for g in parts[2:] if g]
    return lib

def _resolve_go_terms(terms: list[str],
                      library_names: list[str]) -> dict[str, list[str]]:
    libs = []
    for name in library_names:
        if name not in _GO_LIBRARY_CACHE:
            _GO_LIBRARY_CACHE[name] = _fetch_enrichr_library(name)
        libs.append((name, _GO_LIBRARY_CACHE[name]))

    result: dict[str, list[str]] = {}
    for term in terms:
        term_lower = term.lower()
        matched_key = None
        matched_genes = None
        for lib_name, lib in libs:
            if term in lib:
                matched_key = term
                matched_genes = lib[term]
                break
            for key in lib:
                if term_lower in key.lower():
                    matched_key = key
                    matched_genes = lib[key]
                    break
            if matched_key is not None:
                break
        if matched_key is not None:
            result[matched_key] = matched_genes
        else:
            warnings.warn(f"compare_by_go: '{term}' not found in any of {library_names}")
    return result


def _diagonal_fisher(df: pd.DataFrame, mask: pd.Series,
                     background_mask: pd.Series | None = None) -> float:
    if background_mask is None:
        background_mask = pd.Series(True, index=df.index)
    sub  = df[mask]
    rest = df[background_mask & ~mask]
    sa = (sub['normalized_counts_y']  > sub['normalized_counts_x']).sum()
    sb = (sub['normalized_counts_x']  > sub['normalized_counts_y']).sum()
    ra = (rest['normalized_counts_y'] > rest['normalized_counts_x']).sum()
    rb = (rest['normalized_counts_x'] > rest['normalized_counts_y']).sum()
    if sa + sb == 0:
        return 1.0
    _, p = fisher_exact([[sa, sb], [ra, rb]])
    return p


def _fmt_p(p: float) -> str:
    return f"p={p:.2e}"


def _violin_plots(ax_col: list[matplotlib.axes.Axes], data: ad.AnnData, groups: list[str]) -> None:
    """Render a column of scanpy violin plots for a set of cell metadata groups.

    Args:
        ax_col: Ordered list of Matplotlib axes, one per group.
        data: AnnData object to plot.
        groups: obs column names to visualise, matched positionally to ax_col.
    """
    for i, group in enumerate(groups):
        sc.pl.violin(data, group, ax=ax_col[i], stripplot=False, show=False)
        ax_col[i].set_xticks([])


def plot_upsample(datasets: list[ad.AnnData], output_dir: Path, xlim: int = None, ylim: int = None) -> None:
    """Line plot of PreSeq upsampling predictions with confidence intervals.

    Lazily loads results from <output_dir>/<name>_yield.txt into each dataset's
    uns['pred_upsample'] if not already present. Datasets whose output file does
    not yet exist are skipped and reported by name; the plot is omitted entirely
    if no results are available.

    Args:
        datasets: List of AnnData objects with 'name' and 'title' in uns.
        output_dir: Directory where PreSeq wrote its yield files (same path
            passed to upsample()).
        xlim: Upper x-axis limit (total counts). Defaults to None (auto).
        ylim: Upper y-axis limit (projected unique transcripts). Defaults to None (auto).
    """
    ready = []
    for data in datasets:
        if 'pred_upsample' not in data.uns:
            fileout = output_dir / f"{data.uns['name']}_yield.txt"
            if not fileout.is_file():
                print(f"PreSeq output not ready for '{data.uns['title']}' — skipping.")
                continue
            with open(fileout, "r") as f:
                rows = []
                f.readline()
                for line in f:
                    rows.append(list(map(float, line.split())))
            data.uns['pred_upsample'] = np.array(rows)
        ready.append(data)

    if not ready:
        print("No PreSeq results available yet — rerun plot_upsample() when jobs complete.")
        return

    plt.figure(figsize=(16, 4))
    for data in ready:
        pred = data.uns['pred_upsample']
        plt.plot(pred[:, 0], pred[:, 1], label=data.uns['title'])
        plt.fill_between(pred[:, 0], pred[:, 2], pred[:, 3], alpha=0.5)

    plt.legend(loc='upper left')
    if xlim:
        plt.xlim(0, xlim)
    if ylim:
        plt.ylim(0, ylim)
    plt.xlabel('Total Number of Counts')
    plt.ylabel('Projected Number of Unique Transcripts')
    plt.show()


def plot_filtering_metrics(datasets: list[ad.AnnData]) -> None:
    """Stacked bar charts comparing aligned vs unaligned reads and filtered vs rejected counts.

    Args:
        datasets: List of AnnData objects with 'n_processed', 'n_aligned', 'n_raw_counts',
            and 'title' in uns.
    """
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    processed = np.array([data.uns['n_processed'] for data in datasets])
    aligned = np.array([data.uns['n_aligned'] for data in datasets])
    unique = np.array([data.uns['n_unique'] for data in datasets])
    counts = np.array([data.uns['n_raw_counts'] for data in datasets])
    filtered_counts = np.array([data.uns['n_raw_counts_filtered'] for data in datasets])
    labels = [data.uns['title'] for data in datasets]

    unmapped = processed - aligned
    mapped = aligned - unique

    rejected_counts = counts - filtered_counts

    axs[0].bar(labels, unique, label="Uniquely Mapped", color="yellow")
    axs[0].bar(labels, mapped, bottom=unique, label="Mapped", color="orange")
    axs[0].bar(labels, unmapped, bottom=mapped+unique, label="Unmapped", color="red")
    axs[0].set_ylabel('Number of Reads')
    axs[0].legend()

    axs[1].bar(labels, filtered_counts, label="High-Quality Counts", color='green')
    axs[1].bar(labels, rejected_counts, bottom=filtered_counts, label="Filtered Counts", color="blue")
    axs[1].set_ylabel('Number of Counts')
    axs[1].legend()

def plot_nascent_mature_ratios(datasets: list[ad.AnnData], figsize:Tuple[float,float]=(8,6)) -> None:
    """Stacked bar charts comparing the ratio of nascent, ambiguous, and mature reads

    Args:
        datasets: List of AnnData objects with 'percent_nascent', 'percent_ambiguous', 'percent_mature',
        and 'title' in uns.
        figsize: Figure size as (width, height)
    """

    fig, ax = plt.subplots(figsize=figsize)

    percent_nascent = np.array([data.uns['percent_nascent'] for data in datasets])
    percent_ambiguous = np.array([data.uns['percent_ambiguous'] for data in datasets])
    percent_mature = np.array([data.uns['percent_mature'] for data in datasets])
    labels = [data.uns['title'] for data in datasets]

    ax.bar(labels, percent_nascent, label="Nascent", color="yellow")
    ax.bar(labels, percent_ambiguous, bottom=percent_nascent, label="Ambiguous", color="orange")
    ax.bar(labels, percent_mature, bottom=percent_ambiguous+percent_nascent, label="Mature", color="red")
    ax.set_ylabel('Percent Counts')
    ax.legend()

def plot_filtering_metrics_star(datasets: list[ad.AnnData]) -> None:
    """Stacked bar charts comparing STARsolo's filtered vs unfiltered reads.

    Args:
        datasets: List of AnnData objects with 'star_uniquely_mapped', 'star_multimapped', 
        'star_unmapped' and 'title' in uns.
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    mapped = np.array([data.uns['star_uniquely_mapped'] for data in datasets])
    multimapped = np.array([data.uns['star_multimapped'] for data in datasets])
    unmapped = np.array([data.uns['star_unmapped'] for data in datasets])
    labels = [data.uns['title'] for data in datasets]

    ax.bar(labels, mapped, label="Uniquely Mapped", color="yellow")
    ax.bar(labels, multimapped, bottom=mapped, label="Multimapped", color="orange")
    ax.bar(labels, unmapped, bottom=mapped+multimapped, label="Unmapped", color="red")
    ax.set_ylabel('Number of Reads')
    ax.legend()

def plot_cell_metrics(datasets: list[ad.AnnData], groups: list[str], group_names: list[str], figsize: tuple[int, int] = (12, 18)) -> None:
    """Grid of violin plots for multiple cell metrics across datasets.

    Rows correspond to metrics, columns to datasets.

    Args:
        datasets: List of AnnData objects to plot.
        groups: obs column names for each row of the grid.
        group_names: Display labels for the y-axis of each row.
        figsize: Figure size as (width, height) in inches.
    """
    fig, ax = plt.subplots(len(groups), len(datasets), squeeze=False, figsize=(figsize), sharey='row')
    for i, data in enumerate(datasets):
        _violin_plots(ax[:, i], data, groups)

    for i, group in enumerate(group_names):
        ax[i, 0].set_ylabel(group)

    for i, data in enumerate(datasets):
        ax[0, i].set_title(data.uns['title'])
    plt.tight_layout()
    plt.show()


def plot_gene_family_percentages(
    datasets: list[ad.AnnData],
    gene_families: list[str],
    figsize: tuple[int, int] = None,
) -> None:
    """Grid of violin plots showing per-cell percentage of counts from specific gene families.

    Rows correspond to gene families, columns to datasets.

    Each entry in gene_families supports two operators:
      - "+" union:     "ATP+NDUF" matches genes starting with ATP or NDUF.
      - "/" exclusion: "COX/COX-" matches COX genes excluding those starting with COX-.
    Both operators can be combined: "COX+NDUF/COX-" matches (COX or NDUF) minus COX-.
    Mouse (capitalized) and barnyard (HUMAN_/MOUSE_) variants are derived automatically
    for every prefix on both sides of /.

    Args:
        datasets: List of AnnData objects to plot.
        gene_families: Gene family specifiers (e.g. ["MRP", "ATP+NDUF", "COX/COX-"]).
        figsize: Figure size as (width, height). Defaults to (3*n_datasets, 3*n_families).
    """
    def _variants(prefix: str) -> tuple[str, ...]:
        prefix = prefix.lower()
        return (prefix, f"human_{prefix}", f"mouse_{prefix}")

    def _resolve_mask(var_index: pd.Index, family: str) -> pd.Index:
        include_str, *exclude_parts = family.split('/')
        include_prefixes = tuple(v for p in include_str.split('+') for v in _variants(p))
        include_mask = var_index.str.lower().str.startswith(include_prefixes)
        if exclude_parts:
            exclude_prefixes = tuple(v for p in exclude_parts[0].split('+') for v in _variants(p))
            include_mask &= ~var_index.str.lower().str.startswith(exclude_prefixes)
        return var_index[include_mask]

    if figsize is None:
        figsize = (3 * len(datasets), 3 * len(gene_families))

    fig, ax = plt.subplots(len(gene_families), len(datasets), figsize=figsize, sharey='row')

    for j, data in enumerate(datasets):
        for i, family in enumerate(gene_families):
            mask = _resolve_mask(data.var.index, family)

            if len(mask) == 0:
                ax[i, j].text(0.5, 0.5, 'no genes', ha='center', va='center',
                              transform=ax[i, j].transAxes, fontsize=8)
                ax[i, j].set_xticks([])
            else:
                pct = np.nan_to_num(
                    data[:, mask].X.toarray().sum(axis=1) /
                    np.array(data.obs['n_counts']) * 100
                )
                ax[i, j].violinplot(pct, showmedians=True)
                ax[i, j].set_xticks([])

        ax[0, j].set_title(data.uns['title'])

    for i, family in enumerate(gene_families):
        ax[i, 0].set_ylabel(f"{family} %")

    plt.tight_layout()
    plt.show()


def plot_gene_metrics(datasets: list[ad.AnnData], sample_size: int = 1000000) -> None:
    """Violin plots of count-weighted gene length and GC content distributions.

    Args:
        datasets: List of AnnData objects whose var contains 'gene_length' and
            'gc_content' columns.
        sample_size: Number of random draws used to build each violin distribution.
    """
    fig, ax = plt.subplots(2, len(datasets), figsize=(12, 10), sharey='row')

    for i, data in enumerate(datasets):
        gene_lengths = data.var['gene_length'].values
        gc_content = data.var['gc_content'].values
        gene_counts = np.array(data.X.sum(axis=0)).flatten()

        length_mask = (gene_lengths > 1) & (gene_counts > 0)
        gc_mask = (gene_counts > 0) & (gc_content > 0)

        length_weights = gene_counts[length_mask] / gene_counts[length_mask].sum()
        gc_weights = gene_counts[gc_mask] / gene_counts[gc_mask].sum()

        lengths_sampled = np.random.choice(gene_lengths[length_mask], size=sample_size, p=length_weights)
        gc_sampled = np.random.choice(gc_content[gc_mask], size=sample_size, p=gc_weights)

        ax[0, i].violinplot(np.log10(lengths_sampled), showextrema=False)
        ax[0, i].set_title(data.uns['title'])
        ax[0, i].set_xticks([])

        ax[1, i].violinplot(gc_sampled, showextrema=False)
        ax[1, i].set_xticks([])

    ax[0, 0].set_ylabel('Gene Length (log10, weighted by counts)')
    ax[1, 0].set_ylabel('Percent GC Content')

    plt.tight_layout()
    plt.show()


def marker_genes(datasets: list[ad.AnnData], markers: list[str], figsize:Tuple[float,float] = (25,5), plot_title: str = "Marker Genes") -> None:
    """Grid of violin plots of marker gene expression as a percent of total counts per cell for each dataset in `datasets`

    Args:
        datasets: List of AnnData objects with 'title' in uns and 'n_counts' in obs
        markers: List of gene names to plot (must be present in data.var_names).
        figsize: Figure size as (width, height)
        plot_title: The string to be appended to the technology name in the title
            of each plot
    """

    fig, axs = plt.subplots(1, len(datasets), figsize=figsize, sharey=True)

    for ax, data in zip(axs, datasets):
        gene_dist = []
        for gene in markers:
            gene_dist.append(np.nan_to_num(data[:, gene].X.toarray().transpose()[0] / np.array(data.obs['n_counts'].tolist()) * 100))

        ax.violinplot(gene_dist, showmeans=True)
        ax.set_xticks(np.arange(1, len(markers) + 1), markers)
        ax.set_ylabel("")
        ax.set_title(data.uns['title'] + " " +  plot_title)

    axs[0].set_ylabel("% counts")

    plt.tight_layout()
    plt.show()


def collapsed_marker_genes(datasets: list[ad.AnnData], marker_prefixes: list[str], figsize:Tuple[float,float] = (25,5)) -> None:
    """Grid of violin plots of collapsed marker gene family expression per cell
    for datasets in the list 'datasets'

    Sums counts across all genes sharing a common prefix, useful when individual
    gene distinctions are not important (e.g. ribosomal protein families).

    Args:
        datasets: List of AnnData objects with 'title' in uns and 'n_counts' in obs
        marker_prefixes: Gene name prefixes to collapse. Each prefix becomes one
            violin in the plot.
        figsize: Figure size as (width, height)
    """
    fig, axs = plt.subplots(1, len(datasets), figsize=figsize, sharey=True)

    for ax, data in zip(axs, datasets):
        genes = data.var
        gene_dist = []
        for prefix in marker_prefixes:
            mask = genes[genes.index.str.startswith(prefix)].index
            gene_dist.append(np.nan_to_num(data[:, mask].X.toarray().sum(axis=1).transpose() / np.array(data.obs['n_counts'].tolist()) * 100))

        ax.violinplot(gene_dist, showmeans=True)
        ax.set_xticks(np.arange(1, len(marker_prefixes) + 1), marker_prefixes)
        ax.set_ylabel("")
        ax.set_title(data.uns['title'] + " Thymus Marker Genes")
    axs[0].set_ylabel("% counts")

    plt.tight_layout()
    plt.show()


def top_gene_cell_expression(datasets: list[ad.AnnData], figsize: Tuple[float,float] = (25,5)) -> None:
    """Grid of violin plots of the top 10 genes by number of cells expressing them.

    Args:
        datasets: List of AnnData objects with 'title' in uns, 'n_counts' in obs, 
            and'n_cells' in var.
        figsize: Figure size as (width, height)
    """
    fig, axs = plt.subplots(1, len(datasets), figsize=figsize, sharey=True)

    for ax, data in zip(axs, datasets):
        genes = data.var
        top_genes = genes.sort_values(by='n_cells', ascending=False).head(10).index

        gene_dist = []
        for gene in top_genes:
            gene_dist.append(np.nan_to_num(data[:, gene].X.toarray().transpose()[0] / np.array(data.obs['n_counts'].tolist()) * 100))

        ax.violinplot(gene_dist, showmeans=True)
        ax.set_xticks(np.arange(1, len(top_genes) + 1), top_genes)
        ax.set_ylabel("")
        ax.set_title(data.uns['title'] + " Top 10 Genes by Cell Expression")
        

    axs[0].set_ylabel("% counts")

    plt.tight_layout()
    plt.show()


def top_gene_counts(datasets: list[ad.AnnData], figsize: Tuple[float,float] = (25,5)) -> None:
    """Grid of violin plots of the top 10 genes by total percent counts.

    Args:
        datasets: List of AnnData objects with 'title' in uns, 'n_counts' in obs, and
            'percent_counts' in var.
        figsize: Figure size as (width, height)
    """

    fig, axs = plt.subplots(1, len(datasets), figsize=figsize, sharey=True)

    for ax, data in zip(axs, datasets):
        genes = data.var
        top_genes = genes.sort_values(by='percent_counts', ascending=False).head(10).index

        gene_dist = []
        for gene in top_genes:
            gene_dist.append(np.nan_to_num(data[:, gene].X.toarray().transpose()[0] / np.array(data.obs['n_counts'].tolist()) * 100))

        ax.violinplot(gene_dist, showmeans=True)
        ax.set_xticks(np.arange(1, len(top_genes) + 1), top_genes)
        ax.set_ylabel("")
        ax.set_title(data.uns['title'] + " Top 10 Genes by Total Counts")
        

    axs[0].set_ylabel("% counts")

    plt.tight_layout()
    plt.show()


def scatter_genes(ax: matplotlib.axes.Axes, shared_data: pd.DataFrame, data_x: ad.AnnData, data_y: ad.AnnData,
                  c_column: str, xlim: float = None, ylim: float = None,
                  norm: matplotlib.colors.Normalize = None) -> matplotlib.collections.PathCollection:
    """Scatter plot of log-normalized gene counts comparing two datasets, colored by a metric.

    Args:
        ax: Matplotlib axes to draw on.
        shared_data: Output of compare_genes containing 'normalized_counts_x' and
            'normalized_counts_y' columns.
        data_x: AnnData object for the x-axis dataset (used for axis label).
        data_y: AnnData object for the y-axis dataset (used for axis label).
        c_column: Column in shared_data used to color points.
        xlim: Upper x-axis limit. Defaults to None (auto).
        ylim: Upper y-axis limit. Defaults to None (auto).
        norm: Matplotlib normalisation applied to the color mapping.

    Returns:
        The PathCollection returned by ax.scatter, for use with fig.colorbar.
    """
    x_counts = shared_data['normalized_counts_x']
    y_counts = shared_data['normalized_counts_y']

    plot = ax.scatter(x_counts,
                      y_counts,
                      norm=norm,
                      alpha=0.5,
                      s=50,
                      c=shared_data[c_column],
                      cmap='viridis')

    ax.set_xlabel(data_x.uns['title'] + ' Log-normalized Counts')
    ax.set_ylabel(data_y.uns['title'] + ' Log-normalized Counts')
    if xlim:
        ax.set_xlim(0, xlim)
    if ylim:
        ax.set_ylim(0, ylim)

    return plot


def cat_scatter_genes(ax: matplotlib.axes.Axes, shared_data: pd.DataFrame, data_x: ad.AnnData, data_y: ad.AnnData,
                      color: str, label: str = None, xlim: float = None,
                      ylim: float = None) -> matplotlib.collections.PathCollection:
    """Scatter plot of log-normalized gene counts with a single categorical color.

    Args:
        ax: Matplotlib axes to draw on.
        shared_data: Output of compare_genes containing 'normalized_counts_x' and
            'normalized_counts_y' columns.
        data_x: AnnData object for the x-axis dataset (used for axis label).
        data_y: AnnData object for the y-axis dataset (used for axis label).
        color: Matplotlib color string applied to all points.
        label: Legend label for this scatter series. Defaults to None.
        xlim: Upper x-axis limit. Defaults to None (auto).
        ylim: Upper y-axis limit. Defaults to None (auto).

    Returns:
        The PathCollection returned by ax.scatter.
    """
    x_counts = shared_data['normalized_counts_x']
    y_counts = shared_data['normalized_counts_y']

    plot = ax.scatter(x_counts,
                      y_counts,
                      s=50,
                      alpha=0.5,
                      color=color,
                      label=label)

    ax.set_xlabel(data_x.uns['title'] + ' Log-normalized Counts')
    ax.set_ylabel(data_y.uns['title'] + ' Log-normalized Counts')
    if xlim:
        ax.set_xlim(0, xlim)
    if ylim:
        ax.set_ylim(0, ylim)

    return plot


def show_correlation(ax: matplotlib.axes.Axes, shared_data: pd.DataFrame) -> None:
    """Annotate a scatter plot with Pearson, Spearman, and CCC correlation coefficients.

    Also draws a y=x reference line.

    Args:
        ax: Matplotlib axes to annotate.
        shared_data: DataFrame with 'normalized_counts_x' and 'normalized_counts_y' columns.
    """
    x_counts = shared_data['normalized_counts_x']
    y_counts = shared_data['normalized_counts_y']

    pearson_r = pearsonr(x_counts, y_counts).correlation
    spearman_r = spearmanr(x_counts, y_counts).correlation

    x_mean = x_counts.mean()
    y_mean = y_counts.mean()
    x_std = x_counts.std()
    y_std = y_counts.std()
    CCC_r = 2 * pearson_r * x_std * y_std / (x_std**2 + y_std**2 + (x_mean - y_mean)**2)

    lim = max(x_counts.max(), y_counts.max())
    ax.plot([0, lim], [0, lim], color='black', linestyle='--', linewidth=1.5)

    textstr = '\n'.join((
        r'$\mathrm{r}=%.2f$' % (pearson_r,),
        r'$\rho=%.2f$' % (spearman_r,),
        r'$\rho_C=%.2f$' % (CCC_r,),))
    ax.text(0.05, 0.85, textstr, transform=ax.transAxes, fontsize=12)


def _label_genes(ax: matplotlib.axes.Axes, df: pd.DataFrame,
                 n_cooks: int, n_log_ratio: int, min_expr: float = 0.0) -> None:
    eps = 1e-10
    expressed = df[(df['normalized_counts_x'] >= min_expr) | (df['normalized_counts_y'] >= min_expr)]
    log_ratio = np.abs(np.log((expressed['normalized_counts_x'] + eps) / (expressed['normalized_counts_y'] + eps)))
    top_cooks_idx = df.nlargest(n_cooks, 'cooks_distance').index
    top_ratio_idx = log_ratio.nlargest(n_log_ratio).index
    for idx in top_cooks_idx.union(top_ratio_idx):
        row = df.loc[idx]
        ax.annotate(row['gene_name'], (row['normalized_counts_x'], row['normalized_counts_y']),
                    fontsize=7, ha='left', va='bottom', clip_on=True)


def comparison_plotter(compare_dfs: list[pd.DataFrame], comparisons: list[tuple],
                       norm: matplotlib.colors.Normalize,
                       metric: str, metric_name: str,
                       n_cooks: int = 5, n_log_ratio: int = 5, min_expr: float = 0.0) -> None:
    """Four-panel scatter plot comparing log-normalized gene counts colored by a continuous metric.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        norm: Matplotlib normalisation applied to the color mapping.
        metric: Column in each DataFrame used to color points.
        metric_name: Label shown on the colorbar.
    """
    lim = max(max(df['normalized_counts_x'].max(), df['normalized_counts_y'].max())
              for df in compare_dfs)
    fig, axs = plt.subplots(1, 4, figsize=(25, 5))

    for ax, df, pair in zip(axs, compare_dfs, comparisons):
        plot = scatter_genes(ax, df, pair[0], pair[1], metric, xlim=lim, ylim=lim, norm=norm)
        show_correlation(ax, df)
        _label_genes(ax, df, n_cooks, n_log_ratio, min_expr)

    fig.colorbar(plot, label=metric_name)

    plt.tight_layout()
    plt.show()


def compare_by_density(compare_dfs: list[pd.DataFrame],
                       comparisons: list[tuple],
                       n_cooks: int = 5, n_log_ratio: int = 5, min_expr: float = 0.0) -> None:
    """Four-panel comparison scatter plot colored by point density.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log ratio to label.
        min_expr: Genes where both normalized_counts are below this value are excluded
            from the log-ratio ranking.
    """
    c_values = []
    for df in compare_dfs:
        c_values.extend(df['point_density'].tolist())
    norm = LogNorm(min(c_values), max(c_values))

    comparison_plotter(compare_dfs, comparisons, norm, 'point_density', 'Density',
                       n_cooks, n_log_ratio, min_expr)


def compare_by_cooks(compare_dfs: list[pd.DataFrame], comparisons: list[tuple],
                     n_cooks: int = 5, n_log_ratio: int = 5, min_expr: float = 0.0) -> None:
    """Four-panel comparison scatter plot colored by Cook's distance.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log ratio to label.
        min_expr: Genes where both normalized_counts are below this value are excluded
            from the log-ratio ranking.
    """
    c_values = []
    for df in compare_dfs:
        c_values.extend(df['cooks_distance'].tolist())
    norm = LogNorm(min(c_values), max(c_values))

    comparison_plotter(compare_dfs, comparisons, norm, 'cooks_distance', "Cook's Distance",
                       n_cooks, n_log_ratio, min_expr)


def compare_by_length(compare_dfs: list[pd.DataFrame], comparisons: list[tuple],
                      n_cooks: int = 5, n_log_ratio: int = 5, min_expr: float = 0.0) -> None:
    """Four-panel comparison scatter plot colored by gene length.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log ratio to label.
        min_expr: Genes where both normalized_counts are below this value are excluded
            from the log-ratio ranking.
    """
    c_values = []
    for df in compare_dfs:
        df['gene_length'] = df['gene_length'] + 1
        c_values.extend(df['gene_length'].tolist())
    norm = LogNorm(min(c_values), max(c_values))

    comparison_plotter(compare_dfs, comparisons, norm, 'gene_length', 'Gene Length',
                       n_cooks, n_log_ratio, min_expr)


def compare_by_gc(compare_dfs: list[pd.DataFrame], comparisons: list[tuple],
                  n_cooks: int = 5, n_log_ratio: int = 5, min_expr: float = 0.0) -> None:
    """Four-panel comparison scatter plot colored by GC content.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log ratio to label.
        min_expr: Genes where both normalized_counts are below this value are excluded
            from the log-ratio ranking.
    """
    c_values = []
    for df in compare_dfs:
        c_values.extend(df['gc_content'].tolist())
    norm = Normalize(min(c_values), max(c_values))

    comparison_plotter(compare_dfs, comparisons, norm, 'gc_content', 'Percent GC Content',
                       n_cooks, n_log_ratio, min_expr)


def compare_by_type(compare_dfs: list[pd.DataFrame], comparisons: list[tuple],
                    n_cooks: int = 5, n_log_ratio: int = 5, min_expr: float = 0.0) -> None:
    """Four-panel comparison scatter plot with categorical colors by gene biotype.

    Colors: protein-coding (yellow), mtRNA (blue), rRNA (red), lncRNA (green),
    unspecified non-coding (black).

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log ratio to label.
        min_expr: Genes where both normalized_counts are below this value are excluded
            from the log-ratio ranking.
    """
    lim = max(max(df['normalized_counts_x'].max(), df['normalized_counts_y'].max())
              for df in compare_dfs)
    fig, axs = plt.subplots(1, 4, figsize=(25, 5))

    for ax, df, pair in zip(axs, compare_dfs, comparisons):
        cat_scatter_genes(ax, df[~(df['is_mito'] | df['is_ribo'] | df['is_lnc'] | df['is_pseudo'])],
                          pair[0], pair[1], 'black',
                          label='unspecified non-coding', xlim=lim, ylim=lim)
        for col, c, label in zip(['is_pc', 'is_pseudo', 'is_mito', 'is_ribo', 'is_lnc'],
                                  ['yellow', 'cyan', 'blue', 'red', 'green'],
                                  ['protein coding', 'pseudogene', 'mtRNA', 'rRNA', 'lncRNA']):
            cat_scatter_genes(ax, df[df[col]], pair[0], pair[1],
                              c, label=label, xlim=lim, ylim=lim)

        show_correlation(ax, df)
        _label_genes(ax, df, n_cooks, n_log_ratio, min_expr)

    axs[0].legend()

    plt.tight_layout()
    plt.show()


def compare_by_pathway(compare_dfs: list[pd.DataFrame], comparisons: list[tuple],
                       n_cooks: int = 10, n_log_ratio: int = 10, min_expr: float = 0.0) -> None:
    """Four-panel comparison scatter plot highlighting ribosomal and OXPHOS genes.

    All genes are drawn in light grey; ribosomal (RP/MRP) genes are overlaid in
    red and OXPHOS genes (SDH, UQCR, COX, NDUF, ATP, and their mitochondria-encoded
    equivalents) in blue.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log ratio to label.
        min_expr: Genes where both normalized_counts are below this value are excluded
            from the log-ratio ranking.
    """
    lim = max(max(df['normalized_counts_x'].max(), df['normalized_counts_y'].max())
              for df in compare_dfs)
    fig, axs = plt.subplots(1, 4, figsize=(25, 5))

    for ax, df, pair in zip(axs, compare_dfs, comparisons):
        cat_scatter_genes(ax, df, pair[0], pair[1], 'lightgrey', xlim=lim, ylim=lim)
        cat_scatter_genes(ax, df[df['is_ribo']], pair[0], pair[1],
                          'red', label='ribosomal (RP/MRP)', xlim=lim, ylim=lim)
        cat_scatter_genes(ax, df[df['is_oxphos']], pair[0], pair[1],
                          'blue', label='OXPHOS', xlim=lim, ylim=lim)

        show_correlation(ax, df)
        _label_genes(ax, df, n_cooks, n_log_ratio, min_expr)

    axs[0].legend()

    plt.tight_layout()
    plt.show()


def compare_by_enrichment(compare_dfs: list[pd.DataFrame], comparisons: list[tuple],
                           enrichment_path: str | Path,
                           gene_sets: list[str] | None = None,
                           terms: list[str] | None = None,
                           n_cooks: int = 10, n_log_ratio: int = 10, min_expr: float = 0.0,
                           label_enriched: bool = False) -> None:
    """Four-panel comparison scatter plot highlighting genes from GO enrichment terms.

    Each selected term is drawn in a distinct colour (tab10 palette cycling); all
    other genes are drawn in light grey. When a gene belongs to multiple terms it
    takes the colour of the last term rendered. The enrichment CSV must be in the
    format produced by gseapy.enrichr (columns: Gene_set, Term, Genes with semicolons).

    Gene names are normalised before matching so that mouse (e.g. Rpl4) and barnyard
    (e.g. HUMAN_RPL4, MOUSE_Rpl4) datasets correctly map to HGNC symbols in the CSV.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        enrichment_path: Path to an enrichment_results.csv from gseapy.enrichr.
        gene_sets: If provided, restrict to rows whose Gene_set is in this list.
        terms: If provided, restrict to rows whose Term is in this list.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log ratio to label.
        min_expr: Genes where both normalized_counts are below this value are excluded
            from the log-ratio ranking.
        label_enriched: If True, annotate every gene belonging to any selected term
            with its gene name, in addition to the top Cook's/log-ratio labels.
    """
    orthologs = _load_orthologs()
    _norm = lambda name: _normalize_gene_name(name, orthologs)
    # background = genes whose normalized name is a known human HGNC symbol
    human_hgnc = set(orthologs.values())

    enr_df = pd.read_csv(enrichment_path, index_col=0)
    if gene_sets is not None:
        enr_df = enr_df[enr_df['Gene_set'].isin(gene_sets)]
    if terms is not None:
        enr_df = enr_df[enr_df['Term'].isin(terms)]

    per_term_genes: dict[str, set[str]] = {}
    for _, row in enr_df.iterrows():
        t = row['Term']
        if t not in per_term_genes:
            per_term_genes[t] = set()
        per_term_genes[t].update(_norm(g) for g in row['Genes'].split(';'))
    selected_terms = list(per_term_genes.keys())

    available = set()
    for df in compare_dfs:
        available.update(df['gene_name'].apply(_norm))
    for term, gene_set in per_term_genes.items():
        _warn_unmatched(gene_set, available, "compare_by_enrichment", f'"{term}"', "data")

    cmap = matplotlib.colormaps.get_cmap('tab10')
    colors = [cmap(i % 10) for i in range(len(selected_terms))]

    all_enriched: set[str] = set().union(*per_term_genes.values()) if per_term_genes else set()

    lim = max(max(df['normalized_counts_x'].max(), df['normalized_counts_y'].max())
              for df in compare_dfs)
    fig, axs = plt.subplots(1, 4, figsize=(25, 5))

    for ax, df, pair in zip(axs, compare_dfs, comparisons):
        x_label = pair[0].uns['title']
        y_label = pair[1].uns['title']
        norm_names = df['gene_name'].apply(_norm)
        background_mask = norm_names.isin(human_hgnc)
        cat_scatter_genes(ax, df, pair[0], pair[1], 'lightgrey', xlim=lim, ylim=lim)

        for term, color in zip(selected_terms, colors):
            mask = norm_names.isin(per_term_genes[term])
            p = _diagonal_fisher(df, mask, background_mask)
            label = f"{term} ({_fmt_p(p)})"
            cat_scatter_genes(ax, df[mask], pair[0], pair[1],
                              color, label=label, xlim=lim, ylim=lim)
            sub = df[mask]
            n_above = (sub['normalized_counts_y'] > sub['normalized_counts_x']).sum()
            n_below = (sub['normalized_counts_x'] > sub['normalized_counts_y']).sum()
            n_total = len(sub)
            print(f"  {x_label} vs {y_label} | {term}: "
                  f"{n_below}/{n_total} toward {x_label}, "
                  f"{n_above}/{n_total} toward {y_label}, {_fmt_p(p)}")

        ax.plot([0, lim], [0, lim], color='black', linestyle='--', linewidth=1.5)
        _label_genes(ax, df, n_cooks, n_log_ratio, min_expr)

        if label_enriched:
            for idx in df[norm_names.isin(all_enriched)].index:
                row = df.loc[idx]
                ax.annotate(row['gene_name'],
                            (row['normalized_counts_x'], row['normalized_counts_y']),
                            fontsize=7, ha='left', va='bottom', clip_on=True)

        ax.legend(fontsize=7)

    plt.tight_layout()
    plt.show()

def compare_by_de(compare_dfs: list[pd.DataFrame], comparisons: list[tuple],
                   de_results: pd.DataFrame, de_comparing: tuple[str, str],
                   fdr_thresh: float = 0.05,
                   n_cooks: int = 10, n_log_ratio: int = 10, min_expr: float = 0.0,
                   label_de: bool = False) -> None:
    """Overlay differentially-expressed genes onto the four-panel comparison scatter.

    Takes the significant genes (FDR < fdr_thresh) from an edgepy/edgeR result and,
    on each comparison panel, colors the genes higher in each side and reports a
    diagonal Fisher test of whether that set sits preferentially off the identity line
    (i.e. is enriched in one tech). The background is all genes edgeR tested.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        de_results: edgepy result table with 'gene_name', 'logFC', and 'FDR' columns.
        de_comparing: (side1, side2) labels; logFC < 0 is higher in side1, > 0 in side2.
        fdr_thresh: Significance cutoff on FDR for calling a gene differential.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log ratio to label.
        min_expr: Genes where both normalized_counts are below this value are excluded
            from the log-ratio ranking.
        label_de: If True, annotate every differential gene with its name.
    """
    orthologs = _load_orthologs()
    _norm = lambda name: _normalize_gene_name(name, orthologs)

    de_side1, de_side2 = de_comparing
    sig = de_results[de_results['FDR'] < fdr_thresh].copy()
    sig['_norm'] = sig['gene_name'].apply(_norm)
    # background = all genes tested by edgeR (not just significant ones)
    tested_genes = set(de_results['gene_name'].apply(_norm))

    genes_side1 = set(sig.loc[sig['logFC'] < 0, '_norm'])
    genes_side2 = set(sig.loc[sig['logFC'] > 0, '_norm'])

    lim = max(max(df['normalized_counts_x'].max(), df['normalized_counts_y'].max())
              for df in compare_dfs)
    fig, axs = plt.subplots(1, 4, figsize=(25, 5))

    for ax, df, pair in zip(axs, compare_dfs, comparisons):
        x_label = pair[0].uns['title']
        y_label = pair[1].uns['title']
        norm_names = df['gene_name'].apply(_norm)
        background_mask = norm_names.isin(tested_genes)

        cat_scatter_genes(ax, df, pair[0], pair[1], 'lightgrey', xlim=lim, ylim=lim)

        for gene_set, color, base_label in [
            (genes_side1, '#1f77b4', f'Higher in {de_side1} (n={len(genes_side1)})'),
            (genes_side2, '#d62728', f'Higher in {de_side2} (n={len(genes_side2)})'),
        ]:
            mask = norm_names.isin(gene_set)
            p = _diagonal_fisher(df, mask, background_mask)
            label = f"{base_label}"
            cat_scatter_genes(ax, df[mask], pair[0], pair[1],
                              color, label=label, xlim=lim, ylim=lim)
            sub = df[mask]
            n_above = (sub['normalized_counts_y'] > sub['normalized_counts_x']).sum()
            n_below = (sub['normalized_counts_x'] > sub['normalized_counts_y']).sum()
            n_total = len(sub)
            print(f"  {x_label} vs {y_label} | {base_label}: "
                  f"{n_below}/{n_total} ({(n_below/n_total*100):.2f}%) toward {x_label}, "
                  f"{n_above}/{n_total} ({(n_above/n_total*100):.2f}%) toward {y_label}, {_fmt_p(p)}")

        ax.plot([0, lim], [0, lim], color='black', linestyle='--', linewidth=1.5)
        _label_genes(ax, df, n_cooks, n_log_ratio, min_expr)

        if label_de:
            all_de = genes_side1 | genes_side2
            for idx in df[norm_names.isin(all_de)].index:
                row = df.loc[idx]
                ax.annotate(row['gene_name'],
                            (row['normalized_counts_x'], row['normalized_counts_y']),
                            fontsize=7, ha='left', va='bottom', clip_on=True)

        ax.legend(fontsize=7)
    plt.tight_layout()
    plt.show()


def compare_by_go(compare_dfs: list[pd.DataFrame], comparisons: list[tuple],
                   terms: list[str],
                   go_library: str | list[str] | None = None,
                   n_cooks: int = 10, n_log_ratio: int = 10, min_expr: float = 0.0,
                   label_genes: bool = False) -> None:
    """Overlay genes annotated to one or more GO terms onto the comparison scatter.

    Resolves each requested term to its gene set from the Enrichr GO libraries, then on
    every comparison panel colors those genes and reports a diagonal Fisher test of
    whether the set sits preferentially off the identity line. The background is the
    union of all genes in the loaded GO libraries (the annotatable universe).

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        terms: GO term names or GO IDs to look up and overlay, one color each.
        go_library: Enrichr GO library name(s) to search; defaults to the module's
            standard GO libraries when None.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log ratio to label.
        min_expr: Genes where both normalized_counts are below this value are excluded
            from the log-ratio ranking.
        label_genes: If True, annotate the highlighted term genes with their names.

    Raises:
        ValueError: If none of the requested terms match a library entry.
    """
    orthologs = _load_orthologs()
    _norm = lambda name: _normalize_gene_name(name, orthologs)

    library_names = (
        _DEFAULT_GO_LIBRARIES if go_library is None
        else ([go_library] if isinstance(go_library, str) else go_library)
    )
    per_term_genes: dict[str, list[str]] = _resolve_go_terms(terms, library_names)
    if not per_term_genes:
        raise ValueError("No GO terms matched — check term names or GO IDs.")

    all_term_genes = set().union(*[set(g) for g in per_term_genes.values()])

    # background = all genes present in any loaded GO library (the annotatable universe)
    bg_genes: set[str] = set()
    for lib in _GO_LIBRARY_CACHE.values():
        for gene_list in lib.values():
            bg_genes.update(gene_list)

    cmap = matplotlib.colormaps.get_cmap('tab10')
    colors = [cmap(i % 10) for i in range(len(per_term_genes))]

    lim = max(max(df['normalized_counts_x'].max(), df['normalized_counts_y'].max())
              for df in compare_dfs)
    fig, axs = plt.subplots(1, 4, figsize=(25, 5))

    for ax, df, pair in zip(axs, compare_dfs, comparisons):
        x_label = pair[0].uns['title']
        y_label = pair[1].uns['title']
        norm_names = df['gene_name'].apply(_norm)
        background_mask = norm_names.isin(bg_genes)

        cat_scatter_genes(ax, df, pair[0], pair[1], 'lightgrey', xlim=lim, ylim=lim)

        for (term, gene_list), color in zip(per_term_genes.items(), colors):
            gene_set = {_norm(g) for g in gene_list}
            mask = norm_names.isin(gene_set)
            short_label = term.split(' (GO:')[0]
            p = _diagonal_fisher(df, mask, background_mask)
            label = f"{short_label} ({_fmt_p(p)})"
            cat_scatter_genes(ax, df[mask], pair[0], pair[1],
                              color, label=label, xlim=lim, ylim=lim)
            sub = df[mask]
            n_above = (sub['normalized_counts_y'] > sub['normalized_counts_x']).sum()
            n_below = (sub['normalized_counts_x'] > sub['normalized_counts_y']).sum()
            n_total = len(sub)
            print(f"  {x_label} vs {y_label} | {short_label}: "
                  f"{n_below}/{n_total} ({(n_below/n_total*100):.2f}%) toward {x_label}, "
                  f"{n_above}/{n_total} ({(n_above/n_total*100):.2f}%) toward {y_label}, {_fmt_p(p)}")

        ax.plot([0, lim], [0, lim], color='black', linestyle='--', linewidth=1.5)
        _label_genes(ax, df, n_cooks, n_log_ratio, min_expr)

        if label_genes:
            all_norm = {_norm(g) for g in all_term_genes}
            for idx in df[norm_names.isin(all_norm)].index:
                row = df.loc[idx]
                ax.annotate(row['gene_name'],
                            (row['normalized_counts_x'], row['normalized_counts_y']),
                            fontsize=7, ha='left', va='bottom', clip_on=True)

        ax.legend(fontsize=7)
    plt.tight_layout()
    plt.show()


def compare(datasets: list[ad.AnnData],
            comparison_axis: str = "normalized_counts",
            n_cooks: int = 5,
            n_log_ratio: int = 10,
            min_expr: float = 0.0,
) -> tuple[list[str], list[pd.DataFrame]]:
    """Generate all pairwise gene-count comparison plots across four sequencing methods.

    Computes compare_genes for every pair and renders density, Cook's, length, GC,
    and biotype comparison plots. Axis limits are determined automatically from the data.

    Comparison pairs: polyT vs randO, 10X vs polyT, 10X vs randO, 10X vs Parse.

    Args:
        datasets: List of four AnnData objects in order [10X, polyT, randO, Parse].
        comparison_axis: var column to compare (default: 'normalized_counts').
        n_cooks: Number of top genes by Cook's distance to label on each plot.
        n_log_ratio: Number of top genes by absolute log ratio to label on each plot.
        min_expr: Genes where both normalized_counts are below this value are excluded
            from the log-ratio ranking.

    Returns:
        Tuple of (compare_names, compare_dfs) where compare_names is a list of
        pair label strings and compare_dfs is the corresponding list of DataFrames
        from compare_genes.
    """
    data_10x = datasets[0]
    data_polyT = datasets[1]
    data_randO = datasets[2]
    data_parse = datasets[3]

    comparisons = [(data_polyT, data_randO),
                   (data_10x, data_polyT),
                   (data_10x, data_randO),
                   (data_10x, data_parse)]
    compare_names = ['polyT_randO',
                     '10x_polyT',
                     '10x_randO',
                     '10x_parse']

    compare_dfs = []
    for pair in comparisons:
        compare_dfs.append(processing.compare_genes(pair[0], pair[1], comparison_axis))

    compare_by_density(compare_dfs, comparisons, n_cooks, n_log_ratio, min_expr)
    compare_by_cooks(compare_dfs, comparisons, n_cooks, n_log_ratio, min_expr)
    compare_by_length(compare_dfs, comparisons, n_cooks, n_log_ratio, min_expr)
    compare_by_gc(compare_dfs, comparisons, n_cooks, n_log_ratio, min_expr)
    compare_by_type(compare_dfs, comparisons, n_cooks, n_log_ratio, min_expr)
    compare_by_pathway(compare_dfs, comparisons, n_cooks, n_log_ratio, min_expr)

    return compare_names, compare_dfs
