from pathlib import Path
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.axes
import matplotlib.collections
import numpy as np
import scanpy as sc
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from matplotlib.colors import LogNorm, Normalize
from upsetty import Upset
from typing import Tuple
from . import processing


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
    """Scatter plot of gene percent counts comparing two datasets, colored by a metric.

    Args:
        ax: Matplotlib axes to draw on.
        shared_data: Output of compare_genes containing 'percent_counts_x' and
            'percent_counts_y' columns.
        data_x: AnnData object for the x-axis dataset (used for axis label).
        data_y: AnnData object for the y-axis dataset (used for axis label).
        c_column: Column in shared_data used to color points.
        xlim: Upper x-axis limit. Defaults to None (auto).
        ylim: Upper y-axis limit. Defaults to None (auto).
        norm: Matplotlib normalisation applied to the color mapping.

    Returns:
        The PathCollection returned by ax.scatter, for use with fig.colorbar.
    """
    x_percent = shared_data['percent_counts_x']
    y_percent = shared_data['percent_counts_y']

    plot = ax.scatter(x_percent,
                      y_percent,
                      norm=norm,
                      alpha=0.5,
                      s=50,
                      c=shared_data[c_column],
                      cmap='viridis')

    ax.set_xlabel(data_x.uns['title'] + ' Gene Percent Count')
    ax.set_ylabel(data_y.uns['title'] + ' Gene Percent Count')
    if xlim:
        ax.set_xlim(0, xlim)
    if ylim:
        ax.set_ylim(0, ylim)

    return plot


def cat_scatter_genes(ax: matplotlib.axes.Axes, shared_data: pd.DataFrame, data_x: ad.AnnData, data_y: ad.AnnData,
                      color: str, label: str = None, xlim: float = None,
                      ylim: float = None) -> matplotlib.collections.PathCollection:
    """Scatter plot of gene percent counts with a single categorical color.

    Args:
        ax: Matplotlib axes to draw on.
        shared_data: Output of compare_genes containing 'percent_counts_x' and
            'percent_counts_y' columns.
        data_x: AnnData object for the x-axis dataset (used for axis label).
        data_y: AnnData object for the y-axis dataset (used for axis label).
        color: Matplotlib color string applied to all points.
        label: Legend label for this scatter series. Defaults to None.
        xlim: Upper x-axis limit. Defaults to None (auto).
        ylim: Upper y-axis limit. Defaults to None (auto).

    Returns:
        The PathCollection returned by ax.scatter.
    """
    x_percent = shared_data['percent_counts_x']
    y_percent = shared_data['percent_counts_y']

    plot = ax.scatter(x_percent,
                      y_percent,
                      s=50,
                      alpha=0.5,
                      color=color,
                      label=label)

    ax.set_xlabel(data_x.uns['title'] + ' Gene Percent Count')
    ax.set_ylabel(data_y.uns['title'] + ' Gene Percent Count')
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
        shared_data: DataFrame with 'percent_counts_x' and 'percent_counts_y' columns.
    """
    x_percent = shared_data['percent_counts_x']
    y_percent = shared_data['percent_counts_y']

    pearson_r = pearsonr(x_percent, y_percent).correlation
    spearman_r = spearmanr(x_percent, y_percent).correlation

    x_mean = x_percent.mean()
    y_mean = y_percent.mean()
    x_std = x_percent.std()
    y_std = y_percent.std()
    CCC_r = 2 * pearson_r * x_std * y_std / (x_std**2 + y_std**2 + (x_mean - y_mean)**2)

    ax.plot([0, 100], [0, 100], color='black', linestyle='--', linewidth=1.5)

    textstr = '\n'.join((
        r'$\mathrm{r}=%.2f$' % (pearson_r,),
        r'$\rho=%.2f$' % (spearman_r,),
        r'$\rho_C=%.2f$' % (CCC_r,),))
    ax.text(0.05, 0.85, textstr, transform=ax.transAxes, fontsize=12)


def _label_genes(ax: matplotlib.axes.Axes, df: pd.DataFrame,
                 n_cooks: int, n_log_ratio: int, min_pct: float) -> None:
    eps = 1e-10
    expressed = df[(df['percent_counts_x'] >= min_pct) | (df['percent_counts_y'] >= min_pct)]
    log_ratio = np.abs(np.log((expressed['percent_counts_x'] + eps) / (expressed['percent_counts_y'] + eps)))
    top_cooks_idx = df.nlargest(n_cooks, 'cooks_distance').index
    top_ratio_idx = log_ratio.nlargest(n_log_ratio).index
    for idx in top_cooks_idx.union(top_ratio_idx):
        row = df.loc[idx]
        ax.annotate(row['gene_name'], (row['percent_counts_x'], row['percent_counts_y']),
                    fontsize=7, ha='left', va='bottom', clip_on=True)


def comparison_plotter(compare_dfs: list[pd.DataFrame], comparisons: list[tuple],
                       norm: matplotlib.colors.Normalize, lim: float,
                       metric: str, metric_name: str,
                       n_cooks: int = 5, n_log_ratio: int = 5, min_pct: float = 0.1) -> None:
    """Four-panel scatter plot comparing gene percent counts colored by a continuous metric.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        norm: Matplotlib normalisation applied to the color mapping.
        lim: Upper axis limit applied to both x and y axes.
        metric: Column in each DataFrame used to color points.
        metric_name: Label shown on the colorbar.
    """
    fig, axs = plt.subplots(1, 4, figsize=(25, 5))

    for ax, df, pair in zip(axs, compare_dfs, comparisons):
        plot = scatter_genes(ax,
                             df,
                             pair[0],
                             pair[1],
                             metric,
                             xlim=lim,
                             ylim=lim,
                             norm=norm)
        show_correlation(ax, df)
        _label_genes(ax, df, n_cooks, n_log_ratio, min_pct)

    fig.colorbar(plot, label=metric_name)

    plt.tight_layout()
    plt.show()


def compare_by_density(compare_dfs: list[pd.DataFrame], 
                       comparisons: list[tuple], lim: float,
                       n_cooks: int = 5, n_log_ratio: int = 5, min_pct: float = 0.1) -> None:
    """Four-panel comparison scatter plot colored by point density.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        lim: Upper axis limit applied to both x and y axes.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log percent ratio to label.
        min_pct: Genes where both percent_counts are below this value are excluded
            from the log-ratio ranking.
    """
    c_values = []
    for df in compare_dfs:
        c_values.extend(df['point_density'].tolist())
    norm = LogNorm(min(c_values), max(c_values))

    comparison_plotter(compare_dfs, comparisons, norm, lim, 'point_density', 'Density',
                       n_cooks, n_log_ratio, min_pct)


def compare_by_cooks(compare_dfs: list[pd.DataFrame], comparisons: list[tuple], lim: float,
                     n_cooks: int = 5, n_log_ratio: int = 5, min_pct: float = 0.1) -> None:
    """Four-panel comparison scatter plot colored by Cook's distance.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        lim: Upper axis limit applied to both x and y axes.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log percent ratio to label.
        min_pct: Genes where both percent_counts are below this value are excluded
            from the log-ratio ranking.
    """
    c_values = []
    for df in compare_dfs:
        c_values.extend(df['cooks_distance'].tolist())
    norm = LogNorm(min(c_values), max(c_values))

    comparison_plotter(compare_dfs, comparisons, norm, lim, 'cooks_distance', "Cook's Distance",
                       n_cooks, n_log_ratio, min_pct)


def compare_by_length(compare_dfs: list[pd.DataFrame], comparisons: list[tuple], lim: float,
                      n_cooks: int = 5, n_log_ratio: int = 5, min_pct: float = 0.1) -> None:
    """Four-panel comparison scatter plot colored by gene length.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        lim: Upper axis limit applied to both x and y axes.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log percent ratio to label.
        min_pct: Genes where both percent_counts are below this value are excluded
            from the log-ratio ranking.
    """
    c_values = []
    for df in compare_dfs:
        df['gene_length'] = df['gene_length'] + 1
        c_values.extend(df['gene_length'].tolist())
    norm = LogNorm(1, max(c_values))

    comparison_plotter(compare_dfs, comparisons, norm, lim, 'gene_length', 'Gene Length',
                       n_cooks, n_log_ratio, min_pct)


def compare_by_gc(compare_dfs: list[pd.DataFrame], comparisons: list[tuple], lim: float,
                  n_cooks: int = 5, n_log_ratio: int = 5, min_pct: float = 0.1) -> None:
    """Four-panel comparison scatter plot colored by GC content.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        lim: Upper axis limit applied to both x and y axes.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log percent ratio to label.
        min_pct: Genes where both percent_counts are below this value are excluded
            from the log-ratio ranking.
    """
    c_values = []
    for df in compare_dfs:
        c_values.extend(df['gc_content'].tolist())
    norm = Normalize(0, max(c_values))

    comparison_plotter(compare_dfs, comparisons, norm, lim, 'gc_content', 'Percent GC Content',
                       n_cooks, n_log_ratio, min_pct)


def compare_by_type(compare_dfs: list[pd.DataFrame], comparisons: list[tuple], lim: float,
                    n_cooks: int = 5, n_log_ratio: int = 5, min_pct: float = 0.1) -> None:
    """Four-panel comparison scatter plot with categorical colors by gene biotype.

    Colors: protein-coding (yellow), mtRNA (blue), rRNA (red), lncRNA (green),
    unspecified non-coding (black).

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        lim: Upper axis limit applied to both x and y axes.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log percent ratio to label.
        min_pct: Genes where both percent_counts are below this value are excluded
            from the log-ratio ranking.
    """
    fig, axs = plt.subplots(1, 4, figsize=(25, 5))

    for ax, df, pair in zip(axs, compare_dfs, comparisons):
        cat_scatter_genes(ax, df[~(df['is_mito'] | df['is_ribo'] | df['is_lnc'])],
                          pair[0], pair[1], 'black',
                          label='unspecified non-coding', xlim=lim, ylim=lim)
        for col, c, label in zip(['is_pc', 'is_mito', 'is_ribo', 'is_lnc'],
                                  ['yellow', 'blue', 'red', 'green'],
                                  ['protein coding', 'mtRNA', 'rRNA', 'lncRNA']):
            cat_scatter_genes(ax, df[df[col]], pair[0], pair[1],
                              c, label=label, xlim=lim, ylim=lim)

        show_correlation(ax, df)
        _label_genes(ax, df, n_cooks, n_log_ratio, min_pct)

    axs[0].legend()

    plt.tight_layout()
    plt.show()


def compare_by_pathway(compare_dfs: list[pd.DataFrame], comparisons: list[tuple], lim: float,
                       n_cooks: int = 10, n_log_ratio: int = 10, min_pct: float = 0.1) -> None:
    """Four-panel comparison scatter plot highlighting ribosomal and OXPHOS genes.

    All genes are drawn in light grey; ribosomal (RP/MRP) genes are overlaid in
    red and OXPHOS genes (SDH, UQCR, COX, NDUF, ATP, and their mitochondria-encoded
    equivalents) in blue.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        lim: Upper axis limit applied to both x and y axes.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log percent ratio to label.
        min_pct: Genes where both percent_counts are below this value are excluded
            from the log-ratio ranking.
    """
    fig, axs = plt.subplots(1, 4, figsize=(25, 5))

    for ax, df, pair in zip(axs, compare_dfs, comparisons):
        cat_scatter_genes(ax, df, pair[0], pair[1], 'lightgrey', xlim=lim, ylim=lim)
        cat_scatter_genes(ax, df[df['is_ribo']], pair[0], pair[1],
                          'red', label='ribosomal (RP/MRP)', xlim=lim, ylim=lim)
        cat_scatter_genes(ax, df[df['is_oxphos']], pair[0], pair[1],
                          'blue', label='OXPHOS', xlim=lim, ylim=lim)

        show_correlation(ax, df)
        _label_genes(ax, df, n_cooks, n_log_ratio, min_pct)

    axs[0].legend()

    plt.tight_layout()
    plt.show()


def compare_by_enrichment(compare_dfs: list[pd.DataFrame], comparisons: list[tuple], lim: float,
                           enrichment_path: str | Path,
                           gene_sets: list[str] | None = None,
                           terms: list[str] | None = None,
                           n_cooks: int = 10, n_log_ratio: int = 10, min_pct: float = 0.1,
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
        lim: Upper axis limit applied to both x and y axes.
        enrichment_path: Path to an enrichment_results.csv from gseapy.enrichr.
        gene_sets: If provided, restrict to rows whose Gene_set is in this list.
        terms: If provided, restrict to rows whose Term is in this list.
        n_cooks: Number of top genes by Cook's distance to label.
        n_log_ratio: Number of top genes by absolute log percent ratio to label.
        min_pct: Genes where both percent_counts are below this value are excluded
            from the log-ratio ranking.
        label_enriched: If True, annotate every gene belonging to any selected term
            with its gene name, in addition to the top Cook's/log-ratio labels.
    """
    def _norm(name: str) -> str:
        for prefix in ("HUMAN_", "MOUSE_"):
            if name.startswith(prefix):
                return name[len(prefix):].upper()
        return name.upper()

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

    cmap = matplotlib.colormaps.get_cmap('tab10')
    colors = [cmap(i % 10) for i in range(len(selected_terms))]

    all_enriched: set[str] = set().union(*per_term_genes.values()) if per_term_genes else set()

    fig, axs = plt.subplots(1, 4, figsize=(25, 5))

    for ax, df, pair in zip(axs, compare_dfs, comparisons):
        norm_names = df['gene_name'].apply(_norm)
        cat_scatter_genes(ax, df, pair[0], pair[1], 'lightgrey', xlim=lim, ylim=lim)

        for term, color in zip(selected_terms, colors):
            mask = norm_names.isin(per_term_genes[term])
            cat_scatter_genes(ax, df[mask], pair[0], pair[1],
                              color, label=term, xlim=lim, ylim=lim)

        show_correlation(ax, df)
        _label_genes(ax, df, n_cooks, n_log_ratio, min_pct)

        if label_enriched:
            for idx in df[norm_names.isin(all_enriched)].index:
                row = df.loc[idx]
                ax.annotate(row['gene_name'],
                            (row['percent_counts_x'], row['percent_counts_y']),
                            fontsize=7, ha='left', va='bottom', clip_on=True)

    axs[0].legend()

    plt.tight_layout()
    plt.show()

def compare(datasets: list[ad.AnnData], 
            lim: float,
            comparison_axis: str = "percent_counts",
            n_cooks: int = 5, 
            n_log_ratio: int = 10, 
            min_pct: float = 0.1
) -> tuple[list[str], list[pd.DataFrame]]:
    """Generate all pairwise gene-count comparison plots across four sequencing methods.

    Filters each dataset to genes with percent_counts below lim, then computes
    compare_genes for every pair and renders density, Cook's, length, GC, and
    biotype comparison plots.

    Comparison pairs: polyT vs randO, 10X vs polyT, 10X vs randO, 10X vs Parse.

    Args:
        datasets: List of four AnnData objects in order [10X, polyT, randO, Parse].
        lim: Maximum percent_counts value; genes above this threshold are excluded.
        n_cooks: Number of top genes by Cook's distance to label on each plot.
        n_log_ratio: Number of top genes by absolute log percent ratio to label on each plot.
        min_pct: Genes where both percent_counts are below this value are excluded
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
        compare_dfs.append(processing.compare_genes(pair[0], pair[1], lim, comparison_axis))

    compare_by_density(compare_dfs, comparisons, lim, n_cooks, n_log_ratio, min_pct)
    compare_by_cooks(compare_dfs, comparisons, lim, n_cooks, n_log_ratio, min_pct)
    compare_by_length(compare_dfs, comparisons, lim, n_cooks, n_log_ratio, min_pct)
    compare_by_gc(compare_dfs, comparisons, lim, n_cooks, n_log_ratio, min_pct)
    compare_by_type(compare_dfs, comparisons, lim, n_cooks, n_log_ratio, min_pct)
    compare_by_pathway(compare_dfs, comparisons, lim, n_cooks, n_log_ratio, min_pct)

    return compare_names, compare_dfs


def generate_upset(datasets: list[ad.AnnData], gene_info: pd.DataFrame, cell_thresh: int = 10, n_top_genes: int = -1) -> pd.DataFrame:
    """Generate an UpSet plot comparing genes expressed across sequencing methods.

    Args:
        datasets: List of AnnData objects, each with var containing 'n_cells',
            'percent_counts', and 'gene_id', and uns containing 'name'.
        gene_info: DataFrame from query_ensembl with a 'gene_id' column.
        cell_thresh: Minimum number of cells a gene must be detected in to be
            considered expressed.
        n_top_genes: If positive, restricts each dataset to its top n genes by
            percent_counts before building the UpSet membership.

    Returns:
        DataFrame with one boolean column per dataset indicating gene membership,
        plus a 'gene_ids' column.
    """
    contents = pd.DataFrame()

    for data in datasets:
        var = data.var[data.var['n_cells'] > cell_thresh]
        if n_top_genes > 0:
            top_data = var.sort_values(by='percent_counts', ascending=False).head(1000)
            mask = (gene_info['gene_id'].isin(top_data['gene_id'])).tolist()
        else:
            mask = (gene_info['gene_id'].isin(var['gene_id'])).tolist()
        contents[str(data.uns['name'])] = mask

    upset = Upset.generate_plot(contents)
    contents['gene_ids'] = gene_info['gene_id']
    upset.show()
    return contents


def plot_geneset_metrics(gene_sets: list[pd.Series], set_names: list[str], gene_info: pd.DataFrame) -> None:
    """Violin plots of gene length and GC content for each gene set.

    Args:
        gene_sets: List of Series containing gene IDs, one per set to display.
        set_names: Display labels for each gene set column.
        gene_info: DataFrame from query_ensembl with 'gene_id', 'gene_length',
            and 'gc_content' columns.
    """
    fig, axs = plt.subplots(2, len(gene_sets), figsize=(12, 10), sharey='row')
    for ax, data, name in zip(axs[0, :], gene_sets, set_names):
        lengths = gene_info[gene_info['gene_id'].isin(data)]['gene_length'].tolist()
        ax.violinplot(np.log10(lengths), showextrema=False, showmedians=True)
        ax.set_title(name, fontsize=10)
        ax.set_xticks([])

    for ax, data, name in zip(axs[1, :], gene_sets, set_names):
        gcs = gene_info[gene_info['gene_id'].isin(data)]['gc_content'].tolist()
        ax.violinplot(gcs, showextrema=False, showmedians=True)
        ax.set_title(name, fontsize=10)
        ax.set_xticks([])

    axs[0, 0].set_ylabel('Gene Length (by order of magnitude)')
    axs[1, 0].set_ylabel('Percent GC Content')

    plt.tight_layout()
    plt.show()


def plot_genetype_counts(gene_sets: list[pd.Series], set_names: list[str], cols: list[str],
                         col_names: list[str], color: list[str], gene_info: pd.DataFrame) -> None:
    """Bar plots showing gene type counts for each gene set.

    Args:
        gene_sets: List of Series containing gene IDs, one per set to display.
        set_names: Display labels for each panel.
        cols: gene_info boolean columns to sum (e.g. 'is_pc', 'is_mito').
        col_names: Bar labels corresponding to each column in cols.
        color: Bar colors corresponding to each column in cols.
        gene_info: DataFrame from query_ensembl with 'gene_id' and the columns
            listed in cols.
    """
    fig, axs = plt.subplots(1, len(gene_sets), figsize=(25, 5))
    for ax, data, name in zip(axs, gene_sets, set_names):
        sums = []
        for col in cols:
            sums.append(gene_info[col][gene_info['gene_id'].isin(data)].sum())
        ax.bar(col_names, sums, color=color)
        ax.set_title(name, fontsize=10)

    axs[0].set_ylabel('Number of Genes In Intersection')

    plt.tight_layout()
    plt.show()
