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
    filtered_counts = np.array([data.X.sum() for data in datasets])
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


def plot_cell_metrics(datasets: list[ad.AnnData], groups: list[str], group_names: list[str], figsize: tuple[int, int] = (12, 18)) -> None:
    """Grid of violin plots for multiple cell metrics across datasets.

    Rows correspond to metrics, columns to datasets.

    Args:
        datasets: List of AnnData objects to plot.
        groups: obs column names for each row of the grid.
        group_names: Display labels for the y-axis of each row.
        figsize: Figure size as (width, height) in inches.
    """
    fig, ax = plt.subplots(len(groups), len(datasets), figsize=(figsize), sharey='row')
    for i, data in enumerate(datasets):
        _violin_plots(ax[:, i], data, groups)

    for i, group in enumerate(group_names):
        ax[i, 0].set_ylabel(group)

    for i, data in enumerate(datasets):
        ax[0, i].set_title(data.uns['title'])
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


def marker_genes(ax: matplotlib.axes.Axes, data: ad.AnnData, markers: list[str]) -> None:
    """Violin plot of marker gene expression as a percent of total counts per cell.

    Args:
        ax: Matplotlib axes to draw on.
        data: AnnData object with obs containing 'n_counts' and a 'title' key in uns.
        markers: List of gene names to plot (must be present in data.var_names).
    """
    gene_dist = []
    for gene in markers:
        gene_dist.append(np.nan_to_num(data[:, gene].X.toarray().transpose()[0] / np.array(data.obs['n_counts'].tolist()) * 100))

    ax.violinplot(gene_dist, showmeans=True)
    ax.set_xticks(np.arange(1, len(markers) + 1), markers)
    ax.set_ylabel("")
    ax.set_title(data.uns['title'] + " Marker Genes")


def collapsed_marker_genes(ax: matplotlib.axes.Axes, data: ad.AnnData, marker_prefixes: list[str]) -> None:
    """Violin plot of collapsed marker gene family expression per cell.

    Sums counts across all genes sharing a common prefix, useful when individual
    gene distinctions are not important (e.g. ribosomal protein families).

    Args:
        ax: Matplotlib axes to draw on.
        data: AnnData object with obs containing 'n_counts' and a 'title' key in uns.
        marker_prefixes: Gene name prefixes to collapse. Each prefix becomes one
            violin in the plot.
    """
    genes = data.var

    gene_dist = []
    for prefix in marker_prefixes:
        mask = genes[genes.index.str.startswith(prefix)].index
        gene_dist.append(np.nan_to_num(data[:, mask].X.toarray().sum(axis=1).transpose() / np.array(data.obs['n_counts'].tolist()) * 100))

    ax.violinplot(gene_dist, showmeans=True)
    ax.set_xticks(np.arange(1, len(marker_prefixes) + 1), marker_prefixes)
    ax.set_ylabel("")
    ax.set_title(data.uns['title'] + " Thymus Marker Genes")


def top_gene_cell_expression(ax: matplotlib.axes.Axes, data: ad.AnnData) -> None:
    """Violin plot of the top 10 genes by number of cells expressing them.

    Args:
        ax: Matplotlib axes to draw on.
        data: AnnData object with var containing 'n_cells' and obs containing
            'n_counts'. Must have a 'title' key in uns.
    """
    genes = data.var
    top_genes = genes.sort_values(by='n_cells', ascending=False).head(10).index

    gene_dist = []
    for gene in top_genes:
        gene_dist.append(np.nan_to_num(data[:, gene].X.toarray().transpose()[0] / np.array(data.obs['n_counts'].tolist()) * 100))

    ax.violinplot(gene_dist, showmeans=True)
    ax.set_xticks(np.arange(1, len(top_genes) + 1), top_genes)
    ax.set_ylabel("")
    ax.set_title(data.uns['title'] + " Top 10 Genes by Cell Expression")


def top_gene_counts(ax: matplotlib.axes.Axes, data: ad.AnnData) -> None:
    """Violin plot of the top 10 genes by total percent counts.

    Args:
        ax: Matplotlib axes to draw on.
        data: AnnData object with var containing 'percent_counts' and obs containing
            'n_counts'. Must have a 'title' key in uns.
    """
    genes = data.var
    top_genes = genes.sort_values(by='percent_counts', ascending=False).head(10).index

    gene_dist = []
    for gene in top_genes:
        gene_dist.append(np.nan_to_num(data[:, gene].X.toarray().transpose()[0] / np.array(data.obs['n_counts'].tolist()) * 100))

    ax.violinplot(gene_dist, showmeans=True)
    ax.set_xticks(np.arange(1, len(top_genes) + 1), top_genes)
    ax.set_ylabel("")
    ax.set_title(data.uns['title'] + " Top 10 Genes by Total Counts")


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
                      c=color,
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


def comparison_plotter(compare_dfs: list[pd.DataFrame], comparisons: list[tuple],
                       norm: matplotlib.colors.Normalize, lim: float,
                       metric: str, metric_name: str) -> None:
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

    fig.colorbar(plot, label=metric_name)

    plt.tight_layout()
    plt.show()


def compare_by_density(compare_dfs: list[pd.DataFrame], comparisons: list[tuple], lim: float) -> None:
    """Four-panel comparison scatter plot colored by point density.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        lim: Upper axis limit applied to both x and y axes.
    """
    c_values = []
    for df in compare_dfs:
        c_values.extend(df['point_density'].tolist())
    norm = LogNorm(min(c_values), max(c_values))

    comparison_plotter(compare_dfs, comparisons, norm, lim, 'point_density', 'Density')


def compare_by_cooks(compare_dfs: list[pd.DataFrame], comparisons: list[tuple], lim: float) -> None:
    """Four-panel comparison scatter plot colored by Cook's distance.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        lim: Upper axis limit applied to both x and y axes.
    """
    c_values = []
    for df in compare_dfs:
        c_values.extend(df['cooks_distance'].tolist())
    norm = LogNorm(min(c_values), max(c_values))

    comparison_plotter(compare_dfs, comparisons, norm, lim, 'cooks_distance', "Cook's Distance")


def compare_by_length(compare_dfs: list[pd.DataFrame], comparisons: list[tuple], lim: float) -> None:
    """Four-panel comparison scatter plot colored by gene length.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        lim: Upper axis limit applied to both x and y axes.
    """
    c_values = []
    for df in compare_dfs:
        df['gene_length'] = df['gene_length'] + 1
        c_values.extend(df['gene_length'].tolist())
    norm = LogNorm(1, max(c_values))

    comparison_plotter(compare_dfs, comparisons, norm, lim, 'gene_length', 'Gene Length')


def compare_by_gc(compare_dfs: list[pd.DataFrame], comparisons: list[tuple], lim: float) -> None:
    """Four-panel comparison scatter plot colored by GC content.

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        lim: Upper axis limit applied to both x and y axes.
    """
    c_values = []
    for df in compare_dfs:
        c_values.extend(df['gc_content'].tolist())
    norm = Normalize(0, max(c_values))

    comparison_plotter(compare_dfs, comparisons, norm, lim, 'gc_content', 'Percent GC Content')


def compare_by_type(compare_dfs: list[pd.DataFrame], comparisons: list[tuple], lim: float) -> None:
    """Four-panel comparison scatter plot with categorical colors by gene biotype.

    Colors: protein-coding (yellow), mtRNA (blue), rRNA (red), lncRNA (green),
    unspecified non-coding (black).

    Args:
        compare_dfs: List of DataFrames from compare_genes, one per comparison pair.
        comparisons: List of (data_x, data_y) AnnData tuples matching compare_dfs.
        lim: Upper axis limit applied to both x and y axes.
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

    axs[0].legend()

    plt.tight_layout()
    plt.show()


def compare(datasets: list[ad.AnnData], lim: float) -> tuple[list[str], list[pd.DataFrame]]:
    """Generate all pairwise gene-count comparison plots across four sequencing methods.

    Filters each dataset to genes with percent_counts below lim, then computes
    compare_genes for every pair and renders density, Cook's, length, GC, and
    biotype comparison plots.

    Comparison pairs: polyT vs randO, 10X vs polyT, 10X vs randO, 10X vs Parse.

    Args:
        datasets: List of four AnnData objects in order [10X, polyT, randO, Parse].
        lim: Maximum percent_counts value; genes above this threshold are excluded.

    Returns:
        Tuple of (compare_names, compare_dfs) where compare_names is a list of
        pair label strings and compare_dfs is the corresponding list of DataFrames
        from compare_genes.
    """
    data_10x = datasets[0][:, datasets[0].var['percent_counts'] < lim]
    data_polyT = datasets[1][:, datasets[1].var['percent_counts'] < lim]
    data_randO = datasets[2][:, datasets[2].var['percent_counts'] < lim]
    data_parse = datasets[3][:, datasets[3].var['percent_counts'] < lim]

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
        compare_dfs.append(processing.compare_genes(pair[0], pair[1]))

    compare_by_density(compare_dfs, comparisons, lim)
    compare_by_cooks(compare_dfs, comparisons, lim)
    compare_by_length(compare_dfs, comparisons, lim)
    compare_by_gc(compare_dfs, comparisons, lim)
    compare_by_type(compare_dfs, comparisons, lim)

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
