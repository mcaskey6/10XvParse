import matplotlib.pyplot as plt
import matplotlib
import matplotlib.axes
import pandas as pd
import numpy as np
import anndata as ad
from scipy.stats import pearsonr


def plot_counts_by_well(ax: matplotlib.axes.Axes, data: ad.AnnData) -> None:
    """Bar plot of total polyT and randO counts grouped by well.

    Args:
        ax: Matplotlib axes to draw on.
        data: AnnData object whose obs contains 'well', 'polyT_counts', and 'randO_counts'.
    """
    well_counts = data.obs.groupby('well')[['polyT_counts', 'randO_counts']].sum()

    well_counts.plot(kind='bar', ax=ax)
    ax.set_title('Total polyT and randO Counts by Well')
    ax.set_xlabel('Well')
    ax.set_ylabel('Total Counts')
    ax.legend(title='Barcode Type')


def plot_num_cells_by_well(ax: matplotlib.axes.Axes, data: ad.AnnData) -> None:
    """Bar plot of unique cell counts per well.

    Args:
        ax: Matplotlib axes to draw on.
        data: AnnData object whose obs contains a 'well' column.
    """
    well_cell_counts = data.obs.groupby('well').size()

    well_cell_counts.plot(kind='bar', ax=ax)
    ax.set_title('Number of Cells by Well')
    ax.set_xlabel('Well')
    ax.set_ylabel('Number of Cells')


def violin_by_well(ax: matplotlib.axes.Axes, obs: pd.DataFrame) -> None:
    """Violin plot of polyT-to-total count ratios aggregated by well.

    Overlays mean, median, and std as a text annotation.

    Args:
        ax: Matplotlib axes to draw on.
        obs: Cell observation DataFrame with 'well', 'polyT_counts', and 'n_counts' columns.
    """
    well_count_ratios = obs.groupby('well')['polyT_counts'].sum() / obs.groupby('well')['n_counts'].sum()
    well_count_ratios = well_count_ratios[well_count_ratios >= 0]

    ax.violinplot(well_count_ratios.tolist(), showmeans=True, showmedians=False, showextrema=True)
    ax.set_xticks([])
    ax.set_ylabel("3' counts to total counts ratio by well")
    string = "mean: " + str(round(np.mean(well_count_ratios), 2)) + "\nmedian: " + str(round(np.median(well_count_ratios), 2)) + "\nstd: " + str(round(np.std(well_count_ratios), 2))
    ax.text(0.04, 0.87, string, transform=ax.transAxes, fontsize=12)


def violin_by_cell(ax: matplotlib.axes.Axes, obs: pd.DataFrame) -> None:
    """Violin plot of polyT-to-total count ratios per cell.

    Overlays mean, median, and std as a text annotation.

    Args:
        ax: Matplotlib axes to draw on.
        obs: Cell observation DataFrame with 'polyT_counts' and 'randO_counts' columns.
    """
    cell_count_ratios = obs['polyT_counts'] / (obs['polyT_counts'] + obs['randO_counts'])
    cell_count_ratios = cell_count_ratios.fillna(0)

    ax.violinplot(cell_count_ratios.tolist(), showmeans=True, showmedians=False, showextrema=True)
    ax.set_xticks([])
    ax.set_ylabel("3' counts to total counts ratio by cell")
    string = "mean: " + str(round(np.mean(cell_count_ratios), 2)) + "\nmedian: " + str(round(np.median(cell_count_ratios), 2)) + "\nstd: " + str(round(np.std(cell_count_ratios), 2))
    ax.text(0.04, 0.87, string, transform=ax.transAxes, fontsize=12)


def ccc(x: pd.Series, y: pd.Series) -> float:
    """Compute the Concordance Correlation Coefficient (CCC) between two series.

    Args:
        x: First data series.
        y: Second data series.

    Returns:
        CCC value in the range [-1, 1].
    """
    pearson_r = pearsonr(x, y).correlation
    return 2 * pearson_r * x.std() * y.std() / (x.std()**2 + y.std()**2 + (x.mean() - y.mean())**2)


def plot_by_well(ax: matplotlib.axes.Axes, obs: pd.DataFrame, xlim: float = 1.15, ylim: float = 1.15, corr: bool = True) -> None:
    """Scatter plot comparing total polyT vs randO counts per well.

    Args:
        ax: Matplotlib axes to draw on.
        obs: Cell observation DataFrame with 'well', 'polyT_counts', and 'randO_counts' columns.
        xlim: Multiplier applied to the x-axis upper limit relative to the data max.
        ylim: Multiplier applied to the y-axis upper limit relative to the data max.
        corr: If True, annotates the plot with the CCC correlation coefficient.
    """
    polyT_well_counts = obs.groupby('well')['polyT_counts'].sum()
    randO_well_counts = obs.groupby('well')['randO_counts'].sum()
    well_counts_df = pd.concat([polyT_well_counts, randO_well_counts], axis=1).fillna(0)
    well_counts_df = well_counts_df.set_axis(['polyT_counts', 'randO_counts'], axis=1)

    ax.scatter(well_counts_df['randO_counts'], well_counts_df['polyT_counts'])
    ax.set_xlabel('random oligo total counts per well')
    ax.set_ylabel('3\' total counts per well')
    max_lim = max(well_counts_df['randO_counts'].max(), well_counts_df['polyT_counts'].max())
    ax.set_xlim(0, xlim * max_lim)
    ax.set_ylim(0, ylim * max_lim)
    ax.plot([0, 1.05 * max_lim], [0, 1.05 * max_lim], color='black', linestyle='--', linewidth=1)
    if corr:
        CCC_r = ccc(well_counts_df['randO_counts'], well_counts_df['polyT_counts'])
        ax.text(0.05, 0.9, r'$\rho_C=%.2f$' % (CCC_r,), transform=ax.transAxes, fontsize=12)


def plot_by_cell(ax: matplotlib.axes.Axes, obs: pd.DataFrame, xlim: float = 1.15, ylim: float = 1.15, corr: bool = True) -> None:
    """Scatter plot comparing polyT vs randO counts per cell.

    Args:
        ax: Matplotlib axes to draw on.
        obs: Cell observation DataFrame with 'polyT_counts' and 'randO_counts' columns.
        xlim: Multiplier applied to the x-axis upper limit relative to the data max.
        ylim: Multiplier applied to the y-axis upper limit relative to the data max.
        corr: If True, annotates the plot with the CCC correlation coefficient.
    """
    ax.scatter(obs['randO_counts'], obs['polyT_counts'])
    ax.set_xlabel('random oligo total counts per cell')
    ax.set_ylabel('3\' total counts per cell')
    max_lim = max(obs['randO_counts'].max(), obs['polyT_counts'].max())
    ax.set_xlim(0, xlim * max_lim)
    ax.set_ylim(0, ylim * max_lim)
    ax.plot([0, 1.05 * max_lim], [0, 1.05 * max_lim], color='black', linestyle='--', linewidth=1)
    if corr:
        CCC_r = ccc(obs['randO_counts'], obs['polyT_counts'])
        ax.text(0.05, 0.9, r'$\rho_C=%.2f$' % (CCC_r,), transform=ax.transAxes, fontsize=12)
