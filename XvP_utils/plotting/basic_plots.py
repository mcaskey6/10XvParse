import numpy as np
import matplotlib
import matplotlib.axes
import anndata as ad
import scanpy as sc
from . import processing


def scatter_reads(ax: matplotlib.axes.Axes, data: ad.AnnData) -> None:
    """Scatter plot of genes detected vs UMI counts on log-log axes.

    Args:
        ax: Matplotlib axes to draw on.
        data: AnnData object with a 'title' key in uns.
    """
    x = np.asarray(data.X.sum(axis=1))[:, 0]
    y = np.asarray(np.sum(data.X > 0, axis=1))[:, 0]

    ax.scatter(x, y, color="green", alpha=0.25)
    ax.set_xlabel("UMI Counts")
    ax.set_xscale('log')
    ax.set_yscale('log', nonpositive='clip')
    ax.set_title(data.uns['title'] + " Reads")


def knee_plot(ax: matplotlib.axes.Axes, raw_data: ad.AnnData, cutoff: int = 100) -> ad.AnnData:
    """Log-log knee plot with a vertical UMI threshold line.

    Prints the number of cells passing the threshold and returns a filtered dataset.

    Args:
        ax: Matplotlib axes to draw on.
        raw_data: Unfiltered AnnData object with a 'title' key in uns.
        cutoff: Minimum UMI count threshold for retaining a cell.

    Returns:
        AnnData filtered to cells at or above the cutoff.
    """
    knee = np.sort((np.array(raw_data.X.sum(axis=1))).flatten())[::-1]
    cell_set = np.arange(len(knee))
    num_cells = cell_set[knee > cutoff][::-1][0]

    ax.loglog(knee, cell_set, linewidth=5, color="g")
    ax.axvline(x=cutoff, linewidth=3, color="k")
    ax.axhline(y=num_cells, linewidth=3, color="k")
    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("")
    ax.set_title(raw_data.uns['title'] + " Knee Plot")

    print(f"{num_cells:,.0f} cells passed the {cutoff} UMI threshold for {raw_data.uns['title']}")
    data = processing.refilter(raw_data, knee[num_cells])

    return data


def mito_scatter(ax: matplotlib.axes.Axes, data: ad.AnnData) -> None:
    """Scatter plot of UMI counts vs percent mitochondrial gene expression.

    Args:
        ax: Matplotlib axes to draw on.
        data: AnnData object with obs containing 'n_counts' and 'percent_mito',
            and a 'title' key in uns.
    """
    sc.pl.scatter(data, x='n_counts', y='percent_mito', ax=ax, show=False)
    ax.set_title(data.uns['title'])
    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("")
