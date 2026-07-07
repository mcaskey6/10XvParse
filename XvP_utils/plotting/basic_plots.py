import numpy as np
import matplotlib
import matplotlib.axes
import anndata as ad
import scanpy as sc
from . import processing
import matplotlib.pyplot as plt
from typing import Tuple


def scatter_reads(datasets: list[ad.AnnData], figsize: Tuple[float,float]=(20,5)) -> None:
    """Grid of scatter plots of genes detected vs UMI counts on log-log axes.

    Args:
        datasets: A list of anndata objects with 'percent_nascent' in obs
        figsize: Size of the figure
    """

    fig, axs = plt.subplots(1, len(datasets), figsize=figsize, sharey=True, sharex=True)

    for ax, data in zip(axs, datasets):
        x = np.asarray(data.X.sum(axis=1))[:, 0]
        y = np.asarray(np.sum(data.X > 0, axis=1))[:, 0]
        labels = data.obs["percent_nascent"].to_numpy()

        plot = ax.scatter(x, y, c=labels, alpha=0.25)
        ax.set_xlabel("UMI Counts")
        ax.set_xscale('log')
        ax.set_yscale('log', nonpositive='clip')
        ax.set_title(data.uns['title'])

    axs[0].set_ylabel("Genes Detected")
    fig.colorbar(plot, label="Percent Nascent Counts")
    plt.tight_layout()
    plt.show()

def knee_plot(raw_datasets: list[ad.AnnData], cutoffs: list[int], transform: bool = False, figsize: Tuple[float,float]=(27,5)) -> list[ad.AnnData]:
    """Grid of log-log knee plots with a vertical UMI threshold line.

    Prints the number of cells passing the threshold and returns a list of the filtered datasets.

    Args:
        raw_datasets: List of unfiltered AnnData objects with a 'title' key in uns.
        cutoffs: List of minimum UMI count thresholds for retaining a cell for each dataset
        transform: If True, applies CPM normalization and log1p transformation to the filtered data before returning.

    Returns:
        List of AnnData object filtered to cells at or above the specified cutoff and transformed if transform=True.
    """

    fig, axs = plt.subplots(1, len(raw_datasets), figsize=figsize, sharey=True, sharex=True, squeeze=False)
    axs = axs.flatten()

    datasets = []
    for ax, raw_data, cutoff in zip(axs, raw_datasets, cutoffs):
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
        data = processing.refilter(raw_data, knee[num_cells], transform=transform)
        datasets.append(data)

    axs[0].set_ylabel("Barcode Rank (in order of UMI counts)")

    plt.tight_layout()
    plt.show()

    return datasets


def mito_scatter(datasets: list[ad.AnnData], figsize: Tuple[float,float]=(20,5)) -> None:
    """Grid of scatter plots of UMI counts vs percent mitochondrial gene expression.

    Args:
        datasets: A list of anndata objects with 'percent_nascent' in obs
        figsize: Size of the figure
    """
    fig, axs = plt.subplots(1, len(datasets), figsize=figsize, sharey=True,sharex=True)
    for ax, data in zip(axs, datasets):
        sc.pl.scatter(data, x='n_counts', y='percent_mito', ax=ax, show=False)
        ax.set_title(data.uns['title'])
        ax.set_xlabel("UMI Counts")
        ax.set_ylabel("")

    axs[0].set_ylabel("Percent Mitochondrial Counts")

    plt.tight_layout()
    plt.show()
