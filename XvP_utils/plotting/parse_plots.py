import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

# Plot 3' count ratios by well in violin plot
def violinByWell(ax: matplotlib.axes.Axes, obs: pd.DataFrame) -> None:
    well_count_ratios = obs.groupby('well')['polyT_counts'].sum() / obs.groupby('well')['n_counts'].sum()
    well_count_ratios = well_count_ratios[well_count_ratios >= 0]  

    ax.violinplot(well_count_ratios.tolist(), showmeans=True, showmedians=False, showextrema=True)
    ax.set_xticks([])
    ax.set_ylabel("3' counts to total counts ratio by well") 
    string = "mean: " + str(round(np.mean(well_count_ratios),2)) + "\nmedian: " + str(round(np.median(well_count_ratios),2)) + "\nstd: " + str(round(np.std(well_count_ratios),2))
    ax.text(0.04, 0.87, string, transform=ax.transAxes, fontsize=12)

# Plot 3' count ratios by cell in violin plot
def violinByCell(ax: matplotlib.axes.Axes, obs: pd.DataFrame) -> None:
    cell_count_ratios = obs['polyT_counts'] / (obs['polyT_counts'] + obs['randO_counts'])
    cell_count_ratios = cell_count_ratios.fillna(0)

    ax.violinplot(cell_count_ratios.tolist(), showmeans=True, showmedians=False, showextrema=True)
    ax.set_xticks([])
    ax.set_ylabel("3' counts to total counts ratio by cell") 
    string = "mean: " + str(round(np.mean(cell_count_ratios),2)) + "\nmedian: " + str(round(np.median(cell_count_ratios),2)) + "\nstd: " + str(round(np.std(cell_count_ratios),2))
    ax.text(0.04, 0.87, string, transform=ax.transAxes, fontsize=12)

# Calculates the Concordance Correlation Coefficient (CCC) between two array-like objects
def CCC(x, y):
    pearson_r = pearsonr(x,y).correlation
    return 2 * pearson_r * x.std() * y.std() / (x.std()**2 + y.std()**2 + (x.mean() - y.mean())**2)

# Scatter plot comparing polyT to randO by well
def plotByWell(ax: matplotlib.axes.Axes, obs: pd.DataFrame, xlim: float = 1.15, ylim: float = 1.15, corr: bool = True) -> None:
    polyT_well_counts = obs.groupby('well')['polyT_counts'].sum()
    randO_well_counts = obs.groupby('well')['randO_counts'].sum()
    well_counts_df = pd.concat([polyT_well_counts, randO_well_counts], axis=1).fillna(0)
    well_counts_df = well_counts_df.set_axis(['polyT_counts', 'randO_counts'], axis=1)

    ax.scatter(well_counts_df['randO_counts'], well_counts_df['polyT_counts'])
    ax.set_xlabel('random oligo total counts per well')
    ax.set_ylabel('3\' total counts per well')
    max_lim = max(well_counts_df['randO_counts'].max(), well_counts_df['polyT_counts'].max())
    ax.plot([0, 1.05*max_lim], [0, 1.05*max_lim], color='black', linestyle='--', linewidth=1) # x=y line
    if corr:
        CCC_r = CCC(well_counts_df['randO_counts'], well_counts_df['polyT_counts'])
        ax.text(0.05, 0.9, r'$\rho_C=%.2f$' % (CCC_r, ), transform=ax.transAxes, fontsize=12) 

# Scatter plot comparing polyT to randO counts by cell
def plotByCell(ax: matplotlib.axes.Axes, obs: pd.DataFrame, xlim: float = 1.15, ylim: float = 1.15, corr: bool = True):
    ax.scatter(obs['randO_counts'], obs['polyT_counts'])
    ax.set_xlabel('random oligo total counts per cell')
    ax.set_ylabel('3\' total counts per cell')
    max_lim = max(obs['randO_counts'].max(), obs['polyT_counts'].max())
    ax.set_xlim(0, xlim*max_lim)
    ax.set_ylim(0, ylim*max_lim)
    ax.plot([0, 1.05*max_lim], [0, 1.05*max_lim], color='black', linestyle='--', linewidth=1) # x=y line
    if corr:
        CCC_r = CCC(obs['randO_counts'], obs['polyT_counts'])
        ax.text(0.05, 0.9, r'$\rho_C=%.2f$' % (CCC_r, ), transform=ax.transAxes, fontsize=12) 