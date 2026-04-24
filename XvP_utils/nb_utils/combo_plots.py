# Import packages
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from matplotlib.colors import LogNorm, Normalize
from upsetty import Upset
from . import processing

# Generates violin plots of selected cell metadata
def violin_plots(ax_col, data, groups):
    for i, group in enumerate(groups):
        sc.pl.violin(data, group, ax=ax_col[i], stripplot=False, show=False)
        ax_col[i].set_xticks([])

def plot_upsample(datasets: list[ad.AnnData], xlim:int = None, ylim:int = None) -> None:
    plt.figure(figsize=(16, 4))
    for data in datasets:
        pred = data.uns['pred_upsample']
        plt.plot(pred[:, 0], pred[:, 1], label=data.uns['title'])
        plt.fill_between(pred[:, 0], pred[:,2], pred[:,3], alpha=0.5)

    plt.legend(loc='upper left')
    if xlim:
        plt.xlim(0, xlim)
    if ylim:
        plt.ylim(0, ylim)
    plt.xlabel('Total Number of Counts')
    plt.ylabel('Projected Number of Unique Transcripts')
    plt.show()

def plot_filtering_metrics(datasets: list[ad.AnnData]):
    fig, axs = plt.subplots(1,2, figsize = (12,5), sharey=True)

    processed = np.array([data.uns['n_processed'] for data in datasets])
    aligned = np.array([data.uns['n_aligned'] for data in datasets])
    counts = np.array([data.uns['n_raw_counts'] for data in datasets])
    filtered_counts = np.array([data.X.sum() for data in datasets])
    labels = [data.uns['title'] for data in datasets]

    unaligned = processed - aligned
    rejected_counts = counts - filtered_counts

    axs[0].bar(labels, aligned, label = "Aligned", color = "orange")
    axs[0].bar(labels, unaligned, bottom=aligned, label = "Unaligned", color = "red")
    axs[0].set_ylabel('Number of Reads')
    axs[0].legend()

    axs[1].bar(labels, filtered_counts, label = "Filtered", color = 'green')
    axs[1].bar(labels, rejected_counts, bottom = filtered_counts, label = "Rejected", color = "blue")
    axs[1].set_ylabel('Number of Counts')
    axs[1].legend()

def plot_cell_metrics(datasets: list[ad.AnnData], groups: list[str], group_names: list[str], figsize: set[int] = (12,18)):
    fig, ax = plt.subplots(len(groups), len(datasets), figsize=(figsize), sharey='row')
    for i, data in enumerate(datasets):
        violin_plots(ax[:, i], data, groups)

    for i, group in enumerate(group_names):
        ax[i,0].set_ylabel(group)

    for i, data in enumerate(datasets):
        ax[0,i].set_title(data.uns['title'])
    plt.tight_layout()
    plt.show()

def plot_gene_metrics(datasets: ad.AnnData, sample_size:int = 1000000):  
    fig, ax = plt.subplots(2, len(datasets), figsize=(12,10), sharey='row')

    for i, data in enumerate(datasets):
        gene_lengths = data.var['gene_length'].values
        gc_content =  data.var['gc_content'].values
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

# Displays violin plots of the distributions of 
# of the specified marker genes by cell
def marker_genes(ax, data, markers):
    gene_dist = []
    for gene in markers:
        gene_dist.append(np.nan_to_num(data[:,gene].X.toarray().transpose()[0] / np.array(data.obs['n_counts'].tolist()) * 100))

    ax.violinplot(gene_dist, showmeans=True)
    ax.set_xticks(np.arange(1, len(markers) + 1), markers)
    ax.set_ylabel("")

    ax.set_title(data.uns['title'] + " Thymus Marker Genes")

# Displays violin plots of the distributions of 
# of the specified marker genes by cell. Specific to
# similar sets of genes for which it is not important
# to differentiate between so counts can be collapsed
def collapsed_marker_genes(ax, data, marker_prefixes):
    genes = data.var

    gene_dist = []
    for prefix in marker_prefixes:
        mask = genes[genes.index.str.startswith(prefix)].index
        gene_dist.append(np.nan_to_num(data[:,mask].X.toarray().sum(axis=1).transpose() / np.array(data.obs['n_counts'].tolist()) * 100))

    ax.violinplot(gene_dist, showmeans=True)
    ax.set_xticks(np.arange(1, len(marker_prefixes) + 1), marker_prefixes)
    ax.set_ylabel("")

    ax.set_title(data.uns['title'] + " Thymus Marker Genes")

# Displays violin plots of the distributions of the top 10 expressed gene/transcripts
# (by number of cells expressed in) by cell
def top_gene_cell_expression(ax, data):
    genes = data.var
    top_genes = genes.sort_values(by='n_cells',ascending=False).head(10).index

    gene_dist = []
    for gene in top_genes:
        gene_dist.append(np.nan_to_num(data[:,gene].X.toarray().transpose()[0] / np.array(data.obs['n_counts'].tolist()) * 100))

    ax.violinplot(gene_dist, showmeans=True)
    ax.set_xticks(np.arange(1, len(top_genes) + 1), top_genes)
    ax.set_ylabel("")

    ax.set_title(data.uns['title'] + " Top 10 Genes by Cell Expression")

# Displays violin plots of the distributions of the top 10 (by total counts) expressed 
# gene/transcripts by cell
def top_gene_counts(ax, data):
    genes = data.var
    top_genes = genes.sort_values(by='percent_counts',ascending=False).head(10).index

    gene_dist = []
    for gene in top_genes:
        gene_dist.append(np.nan_to_num(data[:,gene].X.toarray().transpose()[0] / np.array(data.obs['n_counts'].tolist()) * 100))

    ax.violinplot(gene_dist, showmeans=True)
    ax.set_xticks(np.arange(1, len(top_genes) + 1), top_genes)
    ax.set_ylabel("")
    ax.set_title(data.uns['title'] + " Top 10 Genes by Total Counts")

# From the output of compare_genes, generates a scatter plot of the gene counts percentages 
# with each method on an axis
def scatter_genes(ax, shared_data, data_x, data_y, c_column,  xlim = None, ylim = None, norm = None):
    x_percent = shared_data['percent_counts_x']
    y_percent = shared_data['percent_counts_y']

    plot = ax.scatter(x_percent, 
                y_percent, 
                norm = norm, 
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

# From the output of compare_genes, generates a scatter plot of the gene counts percentages 
# with each method on an axis
def cat_scatter_genes(ax, shared_data, data_x, data_y, color, label=None, xlim = None, ylim = None):
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

# Displays correlation of percent counts from the output of compare_genes on the plot 
# given by ax
def show_correlation(ax, shared_data):
    x_percent = shared_data['percent_counts_x']
    y_percent = shared_data['percent_counts_y']

    # Calculate correlation
    pearson_r = pearsonr(x_percent, y_percent).correlation
    spearman_r = spearmanr(x_percent, y_percent).correlation
    
    x_mean = x_percent.mean()
    y_mean = y_percent.mean()
    x_std = x_percent.std()
    y_std = y_percent.std()
    CCC_r = 2 * pearson_r * x_std * y_std / (x_std**2 + y_std**2 + (x_mean - y_mean)**2)    

    ax.plot([0, 100], [0, 100], color='black', linestyle='--', linewidth=1.5)  # y=x line

    # Add plot annotations
    textstr = '\n'.join((
    r'$\mathrm{r}=%.2f$' % (pearson_r, ),
    r'$\rho=%.2f$' % (spearman_r, ),
    r'$\rho_C=%.2f$' % (CCC_r, ),))
    ax.text(0.05, 0.85, textstr, transform = ax.transAxes, fontsize = 12)

    return

# Plot a scatter plot of comparing gene percent counts between methods given 
# some metric with which to color the plot
def comparisonPlotter(compare_dfs, comparisons, norm, lim, metric, metric_name):
    fig, axs = plt.subplots(1, 4, figsize = (25, 5))

    for ax, df , pair in zip(axs, compare_dfs,comparisons):
        plot = scatter_genes(ax, 
                            df, 
                            pair[0],
                            pair[1], 
                            metric, 
                            xlim=lim, 
                            ylim=lim, 
                            norm=norm)
        show_correlation(ax,df)

    fig.colorbar(plot, label = metric_name)

    plt.tight_layout()
    plt.show()

# Compare gene percent counts between methods and 
# color by density
def compareByDensity(compare_dfs, comparisons, lim):
    c_values = []
    for df in compare_dfs:
        c_values.extend(df['point_density'].tolist())
    norm = LogNorm(min(c_values), max(c_values))

    comparisonPlotter(compare_dfs,
                      comparisons,
                      norm,
                      lim,
                      'point_density',
                      'Density')

# Compare gene percent counts between methods and 
# color by Cook's Distance
def compareByCooks(compare_dfs, comparisons, lim):
    c_values = []
    for df in compare_dfs:
        c_values.extend(df['cooks_distance'].tolist())
    norm = LogNorm(min(c_values), max(c_values))

    comparisonPlotter(compare_dfs,
                      comparisons,
                      norm,
                      lim,
                      'cooks_distance',
                      "Cook's Distance")
    
# Compare gene percent counts between methods and 
# color by gene length
def compareByLength(compare_dfs, comparisons, lim):
    c_values = []
    for df in compare_dfs:
        df['gene_length'] = df['gene_length']+1
        c_values.extend(df['gene_length'].tolist())
    norm = LogNorm(1, max(c_values))

    comparisonPlotter(compare_dfs,
                      comparisons,
                      norm,
                      lim,
                      'gene_length',
                      'Gene Length')

# Compare gene percent counts between methods and 
# color by GC content
def compareByGC(compare_dfs, comparisons, lim):
    c_values = []
    for df in compare_dfs:
        c_values.extend(df['gc_content'].tolist())
    norm = Normalize(0, max(c_values))

    comparisonPlotter(compare_dfs,
                      comparisons,
                      norm,
                      lim,
                      'gc_content',
                      'Percent GC Content')

# Compare gene percent counts between methods and 
# color by biotype
def compareByType(compare_dfs, comparisons, lim):
    fig, axs = plt.subplots(1, 4, figsize = (25, 5))

    for ax, df , pair in zip(axs, compare_dfs,comparisons):
        cat_scatter_genes(ax, df[~(df['is_mito'] | df['is_ribo'] | df['is_lnc'])], 
                        pair[0], pair[1],'black', 
                        label = 'unspecified non-coding', xlim=lim, ylim=lim)
        for col, c, label in zip(['is_pc', 'is_mito','is_ribo','is_lnc'],
                                ['yellow', 'blue', 'red', 'green'],
                                ['protein coding', 'mtRNA', 'rRNA', 'lncRNA']):
            cat_scatter_genes(ax, df[df[col]], pair[0], pair[1], 
                            c, label=label,xlim=lim, ylim=lim)

        show_correlation(ax,df)
        
    axs[0].legend()

    plt.tight_layout()
    plt.show()

# Generate scatter plots to 
# (1) PolyT vs. RandO (2) 10X vs. PolyT 
# (3) 10X vs. RandO (3) 10X vs. Parse
# compare data in the range specified by lim
def compare(datasets, lim):
    data_10x = datasets[0][:,datasets[0].var['percent_counts'] < lim]
    data_polyT = datasets[1][:,datasets[1].var['percent_counts'] < lim]
    data_randO = datasets[2][:,datasets[2].var['percent_counts'] < lim]
    data_parse = datasets[3][:,datasets[3].var['percent_counts'] < lim]

    comparisons = [(data_polyT,data_randO),
                (data_10x,data_polyT),
                (data_10x,data_randO),
                (data_10x,data_parse)]
    compare_names = ['polyT_randO',
                    '10x_polyT',
                    '10x_randO',
                    '10x_parse']

    compare_dfs = []
    for pair in comparisons:
        compare_dfs.append(processing.compare_genes(pair[0],pair[1]))
    
    #compareByDensity(compare_dfs, comparisons, lim)
    #compareByCooks(compare_dfs, comparisons, lim)
    compareByLength(compare_dfs, comparisons, lim)
    compareByGC(compare_dfs, comparisons, lim)
    compareByType(compare_dfs, comparisons, lim)

    return compare_names, compare_dfs

def generate_upset(datasets: ad.AnnData, gene_info: pd.DataFrame, cell_thresh: int = 10, n_top_genes: int = -1) -> pd.DataFrame:
    contents = pd.DataFrame()

    # Generate upset plot comparing genes expressed in each method
    for data in datasets:
        var = data.var[data.var['n_cells']>cell_thresh]
        if n_top_genes > 0:
            top_data = var.sort_values(by='percent_counts',ascending=False).head(1000)
            mask = (gene_info['gene_id'].isin(top_data['gene_id'])).tolist()
        else:
            mask = (gene_info['gene_id'].isin(var['gene_id'])).tolist()
        contents[str(data.uns['name'])] = mask

    upset = Upset.generate_plot(contents)
    contents['gene_ids'] = gene_info['gene_id']
    upset.show()
    return contents

def plot_geneset_metrics(gene_sets: list[pd.Series], set_names: list[str], gene_info: pd.DataFrame) -> None:
    fig, axs = plt.subplots(2, len(gene_sets), figsize=(12,10), sharey='row')
    for ax, data, name in zip(axs[0,:], gene_sets, set_names):
        lengths = gene_info[gene_info['gene_id'].isin(data)]['gene_length'].tolist()    
        ax.violinplot(np.log10(lengths),showextrema=False, showmedians=True)
        ax.set_title(name, fontsize=10)
        ax.set_xticks([])

    for ax, data, name in zip(axs[1,:], gene_sets, set_names):
        gcs = gene_info[gene_info['gene_id'].isin(data)]['gc_content'].tolist()
        ax.violinplot(gcs,showextrema=False, showmedians=True)
        ax.set_title(name, fontsize=10)
        ax.set_xticks([])

    axs[0,0].set_ylabel('Gene Length (by order of magnitude)')
    axs[1,0].set_ylabel('Percent GC Content')

    plt.tight_layout()
    plt.show()

def plot_genetype_counts(gene_sets: list[pd.Series], set_names: list[str], cols: list[str], col_names: list[str], color: list[str], gene_info: pd.DataFrame) -> None:

    fig, axs = plt.subplots(1, len(gene_sets), figsize=(25,5))
    for ax, data, name in zip(axs, gene_sets, set_names):
        sums = []
        for col in cols:
            sums.append(gene_info[col][gene_info['gene_id'].isin(data)].sum())
        ax.bar(col_names,sums,color=color)
        ax.set_title(name, fontsize=10)

    axs[0].set_ylabel('Number of Genes In Intersection')


    plt.tight_layout()
    plt.show()

def plot_genetype_counts(gene_sets: list[pd.Series], set_names: list[str], cols: list[str], col_names: list[str], color: list[str], gene_info: pd.DataFrame) -> None:
    fig, axs = plt.subplots(1, len(gene_sets), figsize=(25,5))
    for ax, data, name in zip(axs, gene_sets, set_names):
        sums = []
        for col in cols:
            sums.append(gene_info[col][gene_info['gene_id'].isin(data)].sum())
        ax.bar(col_names,sums,color=color)
        ax.set_title(name, fontsize=10)

    axs[0].set_ylabel('Number of Genes In Intersection')


    plt.tight_layout()
    plt.show()