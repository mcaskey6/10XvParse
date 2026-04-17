# Import packages
import json
from pathlib import Path
import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import pandas as pd
from pybiomart import Server
from scipy.stats import pearsonr
import multiprocessing, subprocess
from functools import partial
import os

# Load anndata object and begin initial processing
def initProcessing(data_name: str, kb_dir: Path, data_title: str = None, modified: bool = False) -> ad.AnnData:
    m_string = ""
    if modified:
        m_string = "_modified"

    counts_dir = Path.joinpath(kb_dir, f"counts_unfiltered{m_string}")

    # count in h5ad file from kb-python alignment
    data = ad.read_h5ad(Path.joinpath(counts_dir, f"adata.h5ad"))
    
    # Switch gene ids with gene names
    data.var["gene_id"] = data.var.index.tolist()
    gene_names = []
    with open(Path.joinpath(counts_dir, "cells_x_genes.genes.names.txt"), 'r') as file:
        for line in file:
            gene_names.append(line.strip())
    data.var_names = gene_names

    # Add some metadata
    data.obs['n_genes'] = data.X.astype(bool).sum(axis=1).A1
    data.var['n_cells'] = data.X.astype(bool).sum(axis=0).A1
    data.obs['n_counts'] = data.X.sum(axis=1).A1
    data.var['percent_counts'] = data.X.sum(axis=0).A1/data.X.sum() * 100

    # Add unstructured metadata
    data.uns['name'] = data_name
    if data_title:
        data.uns['title'] = data_title
    
    with open(Path.joinpath(kb_dir, "run_info.json"), 'r') as f:
        run_info = json.load(f)
        data.uns['n_processed'] = run_info['n_processed']
        data.uns['n_aligned'] = run_info['n_pseudoaligned']
    data.uns['n_raw_counts'] = data.X.sum() 

    # Just ensure that there are no zero genes or cells
    sc.pp.filter_genes(data, min_cells=1)
    sc.pp.filter_cells(data, min_genes=1)

    return data

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

# Scatter plot of genes detected vs. UMI counts
def scatter_reads(ax, data):
    x = np.asarray(data.X.sum(axis=1))[:,0]
    y = np.asarray(np.sum(data.X>0, axis=1))[:,0]

    ax.scatter(x, y, color="green", alpha=0.25)
    ax.set_xlabel("UMI Counts")
    ax.set_xscale('log')
    ax.set_yscale('log', nonpositive='clip')
    ax.set_title(data.uns['title'] + " Reads")

# filter based off of UMI threshold (specified by min_counts)
def refilter(raw_data: ad.AnnData, min_counts: int) -> ad.AnnData:
    data = raw_data.copy()
    sc.pp.filter_cells(data, min_counts=min_counts)
    data.obs['n_genes'] = data.X.astype(bool).sum(axis=1).A1
    data.var['n_cells'] = data.X.astype(bool).sum(axis=0).A1
    data.obs['n_counts'] = data.X.sum(axis=1).A1
    data.var['percent_counts'] = data.X.sum(axis=0).A1/data.X.sum() * 100
    
    return data

# Knee plot with threshold axes
def knee_plot(ax: matplotlib.axes.Axes, raw_data: ad.AnnData, cutoff: int = 100) -> ad.AnnData:
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
    data = refilter(raw_data, knee[num_cells])
    
    return data

# Scatter plot for mitochondrial percentage
def mito_scatter(ax, data):
    sc.pl.scatter(data, x='n_counts', y='percent_mito', ax=ax, show=False)
    ax.set_title(data.uns['title'])
    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("")

# Generates violin plots of selected cell metadata
def violin_plots(ax_col, data, groups):
    for i, group in enumerate(groups):
        sc.pl.violin(data, group, ax=ax_col[i], stripplot=False, show=False)
        ax_col[i].set_xticks([])

# Find the projected number of unique transcripts to be found using each method if sampling were to continue
# Using PreSeq (Daley and Smith 2013) based off of the recommendation from this paper: https://www.biorxiv.org/content/10.1101/2024.10.09.615408v1.full.pdf

def upsample_helper(data, output_dir: Path) -> None:
    filename = output_dir / f"{data.uns['name']}_counts.txt"
    fileout = output_dir / f"{data.uns['name']}_yield.txt"
    with open(filename, "w") as f:
        counts = data.X.sum(axis=0).A1
        for count in counts:
            f.write(f"{int(count)}\n")

    subprocess.run([
        "preseq", "lc_extrap", "-o", fileout, "-V", filename
    ], check=True)
    return

def upsample(datasets: list[ad.AnnData], output_dir: Path, overwrite: bool = False) -> list[ad.AnnData]:
    if not output_dir.exists():
        output_dir.mkdir(parents=False)

    yield_files_exist = all((output_dir / f"{data.uns['name']}_yield.txt").is_file() for data in datasets)
    if not yield_files_exist or overwrite:
        with multiprocessing.Pool(processes=4) as pool:
            pool.map(partial(upsample_helper, output_dir = output_dir),datasets)
            print("done")
    
    for data in datasets:    
        fileout = output_dir / f"{data.uns['name']}_yield.txt"
        with open(fileout, "r") as f:
            rows = []
            f.readline()
            for line in f:
                rows.append(list(map(float,line.split())))
        data.uns['pred_upsample'] = np.array(rows)
    
    return datasets

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

# Retrieve gene annotations from ensembl
def queryBiotype(dataset):
    # Query Ensembl for biotypes
    type_result = dataset.query(attributes=[
        'external_gene_name',
        'ensembl_gene_id_version',
        'transcript_biotype'])

    # Clean column names
    type_result.columns = ['gene_name', 'gene_id', 'gene_type']

    # From biotypes to bool cols
    type_result['is_lnc'] = (type_result['gene_type'] == 'lncRNA')
    type_result['is_pc'] = (type_result['gene_type'] == 'protein_coding')
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    type_result['is_mito'] = type_result['gene_name'].str.startswith("MT")
    # ribosomal genes
    type_result['is_ribo'] = type_result['gene_name'].str.startswith(("RPS", "RPL"))

    type_result.drop('gene_type', axis=1, inplace=True)
    return type_result

def queryLength(dataset):
    # Query Ensembl for approximate transcript lengths
    length_result = dataset.query(attributes=[
        'external_gene_name',
        'ensembl_gene_id_version', 
        'ensembl_exon_id',
        'exon_chrom_start',
        'exon_chrom_end'])

    # Clean column names
    length_result.columns = ['gene_name','gene_id', 'exon_id', 'exon_start', 'exon_end']
    l_df = length_result.copy()

    # Drop duplicates to avoid counting shared exons multiple times
    l_df.drop_duplicates(subset=['gene_name','gene_id', 'exon_start', 'exon_end'], inplace=True)

    # Compute exon lengths
    l_df['exon_length'] = l_df['exon_end'] - l_df['exon_start'] + 1

    # Sum exon lengths per gene
    gene_lengths = l_df.groupby(['gene_id','gene_name'])['exon_length'].sum().reset_index()
    gene_lengths.rename(columns={'exon_length': 'gene_length'}, inplace=True)

    return gene_lengths

def queryGC(dataset):
    # Query Ensemble for Gene GC Content
    gc_result = dataset.query(attributes=[
    'ensembl_transcript_id_version',
    'external_gene_name',
    'ensembl_gene_id_version',
    'percentage_gene_gc_content'])

    # Clean column names
    gc_result.columns = ['transcript_id', 'gene_name', 'gene_id', 'gc_content']
    # Drop duplicates to avoid counting shared exons multiple times
    gc_result.drop_duplicates(subset=['gene_name','gene_id','transcript_id','gc_content'], inplace=True)
    gc_result.drop('transcript_id', axis=1, inplace=True)

    # average gc_content per trancript
    gc_content = gc_result.groupby(['gene_id','gene_name'])['gc_content'].mean().reset_index()

    return gc_content
    

def queryEnsembl(dir: Path, species: str, overwrite: bool = False) -> None:
    dir = dir / "gene_data"
    if not dir.exists():
        dir.mkdir(parents=False)   
    path = dir / "gene_attributes.csv"

    if os.path.exists(path) and not overwrite:
        gene_info = pd.read_csv(path, index_col = [0])
    else:
        # Connect to server
        server = Server(host='http://ensembl.org')
        dataset = server.marts['ENSEMBL_MART_ENSEMBL'] \
                        .datasets[f'{species}_gene_ensembl']

        type_result = queryBiotype(dataset)
        length_result = queryLength(dataset)
        gc_result = queryGC(dataset)

        gene_info = pd.merge(length_result, type_result, on=['gene_name','gene_id'])
        gene_info = pd.merge(gene_info, gc_result, on=['gene_name','gene_id'])
            
        gene_info.drop_duplicates(subset = ['gene_name', 'gene_id'], inplace=True)
        gene_info.to_csv(path)
        
    return gene_info

# Add cell metrics to anndata object given Ensembl gene metadata
def add_cell_metrics(data, gene_info):
    # retrieve gene metadata
    lnc_result = gene_info["gene_id"][gene_info['is_lnc']].tolist()
    pc_result = gene_info["gene_id"][gene_info['is_pc']].tolist()
    gene_lengths = gene_info[['gene_id', 'gene_length']].drop_duplicates()
    gc_content = gene_info[['gene_id', 'gc_content']].drop_duplicates()

    lncRNA_genes = set(data.var["gene_id"].tolist()).intersection(set(lnc_result))
    pc_genes = set(data.var["gene_id"].tolist()).intersection(set(pc_result))
    
    # Identify lncRNA genes
    data.var["is_lnc"] = np.full(len(data.var_names), False)
    data.var.loc[data.var["gene_id"].isin(list(lncRNA_genes)), ["is_lnc"]] = True

    # Identify protein-coding genes
    data.var["is_pc"] = np.full(len(data.var_names), False)
    data.var.loc[data.var["gene_id"].isin(list(pc_genes)), ["is_pc"]] = True

    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    data.var["is_mito"] = data.var_names.str.startswith("Mt")
    # ribosomal genes, all caps for human, first-letter capitalized for mouse
    data.var["is_ribo"] = data.var_names.str.startswith(("Rps", "Rps"))

    pc_counts = data[:, data.var['is_pc']].X.sum(axis=1)
    mito_counts = data[:, data.var['is_mito']].X.sum(axis=1)
    ribo_counts = data[:, data.var['is_ribo']].X.sum(axis=1)
    lnc_counts = data[:, data.var['is_lnc']].X.sum(axis=1)

    # Calculate total counts per cell
    total_counts = data.X.sum(axis=1)

    # Calculate percent mitochondrial and ribosomal gene expression per cell
    data.obs['percent_pc'] = np.array(pc_counts / total_counts * 100).flatten()
    data.obs['percent_mito'] = np.array(mito_counts / total_counts * 100).flatten()
    data.obs['percent_ribo'] = np.array(ribo_counts / total_counts * 100).flatten()
    data.obs['percent_lnc'] = np.array(lnc_counts / total_counts * 100).flatten()

    # Calculate gene lengths based on exon lengths
    index = data.var.index
    if not 'gene_length' in data.var.columns:
        data.var = data.var.merge(gene_lengths,how='left',on=['gene_id']).fillna(1)
    if not 'gc_content' in data.var.columns:
        data.var = data.var.merge(gc_content,how='left',on=['gene_id']).fillna(0)
    data.var.set_index(index,inplace=True)

def updateGeneInfo(gene_info, datasets, path):
    for data in datasets:
        gene_info = gene_info.merge(data.var[['gene_id','n_cells','percent_counts']], on = ['gene_id'], how = 'left')
        gene_info.rename(columns={'n_cells':data.uns['name']+'_n_cells', 
                                'percent_counts':data.uns['name']+'_percent_counts'}, 
                                inplace=True)
    gene_info = gene_info.dropna(subset = ['10x_n_cells', 'polyT_n_cells', 'randO_n_cells', 'parse_n_cells'], how = 'all')
    gene_info = gene_info.fillna(0)
    file_path = path / 'gene_data/gene_comparisons.csv' 
    gene_info.to_csv(file_path)
    return gene_info

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