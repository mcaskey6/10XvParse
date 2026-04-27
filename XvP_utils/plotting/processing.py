from pathlib import Path
import anndata as ad
import scanpy as sc
import numpy as np
import json
import subprocess
import multiprocessing
from functools import partial
import os
from pybiomart import Server
import scrublet as scr
import statsmodels.api as sm
from scipy.stats import gaussian_kde
import pandas as pd

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

# filter based off of UMI threshold (specified by min_counts)
def refilter(raw_data: ad.AnnData, min_counts: int) -> ad.AnnData:
    data = raw_data.copy()
    sc.pp.filter_cells(data, min_counts=min_counts)
    data.obs['n_genes'] = data.X.astype(bool).sum(axis=1).A1
    data.var['n_cells'] = data.X.astype(bool).sum(axis=0).A1
    data.obs['n_counts'] = data.X.sum(axis=1).A1
    data.var['percent_counts'] = data.X.sum(axis=0).A1/data.X.sum() * 100
    
    return data

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
        output_dir.mkdir(parents=True)

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

    # mitochondrial genes, "MT-" for human, "mt-" for mouse
    data.var["is_mito"] = data.var_names.str.startswith(("mt", "MT"))
    # ribosomal genes, all caps for human, first-letter capitalized for mouse
    data.var["is_ribo"] = data.var_names.str.startswith(("Rps", "Rpl", "RPS", "RPL"))

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
    subset_cols = [col for col in gene_info.columns if col.endswith('_n_cells')]
    gene_info = gene_info.dropna(subset = subset_cols, how = 'all')
    gene_info = gene_info.fillna(0)
    file_path = path / 'gene_data/gene_comparisons.csv' 
    gene_info.to_csv(file_path)
    return gene_info

def detect_doublets(datasets: list[ad.AnnData]) -> scr.Scrublet:
    # Detect doublets Scrublet (by sample is optional)
    def doublet_detection(data):
        scrub = scr.Scrublet(data.X, random_state = 42)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        data.obs['doublet_score'] = doublet_scores
        data.obs['predicted_doublet'] = predicted_doublets
        return data, scrub

    scrubs = []
    for data in datasets:
        data, scrub = doublet_detection(data)
        scrubs.append(scrub)
    return scrubs

# Compares the gene count percentatages between two methods
def compare_genes(data_x, data_y):
    x_var = data_x.var.reset_index(names='gene_name')
    x_var.drop('n_cells', axis=1,inplace=True)

    y_var = data_y.var.reset_index(names='gene_name')
    y_var.drop('n_cells',axis=1,inplace=True)

    shared_data = pd.merge(x_var, y_var, on=['gene_name','gene_id','gene_length', 'gc_content'], how='outer')
    shared_data.fillna(0,inplace=True)

    for col in ['is_lnc','is_mito','is_ribo','is_pc']:
        shared_data[col] = shared_data[col+"_x"] | shared_data[col+"_y"]
        shared_data.drop(col+"_x",axis=1,inplace=True)
        shared_data.drop(col+"_y",axis=1,inplace=True)

    #fit linear regression model
    model = sm.OLS(shared_data['percent_counts_y'], shared_data['percent_counts_x']).fit() 
    np.set_printoptions(suppress=True) # suppress scientific notation
    shared_data['cooks_distance'] = model.get_influence().cooks_distance[0] + 1e-10  # add small value to avoid log(0) issues

    #calculate point density
    xy = np.vstack([shared_data['percent_counts_x'].to_numpy().flatten(), shared_data['percent_counts_y'].to_numpy().flatten()])
    shared_data['point_density'] = gaussian_kde(xy)(xy)

    return shared_data

def mergeByCooks(gene_info, compare_name, compare_df):
    cooks_df = compare_df[['gene_id', 'cooks_distance']]
    gene_info = gene_info.merge(cooks_df, how = 'left', on = ['gene_id'])
    gene_info = gene_info.rename(columns = {'cooks_distance': compare_name + '_distance'})
    return gene_info