from pathlib import Path
from typing import Any
import anndata as ad
import scanpy as sc
import numpy as np
import json
import subprocess
import os
from pybiomart import Server
import scrublet as scr
import statsmodels.api as sm
from scipy.stats import gaussian_kde
import pandas as pd


def init_processing(data_name: str, kb_dir: Path, data_title: str = None, modified: bool = False) -> ad.AnnData:
    """Load a kb-python h5ad output and initialize standard metadata fields.

    Switches gene indices from Ensembl IDs to gene names, computes per-cell
    and per-gene summary statistics, and attaches alignment run info from
    run_info.json.

    Args:
        data_name: Short identifier stored in adata.uns['name'].
        kb_dir: Path to the kb-python output directory containing run_info.json
            and the counts_unfiltered subdirectory.
        data_title: Human-readable label stored in adata.uns['title']. Defaults
            to None (not set).
        modified: If True, loads from the 'counts_unfiltered_modified' subdirectory
            instead of 'counts_unfiltered'.

    Returns:
        AnnData object with gene names as var_names, obs columns 'n_genes' and
        'n_counts', var columns 'gene_id', 'n_cells', and 'percent_counts', and
        uns keys 'name', 'title', 'n_processed', 'n_aligned', and 'n_raw_counts'.
    """
    m_string = ""
    if modified:
        m_string = "_modified"

    counts_dir = Path.joinpath(kb_dir, f"counts_unfiltered{m_string}")

    data = ad.read_h5ad(Path.joinpath(counts_dir, f"adata.h5ad"))

    data.var["gene_id"] = data.var.index.tolist()
    gene_names = []
    with open(Path.joinpath(counts_dir, "cells_x_genes.genes.names.txt"), 'r') as file:
        for line in file:
            gene_names.append(line.strip())
    data.var_names = gene_names

    data.obs['n_genes'] = data.X.astype(bool).sum(axis=1).A1
    data.var['n_cells'] = data.X.astype(bool).sum(axis=0).A1
    data.obs['n_counts'] = data.X.sum(axis=1).A1
    data.var['percent_counts'] = data.X.sum(axis=0).A1 / data.X.sum() * 100

    data.uns['name'] = data_name
    if data_title:
        data.uns['title'] = data_title

    with open(Path.joinpath(kb_dir, "run_info.json"), 'r') as f:
        run_info = json.load(f)
        data.uns['n_processed'] = run_info['n_processed']
        data.uns['n_aligned'] = run_info['n_pseudoaligned']
    data.uns['n_raw_counts'] = data.X.sum()

    sc.pp.filter_genes(data, min_cells=1)
    sc.pp.filter_cells(data, min_genes=1)

    return data


def refilter(raw_data: ad.AnnData, min_counts: int) -> ad.AnnData:
    """Filter cells by a minimum UMI count threshold and recompute metadata.

    Args:
        raw_data: AnnData object to filter (not modified in place).
        min_counts: Minimum number of UMI counts required to retain a cell.

    Returns:
        Filtered copy of raw_data with updated 'n_genes', 'n_counts',
        'n_cells', and 'percent_counts' fields.
    """
    data = raw_data.copy()
    sc.pp.filter_cells(data, min_counts=min_counts)
    data.obs['n_genes'] = data.X.astype(bool).sum(axis=1).A1
    data.var['n_cells'] = data.X.astype(bool).sum(axis=0).A1
    data.obs['n_counts'] = data.X.sum(axis=1).A1
    data.var['percent_counts'] = data.X.sum(axis=0).A1 / data.X.sum() * 100

    return data


def upsample(datasets: list[ad.AnnData], output_dir: Path, overwrite: bool = False) -> None:
    """Launch PreSeq lc_extrap as background processes for each dataset.

    Writes per-gene count files to disk then starts lc_extrap via Popen without
    blocking. Output lands in <output_dir>/<name>_yield.txt when each job
    finishes. Call plot_upsample() to visualise; it loads results lazily from
    disk and skips any dataset whose file is not yet ready.

    Args:
        datasets: List of AnnData objects to process.
        output_dir: Directory for count input and PreSeq output files. Created if absent.
        overwrite: If True, re-runs PreSeq even when output files already exist.
    """
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    launched = []
    for data in datasets:
        fileout = output_dir / f"{data.uns['name']}_yield.txt"
        if fileout.is_file() and not overwrite:
            continue

        filename = output_dir / f"{data.uns['name']}_counts.txt"
        with open(filename, "w") as f:
            counts = data.X.sum(axis=0).A1
            for count in counts:
                f.write(f"{int(count)}\n")

        subprocess.Popen(
            ["preseq", "lc_extrap", "-o", str(fileout), "-V", str(filename)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        launched.append(data.uns['name'])

    if launched:
        print(f"PreSeq running in background for: {', '.join(launched)}")
    else:
        print("PreSeq output already exists for all datasets.")


def _query_biotype(dataset: Any) -> pd.DataFrame:
    """Query Ensembl for gene biotype annotations.

    Retrieves transcript biotypes and converts them to boolean indicator columns
    for lncRNA, protein-coding, mitochondrial, and ribosomal genes.

    Args:
        dataset: pybiomart Dataset object connected to an Ensembl mart.

    Returns:
        DataFrame with columns 'gene_name', 'gene_id', 'is_lnc', 'is_pc',
        'is_mito', and 'is_ribo'.
    """
    type_result = dataset.query(attributes=[
        'external_gene_name',
        'ensembl_gene_id_version',
        'transcript_biotype'])

    type_result.columns = ['gene_name', 'gene_id', 'gene_type']

    type_result['is_lnc'] = (type_result['gene_type'] == 'lncRNA')
    type_result['is_pc'] = (type_result['gene_type'] == 'protein_coding')
    type_result['is_mito'] = type_result['gene_name'].str.startswith("MT")
    type_result['is_ribo'] = type_result['gene_name'].str.startswith(("RPS", "RPL"))

    type_result.drop('gene_type', axis=1, inplace=True)
    return type_result


def _query_length(dataset: Any) -> pd.DataFrame:
    """Query Ensembl for approximate gene lengths based on summed exon coordinates.

    Deduplicates shared exons before summing to avoid double-counting.

    Args:
        dataset: pybiomart Dataset object connected to an Ensembl mart.

    Returns:
        DataFrame with columns 'gene_id', 'gene_name', and 'gene_length' (bp).
    """
    length_result = dataset.query(attributes=[
        'external_gene_name',
        'ensembl_gene_id_version',
        'ensembl_exon_id',
        'exon_chrom_start',
        'exon_chrom_end'])

    length_result.columns = ['gene_name', 'gene_id', 'exon_id', 'exon_start', 'exon_end']
    l_df = length_result.copy()

    l_df.drop_duplicates(subset=['gene_name', 'gene_id', 'exon_start', 'exon_end'], inplace=True)
    l_df['exon_length'] = l_df['exon_end'] - l_df['exon_start'] + 1

    gene_lengths = l_df.groupby(['gene_id', 'gene_name'])['exon_length'].sum().reset_index()
    gene_lengths.rename(columns={'exon_length': 'gene_length'}, inplace=True)

    return gene_lengths


def _query_gc(dataset: Any) -> pd.DataFrame:
    """Query Ensembl for mean GC content percentage per gene.

    Averages across transcripts after deduplicating per-transcript entries.

    Args:
        dataset: pybiomart Dataset object connected to an Ensembl mart.

    Returns:
        DataFrame with columns 'gene_id', 'gene_name', and 'gc_content' (percent).
    """
    gc_result = dataset.query(attributes=[
        'ensembl_transcript_id_version',
        'external_gene_name',
        'ensembl_gene_id_version',
        'percentage_gene_gc_content'])

    gc_result.columns = ['transcript_id', 'gene_name', 'gene_id', 'gc_content']
    gc_result.drop_duplicates(subset=['gene_name', 'gene_id', 'transcript_id', 'gc_content'], inplace=True)
    gc_result.drop('transcript_id', axis=1, inplace=True)

    gc_content = gc_result.groupby(['gene_id', 'gene_name'])['gc_content'].mean().reset_index()

    return gc_content


def query_ensembl(dir: Path, species: str, overwrite: bool = False) -> pd.DataFrame:
    """Retrieve and cache gene annotations (biotype, length, GC content) from Ensembl.

    If a cached CSV already exists at <dir>/gene_data/gene_attributes.csv and
    overwrite is False, the cached file is returned without querying the server.

    Args:
        dir: Base directory under which gene_data/ is created.
        species: Ensembl species prefix used to select the gene dataset
            (e.g. 'hsapiens' or 'mmusculus').
        overwrite: If True, re-queries Ensembl even when the cache file exists.

    Returns:
        DataFrame with columns 'gene_id', 'gene_name', 'gene_length', 'gc_content',
        'is_lnc', 'is_pc', 'is_mito', and 'is_ribo'.
    """
    dir = dir / "gene_data"
    if not dir.exists():
        dir.mkdir(parents=False)
    path = dir / "gene_attributes.csv"

    if os.path.exists(path) and not overwrite:
        gene_info = pd.read_csv(path, index_col=[0])
    else:
        server = Server(host='http://ensembl.org')
        dataset = server.marts['ENSEMBL_MART_ENSEMBL'] \
                        .datasets[f'{species}_gene_ensembl']

        type_result = _query_biotype(dataset)
        length_result = _query_length(dataset)
        gc_result = _query_gc(dataset)

        gene_info = pd.merge(length_result, type_result, on=['gene_name', 'gene_id'])
        gene_info = pd.merge(gene_info, gc_result, on=['gene_name', 'gene_id'])

        gene_info.drop_duplicates(subset=['gene_name', 'gene_id'], inplace=True)
        gene_info.to_csv(path)

    return gene_info


def add_cell_metrics(data: ad.AnnData, gene_info: pd.DataFrame) -> None:
    """Add per-cell biotype percentage metrics and per-gene length/GC annotations.

    Computes percent_pc, percent_mito, percent_ribo, and percent_lnc for each
    cell and merges gene_length and gc_content into data.var if not already present.

    Args:
        data: AnnData object to annotate in place. var must contain 'gene_id'.
        gene_info: DataFrame from query_ensembl with 'gene_id', 'is_lnc', 'is_pc',
            'gene_length', and 'gc_content' columns.
    """
    lnc_result = gene_info["gene_id"][gene_info['is_lnc']].tolist()
    pc_result = gene_info["gene_id"][gene_info['is_pc']].tolist()
    gene_lengths = gene_info[['gene_id', 'gene_length']].drop_duplicates()
    gc_content = gene_info[['gene_id', 'gc_content']].drop_duplicates()

    lncRNA_genes = set(data.var["gene_id"].tolist()).intersection(set(lnc_result))
    pc_genes = set(data.var["gene_id"].tolist()).intersection(set(pc_result))

    data.var["is_lnc"] = np.full(len(data.var_names), False)
    data.var.loc[data.var["gene_id"].isin(list(lncRNA_genes)), ["is_lnc"]] = True

    data.var["is_pc"] = np.full(len(data.var_names), False)
    data.var.loc[data.var["gene_id"].isin(list(pc_genes)), ["is_pc"]] = True

    data.var["is_mito"] = data.var_names.str.startswith(("mt", "MT"))
    data.var["is_ribo"] = data.var_names.str.startswith(("Rps", "Rpl", "RPS", "RPL"))

    pc_counts = data[:, data.var['is_pc']].X.sum(axis=1)
    mito_counts = data[:, data.var['is_mito']].X.sum(axis=1)
    ribo_counts = data[:, data.var['is_ribo']].X.sum(axis=1)
    lnc_counts = data[:, data.var['is_lnc']].X.sum(axis=1)

    total_counts = data.X.sum(axis=1)

    data.obs['percent_pc'] = np.array(pc_counts / total_counts * 100).flatten()
    data.obs['percent_mito'] = np.array(mito_counts / total_counts * 100).flatten()
    data.obs['percent_ribo'] = np.array(ribo_counts / total_counts * 100).flatten()
    data.obs['percent_lnc'] = np.array(lnc_counts / total_counts * 100).flatten()

    index = data.var.index
    if 'gene_length' not in data.var.columns:
        data.var = data.var.merge(gene_lengths, how='left', on=['gene_id']).fillna(1)
    if 'gc_content' not in data.var.columns:
        data.var = data.var.merge(gc_content, how='left', on=['gene_id']).fillna(0)
    data.var.set_index(index, inplace=True)


def update_gene_info(gene_info: pd.DataFrame, datasets: list[ad.AnnData], path: Path) -> pd.DataFrame:
    """Merge per-dataset n_cells and percent_counts into gene_info and save to CSV.

    Args:
        gene_info: Base gene annotation DataFrame from query_ensembl.
        datasets: List of AnnData objects whose var contains 'gene_id', 'n_cells',
            and 'percent_counts'.
        path: Base directory under which gene_data/gene_comparisons.csv is written.

    Returns:
        Updated gene_info DataFrame with per-dataset n_cells and percent_counts
        columns appended, rows with all-NaN n_cells dropped, and remaining NaNs
        filled with 0.
    """
    for data in datasets:
        gene_info = gene_info.merge(data.var[['gene_id', 'n_cells', 'percent_counts']], on=['gene_id'], how='left')
        gene_info.rename(columns={'n_cells': data.uns['name'] + '_n_cells',
                                  'percent_counts': data.uns['name'] + '_percent_counts'},
                         inplace=True)
    subset_cols = [col for col in gene_info.columns if col.endswith('_n_cells')]
    gene_info = gene_info.dropna(subset=subset_cols, how='all')
    gene_info = gene_info.fillna(0)
    file_path = path / 'gene_data/gene_comparisons.csv'
    gene_info.to_csv(file_path)
    return gene_info


def detect_doublets(datasets: list[ad.AnnData]) -> list[scr.Scrublet]:
    """Score and flag doublets in each dataset using Scrublet.

    Results are stored in each AnnData's obs as 'doublet_score' and
    'predicted_doublet'.

    Args:
        datasets: List of AnnData objects to process.

    Returns:
        List of fitted Scrublet objects, one per dataset.
    """
    def doublet_detection(data):
        scrub = scr.Scrublet(data.X, random_state=42)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        data.obs['doublet_score'] = doublet_scores
        data.obs['predicted_doublet'] = predicted_doublets
        return data, scrub

    scrubs = []
    for data in datasets:
        data, scrub = doublet_detection(data)
        scrubs.append(scrub)
    return scrubs


def compare_genes(data_x: ad.AnnData, data_y: ad.AnnData) -> pd.DataFrame:
    """Compare gene percent counts between two datasets and compute Cook's distance.

    Merges var DataFrames on shared gene identifiers, fits an OLS regression of
    y percent counts on x percent counts, and calculates Cook's distance and
    point density for each gene.

    Args:
        data_x: AnnData object for the x-axis dataset.
        data_y: AnnData object for the y-axis dataset.

    Returns:
        DataFrame with columns 'percent_counts_x', 'percent_counts_y',
        'cooks_distance', 'point_density', and gene metadata columns. Genes
        present in only one dataset are included with 0 for the missing values.
    """
    x_var = data_x.var.reset_index(names='gene_name')
    x_var.drop('n_cells', axis=1, inplace=True)

    y_var = data_y.var.reset_index(names='gene_name')
    y_var.drop('n_cells', axis=1, inplace=True)

    shared_data = pd.merge(x_var, y_var, on=['gene_name', 'gene_id', 'gene_length', 'gc_content'], how='outer')
    shared_data.fillna(0, inplace=True)

    for col in ['is_lnc', 'is_mito', 'is_ribo', 'is_pc']:
        shared_data[col] = shared_data[col + "_x"] | shared_data[col + "_y"]
        shared_data.drop(col + "_x", axis=1, inplace=True)
        shared_data.drop(col + "_y", axis=1, inplace=True)

    model = sm.OLS(shared_data['percent_counts_y'], shared_data['percent_counts_x']).fit()
    np.set_printoptions(suppress=True)
    shared_data['cooks_distance'] = model.get_influence().cooks_distance[0] + 1e-10

    xy = np.vstack([shared_data['percent_counts_x'].to_numpy().flatten(), shared_data['percent_counts_y'].to_numpy().flatten()])
    shared_data['point_density'] = gaussian_kde(xy)(xy)

    return shared_data


def merge_by_cooks(gene_info: pd.DataFrame, compare_name: str, compare_df: pd.DataFrame) -> pd.DataFrame:
    """Merge Cook's distances from a comparison DataFrame into gene_info.

    Args:
        gene_info: Base gene annotation DataFrame from query_ensembl.
        compare_name: Label prepended to the resulting column name
            (e.g. '10x_polyT' → '10x_polyT_distance').
        compare_df: Output of compare_genes containing 'gene_id' and
            'cooks_distance' columns.

    Returns:
        Updated gene_info with a new '<compare_name>_distance' column.
    """
    cooks_df = compare_df[['gene_id', 'cooks_distance']]
    gene_info = gene_info.merge(cooks_df, how='left', on=['gene_id'])
    gene_info = gene_info.rename(columns={'cooks_distance': compare_name + '_distance'})
    return gene_info
