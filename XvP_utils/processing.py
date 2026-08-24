from pathlib import Path
import re
import warnings
import anndata as ad
import scanpy as sc
import numpy as np
import json
import subprocess
import os
from Bio import SeqIO
import scrublet as scr
import statsmodels.api as sm
from scipy.stats import gaussian_kde, entropy
import pandas as pd
import urllib.request
import urllib.parse
import io
import scclr

_ORTHOLOGS: dict[str, str] | None = None

def _load_orthologs(cache_dir: Path | None = None) -> dict[str, str]:
    global _ORTHOLOGS
    if _ORTHOLOGS is not None:
        return _ORTHOLOGS

    if cache_dir is None:
        cache_dir = Path(__file__).parent
    cache_path = cache_dir / "orthologs.csv"

    if cache_path.exists():
        df = pd.read_csv(cache_path)
    else:
        xml = ('<?xml version="1.0" encoding="UTF-8"?>'
               '<!DOCTYPE Query>'
               '<Query virtualSchemaName="default" formatter="TSV" header="1"'
               ' uniqueRows="1" count="" datasetConfigVersion="0.6">'
               '<Dataset name="mmusculus_gene_ensembl" interface="default">'
               '<Attribute name="external_gene_name"/>'
               '<Attribute name="hsapiens_homolog_associated_gene_name"/>'
               '</Dataset></Query>')
        url = "https://www.ensembl.org/biomart/martservice?query=" + urllib.parse.quote(xml)
        with urllib.request.urlopen(url, timeout=60) as resp:
            raw = resp.read().decode("utf-8")
        df = pd.read_csv(io.StringIO(raw), sep="\t")
        df.columns = ["mouse_symbol", "human_symbol"]
        df = df.dropna(subset=["mouse_symbol", "human_symbol"])
        df = df[(df["mouse_symbol"] != "") & (df["human_symbol"] != "")]
        df = df.drop_duplicates()
        df.to_csv(cache_path, index=False)

    _ORTHOLOGS = dict(zip(df["mouse_symbol"], df["human_symbol"]))
    return _ORTHOLOGS


def _normalize_gene_name(name: str, orthologs: dict[str, str]) -> str:
    for prefix in ("HUMAN_", "MOUSE_"):
        if name.startswith(prefix):
            name = name[len(prefix):]
            break
    name = re.sub(r"\s*\(.*$", "", name).strip()
    if name in orthologs:
        return orthologs[name].upper()
    return name.upper()


def init_processing(data_name: str, assay: str, project_dir: str, analysis_name: str, data_title: str = None, type: str = "", modified: bool = False, sampled: bool = True) -> ad.AnnData:
    """Load a kb-python h5ad output and initialize standard metadata fields.

    Switches gene indices from Ensembl IDs to gene names, computes per-cell
    and per-gene summary statistics, and attaches alignment run info from
    kb_python's run_info.json and STARsolo's Log.final.out.

    Args:
        data_name: Short identifier stored in adata.uns['name'].
        project_dir: Base directory containing the Data/ subdirectory with kb-python and STAR outputs.
        analysis_name: Subdirectory under Data/ containing the assay-specific outputs.
        data_title: Optional longer title stored in adata.uns['title'] for display in plots
        sample: Optional sample identifier to distinguish different runs of the same analysis (e.g. "H2"). 
                If provided, looks for data under assay+sample (e.g. "10x_H2") and appends sample to the 
                name and title.
        type: Optional string to distinguish different comparison within the same analysis (e.g. "standard", "mini")
              If provided, looks for data under assay/type (e.g. "10x/standard").
        modified: If True, looks for data under counts_unfiltered_modified/ instead of counts_unfiltered/.
        sampled: If True, looks for data under sampled_{data_name}_out/ instead of {data_name}_out/. 
                 This is used to distinguish the downsampled datasets from the full datasets.
        
    Returns:
        AnnData object with gene names as var_names, obs columns 'n_genes' and
        'n_counts', var columns 'gene_id', 'n_cells', and 'percent_counts', and
        uns keys 'name', 'title', 'n_processed', 'n_aligned', and 'n_raw_counts'.
    """
    def get_kb_dir(project_dir: Path, data_name: str, assay: str, type: str = None) -> Path:
        type_string = f"_{type}" if type else ""
        if sampled:
            kb_dir = project_dir / "Data" / analysis_name / assay / "kb_python" / f"sampled_{data_name}{type_string}_out"
        else:
            kb_dir = project_dir / "Data" / analysis_name / assay / "kb_python" / f"{data_name}_out"
        return kb_dir
    
    def get_star_dir(project_dir: Path, assay: str, type: str = None) -> Path:
        if type is None:
            star_dir = project_dir / "Data" / analysis_name / assay / "STARsolo"
        else:
            star_dir = project_dir / "Data" / analysis_name / assay / "STARsolo" / type
        return star_dir

    m_string = ""
    if modified:
        m_string = "_modified"
    
    kb_dir = get_kb_dir(Path(project_dir), data_name, assay, type)
    star_dir = get_star_dir(Path(project_dir), assay, type)


    counts_dir = Path.joinpath(kb_dir, f"counts_unfiltered{m_string}")

    data = ad.read_h5ad(Path.joinpath(counts_dir, f"adata.h5ad"))

    data.var["gene_id"] = data.var.index.tolist()
    gene_names = []
    with open(Path.joinpath(counts_dir, "cells_x_genes.genes.names.txt"), 'r') as file:
        for line in file:
            gene_names.append(line.strip())
    data.var_names = gene_names

    sc.pp.filter_genes(data, min_cells=1)
    sc.pp.filter_cells(data, min_genes=1)

    data.uns['name'] = data_name
    if data_title:
        data.uns['title'] = data_title

    with open(Path.joinpath(kb_dir, "run_info.json"), 'r') as f:
        kb_run_info = json.load(f)
    data.uns['n_processed'] = kb_run_info['n_processed']
    data.uns['n_aligned'] = kb_run_info['n_pseudoaligned']
    data.uns['n_unique'] = kb_run_info['n_unique']
    data.uns['n_raw_counts'] = data.X.sum()

    with open(Path.joinpath(star_dir, f"{data_name}_Log.final.out"), 'r') as f:
        star_dict = {}
        for line in f:
            line = line.strip()
            if not "|" in line:
                continue
            key, value = [part.strip() for part in line.split("|", 1)]
            star_dict[key] = value
    data.uns['star_uniquely_mapped'] = int(star_dict.get("Uniquely mapped reads number"))
    data.uns['star_multimapped'] = int(star_dict.get("Number of reads mapped to multiple loci")) \
                                    + int(star_dict.get("Number of reads mapped to too many loci"))
    data.uns['star_unmapped'] = int(star_dict.get("Number of input reads")) \
                                - int(data.uns['star_uniquely_mapped']) \
                                - int(data.uns['star_multimapped'])
    
    data.obs['n_genes'] = data.X.astype(bool).sum(axis=1).A1
    data.var['n_cells'] = data.X.astype(bool).sum(axis=0).A1
    data.obs['n_counts'] = data.X.sum(axis=1).A1
    data.var['percent_counts'] = data.X.sum(axis=0).A1 / data.uns["n_raw_counts"] * 100
    data.uns['n_genes'] = data.X.astype(bool).sum()

    data.uns['percent_nascent'] = data.layers['nascent'].sum() / data.uns["n_raw_counts"] * 100 
    data.uns['percent_ambiguous'] = data.layers['ambiguous'].sum() / data.uns["n_raw_counts"] * 100
    data.uns['percent_mature'] = data.layers['mature'].sum() / data.uns["n_raw_counts"] * 100
    data.var['percent_nascent'] = data.layers['nascent'].sum(axis=0).A1 / data.X.sum(axis=0).A1 * 100
    data.obs['percent_nascent'] = data.layers['nascent'].sum(axis=1).A1 / data.obs['n_counts'] * 100

    return data


def refilter(raw_data: ad.AnnData, min_counts: int, transform: bool = False) -> ad.AnnData:
    """Filter cells by a minimum UMI count threshold and recompute metadata.

    Args:
        raw_data: AnnData object to filter (not modified in place).
        min_counts: Minimum number of UMI counts required to retain a cell.
        transform: If True, recompute the X matrix as log-CPM (log1p of counts per million)
                   after filtering.

    Returns:
        Filtered copy of raw_data with updated 'n_genes', 'n_counts',
        'n_cells', and 'percent_counts' fields.
    """
    data = raw_data.copy()
    sc.pp.filter_cells(data, min_counts=min_counts)
    sc.pp.filter_genes(data, min_cells=1)

    data.uns['n_raw_counts_filtered'] = data.X.sum()
    data.uns['n_genes'] = data.X.astype(bool).sum()

    data.uns['percent_nascent'] = data.layers['nascent'].sum() / data.uns["n_raw_counts_filtered"] * 100 
    data.uns['percent_ambiguous'] = data.layers['ambiguous'].sum() / data.uns["n_raw_counts_filtered"] * 100
    data.uns['percent_mature'] = data.layers['mature'].sum() / data.uns["n_raw_counts_filtered"] * 100
    data.var['percent_nascent'] = data.layers['nascent'].sum(axis=0).A1 / data.X.sum(axis=0).A1 * 100
    data.obs['percent_nascent'] = data.layers['nascent'].sum(axis=1).A1 / data.X.sum(axis=1).A1 * 100
    data.var['percent_counts'] = data.X.sum(axis=0).A1 / data.uns['n_raw_counts_filtered'] 
    
    data.var['normalized_counts'] = np.log((data.X.sum(axis=0).A1 / data.uns["n_raw_counts_filtered"] * 1e6) + 1)

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


def _warn_unmatched(left_names: set, right_names: set,
                    caller: str, left_label: str, right_label: str) -> None:
    """Warn about gene names present in one set but not the other."""
    missing = left_names - right_names
    if missing:
        names = sorted(missing)
        warnings.warn(
            f"{caller}: {len(missing)} genes in {left_label} "
            f"but missing from {right_label}:\n  {names}"
        )


def _biotype_from_gtf(gtf_path: Path) -> pd.DataFrame:
    """Extract gene biotype annotations from a GTF file.

    Reads gene-level rows only and converts gene_biotype to boolean indicator
    columns for lncRNA, protein-coding, mitochondrial, and ribosomal genes.

    Args:
        gtf_path: Path to the GTF file (plain text or .gz not supported here).

    Returns:
        DataFrame with columns 'gene_name', 'gene_id', 'is_lnc', 'is_pc',
        'is_mito', and 'is_ribo'. gene_id includes the version suffix (e.g.
        ENSMUSG00000104478.2) to match the t2g.txt format.
    """
    _attr = re.compile(r'(\w+) "([^"]+)"')
    records = []
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[2] != 'gene':
                continue
            attrs = dict(_attr.findall(fields[8]))
            gene_id = attrs.get('gene_id', '')
            version = attrs.get('gene_version', '')
            gene_name = attrs.get('gene_name', gene_id)
            biotype = attrs.get('gene_biotype', '')
            if version:
                gene_id = f"{gene_id}.{version}"
            records.append({'gene_name': gene_name, 'gene_id': gene_id, 'gene_biotype': biotype})

    df = pd.DataFrame(records)
    df['is_lnc'] = df['gene_biotype'] == 'lncRNA'
    df['is_pc'] = df['gene_biotype'] == 'protein_coding'
    df['is_mito'] = df['gene_name'].str.startswith(("mt-", "MT-", "MOUSE_mt-", "HUMAN_MT-"))
    df['is_ribo'] = df['gene_name'].str.startswith((
        "Rps", "Rpl", "Mrps", "Mrpl",
        "RPS", "RPL", "MRPS", "MRPL",
        "MOUSE_Rps", "MOUSE_Rpl", "MOUSE_Mrps", "MOUSE_Mrpl",
        "HUMAN_RPS", "HUMAN_RPL", "HUMAN_MRPS", "HUMAN_MRPL",
    ))
    df['is_oxphos'] = df['gene_name'].str.startswith((
        "Sdh",  "Uqcr", "Cox",  "Nduf", "Atp5",
        "SDH",  "UQCR", "COX",  "NDUF", "ATP5",
        "MOUSE_Sdh",  "MOUSE_Uqcr", "MOUSE_Cox",  "MOUSE_Nduf", "MOUSE_Atp5",
        "HUMAN_SDH",  "HUMAN_UQCR", "HUMAN_COX",  "HUMAN_NDUF", "HUMAN_ATP5",
        # "mt-Nd", "mt-Co", "mt-Atp", "mt-Cytb",
        # "MT-ND", "MT-CO", "MT-ATP", "MT-CYB",
        # "MOUSE_mt-Nd", "MOUSE_mt-Co", "MOUSE_mt-Atp", "MOUSE_mt-Cytb",
        # "HUMAN_MT-ND", "HUMAN_MT-CO", "HUMAN_MT-ATP", "HUMAN_MT-CYB",
    ))
    df['is_pseudo'] = df['gene_biotype'].str.contains('pseudogene', case=False)
    df.drop('gene_biotype', axis=1, inplace=True)
    return df


def annotate_gene_set(gene_info: pd.DataFrame, csv_path: Path,
                      gene_set: str, col_name: str,
                      term: str | None = None) -> None:
    """Add a boolean column to gene_info by matching against a named gene set CSV.

    Reads all rows in the CSV whose Gene_set column equals gene_set, collects
    every gene name across those rows, and matches case-insensitively (stripping
    HUMAN_/MOUSE_ prefixes) against gene_info['gene_name'].

    Args:
        gene_info: DataFrame to annotate in place. Must contain 'gene_name'.
        csv_path: Path to the gene sets CSV (columns: Gene_set, Term, Genes
            with semicolon-separated names).
        gene_set: Value to match in the Gene_set column (e.g. "Transcription Factors").
        col_name: Name of the boolean column to add (e.g. "is_tf").
        term: If provided, further restrict to rows whose Term column equals
            this value (e.g. "Human" or "Mouse").
    """
    def _norm(name: str) -> str:
        for prefix in ("HUMAN_", "MOUSE_"):
            if name.startswith(prefix):
                return name[len(prefix):].upper()
        return name.upper()

    enr_df = pd.read_csv(csv_path, index_col=0)
    enr_df = enr_df[enr_df['Gene_set'] == gene_set]
    if term is not None:
        enr_df = enr_df[enr_df['Term'] == term]

    ref_names: set[str] = set()
    for _, row in enr_df.iterrows():
        ref_names.update(_norm(g) for g in row['Genes'].split(';'))

    norm_info = gene_info['gene_name'].apply(_norm)
    gene_info[col_name] = norm_info.isin(ref_names)

    label = f'"{gene_set}"' + (f' / "{term}"' if term else "")
    _warn_unmatched(ref_names, set(norm_info), "annotate_gene_set", label, "gene_info")

    n = gene_info[col_name].sum()
    print(f"annotate_gene_set: {n}/{len(ref_names)} {label} genes matched → '{col_name}'")


def _query_from_fasta(cdna_fasta: Path, t2g: Path) -> pd.DataFrame:
    """Compute per-gene median transcript length and GC content from a kb-python cDNA FASTA.

    Parses every transcript in the cDNA FASTA, maps transcripts to genes via the
    t2g file, and takes the median spliced length and GC content across all
    isoforms of each gene. This is more complete and accurate than querying
    Ensembl because every gene present in the kallisto index is guaranteed to
    have an entry.

    Gene names that are themselves Ensembl gene IDs (e.g. unnamed genes) are
    normalised to the gene_id so that the (gene_id, gene_name) grouping is stable.

    Args:
        cdna_fasta: Path to the kb-python cdna.fasta file used to build the kallisto index.
        t2g: Path to the transcript-to-gene mapping file (t2g.txt).

    Returns:
        DataFrame with columns 'gene_id', 'gene_name', 'gene_length', and 'gc_content'
        (one row per gene).
    """
    tx_records = {}
    for record in SeqIO.parse(cdna_fasta, "fasta"):
        seq = str(record.seq).upper()
        length = len(seq)
        gc = (seq.count('G') + seq.count('C')) / length * 100 if length > 0 else 0.0
        tx_records[record.id] = {"length": length, "gc_content": gc}

    t2g_df = pd.read_csv(t2g, sep="\t", header=None,
                         usecols=[0, 1, 2], names=["transcript_id", "gene_id", "gene_name"])
    t2g_df["gene_name"] = t2g_df["gene_name"].fillna(t2g_df["gene_id"])
    ens_mask = t2g_df["gene_name"].str.match(r'^ENS[A-Z]*G\d')
    t2g_df.loc[ens_mask, "gene_name"] = t2g_df.loc[ens_mask, "gene_id"]
    t2g_df["gene_length"] = t2g_df["transcript_id"].map(
        lambda x: tx_records.get(x, {}).get("length"))
    t2g_df["gc_content"] = t2g_df["transcript_id"].map(
        lambda x: tx_records.get(x, {}).get("gc_content"))

    return t2g_df.groupby(["gene_id", "gene_name"])[["gene_length", "gc_content"]].median().reset_index()


def query_ensembl(project_dir: Path, index_dir: Path, species: str, overwrite: bool = False) -> pd.DataFrame:
    """Retrieve and cache gene annotations (biotype, length, GC content).

    Biotype is read from the GTF file in index_dir. Length and GC content are
    computed directly from the kb-python cDNA FASTA, guaranteeing complete
    coverage for every gene in the kallisto index.

    Results are cached per species at <project_dir>/Notebooks/gene_info/<species>/gene_attributes.csv.

    Args:
        project_dir: Root project directory.
        index_dir: Path to the kb-python kallisto index directory containing
            ref.gtf, cdna.fasta, and t2g.txt.
        species: Species identifier (e.g. "human", "mouse") used for cache path.
        overwrite: If True, re-reads even when the cache file exists.

    Returns:
        DataFrame with columns 'gene_id', 'gene_name', 'gene_length', 'gc_content',
        'is_lnc', 'is_pc', 'is_mito', and 'is_ribo'.
    """
    cache_dir = project_dir / "Notebooks" / "gene_info" / species
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / "gene_attributes.csv"

    if os.path.exists(path) and not overwrite:
        gene_info = pd.read_csv(path, index_col=[0])
    else:
        type_result = _biotype_from_gtf(index_dir / "ref.gtf")
        fasta_result = _query_from_fasta(index_dir / "cdna.fasta", index_dir / "t2g.txt")

        gene_info = pd.merge(fasta_result, type_result.drop(columns='gene_name'),
                             on='gene_id', how='left')
        bool_cols = gene_info.select_dtypes(include='bool').columns
        for col in type_result.columns.difference(['gene_id', 'gene_name']):
            if col not in bool_cols:
                continue
            gene_info[col] = gene_info[col].fillna(False)
        gene_info.drop_duplicates(subset='gene_id', inplace=True)

        fasta_ids = set(fasta_result['gene_id'])
        gtf_ids = set(type_result['gene_id'])
        _warn_unmatched(fasta_ids, gtf_ids, "query_ensembl", "cDNA FASTA", "GTF")

        gene_info.to_csv(path)

    return gene_info


def query_ensembl_combined(project_dir: Path, index_dir: Path, species: str, overwrite: bool = False) -> pd.DataFrame:
    """Retrieve and cache gene annotations for a combined human/mouse reference.

    Length and GC content are computed from the combined kb-python cDNA FASTA
    (which already contains HUMAN_/MOUSE_ prefixes). Biotype annotations are
    read from the combined ref.gtf, which also carries HUMAN_/MOUSE_ prefixes.

    Results are cached per species at <project_dir>/Notebooks/gene_info/<species>/gene_attributes.csv.

    Args:
        project_dir: Root project directory.
        index_dir: Path to the combined kb-python kallisto index directory containing
            ref.gtf, cdna.fasta, and t2g.txt.
        species: Species identifier (e.g. "human_mouse") used for cache path.
        overwrite: If True, re-reads even when the cache file exists.

    Returns:
        DataFrame with columns 'gene_id', 'gene_name', 'gene_length', 'gc_content',
        'is_lnc', 'is_pc', 'is_mito', and 'is_ribo'.
    """
    cache_dir = project_dir / "Notebooks" / "gene_info" / species
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / "gene_attributes.csv"

    if os.path.exists(path) and not overwrite:
        return pd.read_csv(path, index_col=[0])

    fasta_result = _query_from_fasta(index_dir / "cdna.fasta", index_dir / "t2g.txt")
    type_result = _biotype_from_gtf(index_dir / "ref.gtf")

    combined = pd.merge(fasta_result, type_result.drop(columns='gene_name'),
                        on='gene_id', how='left')
    bool_cols = combined.select_dtypes(include='bool').columns
    for col in type_result.columns.difference(['gene_id', 'gene_name']):
        if col not in bool_cols:
            continue
        combined[col] = combined[col].fillna(False)
    combined.drop_duplicates(subset='gene_id', inplace=True)

    fasta_ids = set(fasta_result['gene_id'])
    gtf_ids = set(type_result['gene_id'])
    _warn_unmatched(fasta_ids, gtf_ids, "query_ensembl_combined", "cDNA FASTA", "GTF")

    combined.to_csv(path)
    return combined


def compute_gene_metrics(
    project_dir: Path,
    index_dir: Path,
    species: str,
    overwrite: bool = False,
) -> pd.DataFrame:
    """Compute three pairs of gene length and GC content metrics and cache to CSV.

    The three methods capture different biological definitions of a gene:
      - Genomic span: distance from gene start to end (including introns).
      - Median transcript: median spliced length and GC across all isoforms.
      - Exon union: total bp in the merged non-redundant exon intervals.

    GC content for genomic span and exon union is extracted from the genome
    FASTA via pysam (0-based half-open coordinates). Median transcript GC
    is derived from the cDNA FASTA.

    Results are cached per species at <project_dir>/Notebooks/gene_info/<species>/gene_metrics_detailed.csv.

    Args:
        project_dir: Root project directory.
        index_dir: Path to the kb-python kallisto index directory containing
            ref.gtf, cdna.fasta, t2g.txt, and ref.fa.
        species: Species identifier (e.g. "human", "mouse") used for cache path.
        overwrite: If True, re-computes even when the cache file exists.

    Returns:
        DataFrame with columns 'gene_id', 'gene_name', 'length_genomic',
        'gc_genomic', 'length_median_tx', 'gc_median_tx', 'length_exon_union',
        and 'gc_exon_union'.
    """
    import pysam

    cache_dir = project_dir / "Notebooks" / "gene_info" / species
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / "gene_metrics_detailed.csv"

    if path.exists() and not overwrite:
        return pd.read_csv(path, index_col=0)

    gtf_path = index_dir / "ref.gtf"
    cdna_fasta = index_dir / "cdna.fasta"
    t2g = index_dir / "t2g.txt"
    fasta_path = index_dir / "ref.fa"

    _attr = re.compile(r'(\w+) "([^"]+)"')
    gene_meta = {}   # (gene_id, gene_name) -> (chrom, start, end)
    exon_lists = {}  # (gene_id, gene_name) -> [(chrom, start, end), ...]

    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            feat = fields[2]
            if feat not in ('gene', 'exon'):
                continue
            attrs = dict(_attr.findall(fields[8]))
            gid = attrs.get('gene_id', '')
            ver = attrs.get('gene_version', '')
            gname = attrs.get('gene_name', gid)
            if ver:
                gid = f"{gid}.{ver}"
            key = (gid, gname)
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            if feat == 'gene':
                gene_meta[key] = (chrom, start, end)
            else:
                exon_lists.setdefault(key, []).append((chrom, start, end))

    # Index genome FASTA if needed
    fai = Path(str(fasta_path) + ".fai")
    if not fai.exists():
        print("Indexing genome FASTA (one-time)…")
        pysam.faidx(str(fasta_path))
    fasta = pysam.FastaFile(str(fasta_path))

    def _gc(seq: str) -> float:
        seq = seq.upper()
        n = len(seq)
        return (seq.count('G') + seq.count('C')) / n * 100 if n > 0 else 0.0

    def _merge_intervals(intervals):
        """Merge a list of (chrom, start, end) 1-based inclusive intervals."""
        sorted_ivs = sorted(intervals, key=lambda x: x[1])
        merged = [sorted_ivs[0]]
        for chrom, s, e in sorted_ivs[1:]:
            if s <= merged[-1][2] + 1:
                merged[-1] = (merged[-1][0], merged[-1][1], max(merged[-1][2], e))
            else:
                merged.append((chrom, s, e))
        return merged

    genomic_rows = []
    exon_union_rows = []

    for key, (chrom, start, end) in gene_meta.items():
        gid, gname = key
        length = end - start + 1
        try:
            seq = fasta.fetch(chrom, start - 1, end)  # pysam: 0-based half-open
            gc = _gc(seq)
        except (KeyError, ValueError):
            gc = float('nan')
        genomic_rows.append({'gene_id': gid, 'gene_name': gname,
                              'length_genomic': length, 'gc_genomic': gc})

        exons = exon_lists.get(key, [(chrom, start, end)])
        merged = _merge_intervals(exons)
        union_len = sum(e - s + 1 for _, s, e in merged)
        seqs = []
        for chrom_e, s, e in merged:
            try:
                seqs.append(fasta.fetch(chrom_e, s - 1, e))
            except (KeyError, ValueError):
                pass
        gc_union = _gc(''.join(seqs)) if seqs else float('nan')
        exon_union_rows.append({'gene_id': gid, 'gene_name': gname,
                                'length_exon_union': union_len, 'gc_exon_union': gc_union})

    fasta.close()

    # Median transcript length and GC from cDNA FASTA
    median_tx = _query_from_fasta(cdna_fasta, t2g).rename(
        columns={"gene_length": "length_median_tx", "gc_content": "gc_median_tx"})

    df_genomic = pd.DataFrame(genomic_rows)
    df_exon = pd.DataFrame(exon_union_rows)

    genomic_exon = df_genomic.drop(columns='gene_name').merge(
        df_exon.drop(columns='gene_name'), on='gene_id')
    result = median_tx.merge(genomic_exon, on='gene_id', how='left')
    result.to_csv(path)
    return result


def add_cell_metrics(data: ad.AnnData, gene_info: pd.DataFrame) -> None:
    """Add per-cell biotype percentage metrics and per-gene length/GC annotations.

    Computes percent_pc, percent_mito, percent_ribo, and percent_lnc for each
    cell and merges gene_length and gc_content into data.var if not already present.

    Args:
        data: AnnData object to annotate in place. var must contain 'gene_id'.
        gene_info: DataFrame from query_ensembl with 'gene_id', 'is_lnc', 'is_pc',
            'gene_length', and 'gc_content' columns.
    """
    data_ids = set(data.var["gene_id"].tolist())
    info_ids = set(gene_info["gene_id"].tolist())
    _warn_unmatched(data_ids, info_ids, "add_cell_metrics", "adata", "gene_info")

    lnc_result = gene_info["gene_id"][gene_info['is_lnc']].tolist()
    pc_result = gene_info["gene_id"][gene_info['is_pc']].tolist()
    mito_result = gene_info["gene_id"][gene_info['is_mito']].tolist()
    ribo_result = gene_info["gene_id"][gene_info['is_ribo']].tolist()
    oxphos_result = gene_info["gene_id"][gene_info['is_oxphos']].tolist()
    pseudo_result = gene_info["gene_id"][gene_info['is_pseudo']].tolist()
    gene_lengths = gene_info[['gene_id', 'gene_length']].drop_duplicates()
    gc_content = gene_info[['gene_id', 'gc_content']].drop_duplicates()

    lncRNA_genes = set(data.var["gene_id"].tolist()).intersection(set(lnc_result))
    pc_genes = set(data.var["gene_id"].tolist()).intersection(set(pc_result))
    mito_genes = set(data.var["gene_id"].tolist()).intersection(set(mito_result))
    ribo_genes = set(data.var["gene_id"].tolist()).intersection(set(ribo_result))
    pseudo_genes = set(data.var["gene_id"].tolist()).intersection(set(pseudo_result))
    oxphos_genes = set(data.var["gene_id"].tolist()).intersection(set(oxphos_result))

    data.var["is_lnc"] = np.full(len(data.var_names), False)
    data.var.loc[data.var["gene_id"].isin(list(lncRNA_genes)), ["is_lnc"]] = True

    data.var["is_pc"] = np.full(len(data.var_names), False)
    data.var.loc[data.var["gene_id"].isin(list(pc_genes)), ["is_pc"]] = True

    data.var["is_mito"] = np.full(len(data.var_names), False)
    data.var.loc[data.var["gene_id"].isin(list(mito_genes)), ["is_mito"]] = True

    data.var["is_ribo"] = np.full(len(data.var_names), False)
    data.var.loc[data.var["gene_id"].isin(list(ribo_genes)), ["is_ribo"]] = True

    data.var["is_pseudo"] = np.full(len(data.var_names), False)
    data.var.loc[data.var["gene_id"].isin(list(pseudo_genes)), ["is_pseudo"]] = True

    data.var["is_oxphos"] = np.full(len(data.var_names), False)
    data.var.loc[data.var["gene_id"].isin(list(oxphos_genes)), ["is_oxphos"]] = True

    
    pc_counts = data[:, data.var['is_pc']].X.sum(axis=1)
    mito_counts = data[:, data.var['is_mito']].X.sum(axis=1)
    ribo_counts = data[:, data.var['is_ribo']].X.sum(axis=1)
    lnc_counts = data[:, data.var['is_lnc']].X.sum(axis=1)
    oxphos_counts = data[:, data.var['is_oxphos']].X.sum(axis=1)

    total_counts = data.X.sum(axis=1)

    data.obs['percent_pc'] = np.array(pc_counts / total_counts * 100).flatten()
    data.obs['percent_mito'] = np.array(mito_counts / total_counts * 100).flatten()
    data.obs['percent_ribo'] = np.array(ribo_counts / total_counts * 100).flatten()
    data.obs['percent_lnc'] = np.array(lnc_counts / total_counts * 100).flatten()
    data.obs['percent_oxphos'] = np.array(oxphos_counts / total_counts * 100).flatten()

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
    gene_info.to_csv(path)
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

# Export bulk counts for H2 to compare with bulk RNA-seq data.
def export_bulk_counts(datasets: list[ad.AnnData], sample: str = None):
    """Aggregate counts across all cells in each dataset and export to a TSV file for comparison with edgeR. 
    Output file is bulk_counts/bulk_counts.tsv, or bulk_counts/<sample>_bulk_counts.tsv if a sample name is provided.

    Args:
        datasets: List of AnnData objects to process. Expects datasets[0] to be 10x and datasets[3] to be Parse.
        sample: Optional sample identifier to append to output file. If None, defaults to "bulk_counts.tsv". 
                If provided, output file is "<sample>_bulk_counts.tsv".
    """

    if sample:
        sample_str = f"_{sample}"
    else:        
        sample_str = ""

    bulk_10x_df = pd.DataFrame({
        "gene_id":   datasets[0].var["gene_id"].values,
        "gene_name": datasets[0].var_names,
        f"tenx":   np.asarray(datasets[0].X.sum(axis=0)).flatten().astype(int)})
    
    bulk_parse_df = pd.DataFrame({
        "gene_id":   datasets[3].var["gene_id"].values,
        "gene_name": datasets[3].var_names,
        f"parse":  np.asarray(datasets[3].X.sum(axis=0)).flatten().astype(int)})
    
    bulk_df = pd.merge(bulk_10x_df, bulk_parse_df, on=["gene_id", "gene_name"], how="outer")
    bulk_df.fillna(0, inplace=True)

    bulk_counts_dir = Path("bulk_counts")
    bulk_counts_dir.mkdir(exist_ok=True)

    if sample:
        bulk_counts_file = bulk_counts_dir / f"{sample}_bulk_counts.tsv"
    else:
        bulk_counts_file = bulk_counts_dir / f"bulk_counts.tsv"

    bulk_df.to_csv(bulk_counts_file, sep="\t", index=False)

    print(f"Exported {len(datasets[0].var)} genes")
    print(f"\ttenx{sample_str} total counts: {bulk_10x_df['tenx'].sum():,}  ({datasets[0].n_obs:,} cells)")
    print(f"\tparse{sample_str} total counts: {bulk_parse_df['parse'].sum():,}  ({datasets[3].n_obs:,} cells)")


def compare_genes(data_x: ad.AnnData, data_y: ad.AnnData, comparison_axis: str="normalized_counts") -> pd.DataFrame:
    """Compare log-normalized bulk gene counts between two datasets and compute Cook's distance.

    Merges var DataFrames on shared gene identifiers, fits an OLS regression of
    y normalized counts on x normalized counts, and calculates Cook's distance and
    point density for each gene.

    Args:
        data_x: AnnData object for the x-axis dataset.
        data_y: AnnData object for the y-axis dataset.
        comparison_axis: var column to compare (default: 'normalized_counts').

    Returns:
        DataFrame with columns "{comparison_axis}_x", "{comparison_axis}_y",
        'cooks_distance', 'point_density', and gene metadata columns. Genes
        present in only one dataset are included with 0 for the missing values.
    """
    x_var = data_x.var.reset_index(names='gene_name')
    x_var.drop('n_cells', axis=1, inplace=True)

    y_var = data_y.var.reset_index(names='gene_name')
    y_var.drop('n_cells', axis=1, inplace=True)

    shared_data = pd.merge(x_var, y_var, on=['gene_name', 'gene_id', 'gene_length', 'gc_content'], how='outer')
    shared_data.fillna(0, inplace=True)

    for col in ['is_lnc', 'is_mito', 'is_ribo', 'is_pc', 'is_oxphos', 'is_pseudo']:
        shared_data[col] = shared_data[col + "_x"] | shared_data[col + "_y"]
        shared_data.drop(col + "_x", axis=1, inplace=True)
        shared_data.drop(col + "_y", axis=1, inplace=True)

    model = sm.OLS(shared_data[f"{comparison_axis}_y"], shared_data[f"{comparison_axis}_x"]).fit()
    np.set_printoptions(suppress=True)
    shared_data['cooks_distance'] = model.get_influence().cooks_distance[0] + 1e-10

    xy = np.vstack([shared_data[f"{comparison_axis}_x"].to_numpy().flatten(), shared_data[f"{comparison_axis}_y"].to_numpy().flatten()])
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
