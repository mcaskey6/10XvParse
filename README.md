# A Comparison of Parse Biosciences Evercode WT to 10x Genomics Chromium Single Cell 3' 

This repository reproduces and compares published single-cell RNA-seq datasets generated with Parse Biosciences Evercode and 10x Genomics Chromium Single Cell 3', using a unified preprocessing pipeline based on `kb-python` (kallisto|bustools).

## Repository Structure

```
10XvParse/
├── Configs/          # Per-analysis YAML config files (SRA accessions, reference URLs, technology strings)
├── Data/             # Downloaded FASTQs, pseudoalignment outputs, and H5AD files (gitignored)
├── Index/            # kallisto indices, organized by species
├── Logs/             # Log files from preprocessing runs
├── Notebooks/        # Jupyter notebooks for downstream analysis and plotting
├── Scripts/          # Python entry-point scripts, one per analysis
└── XvP_utils/        # Shared utility package
    ├── preprocessing/ # Download, subsampling, pseudoalignment pipelines
    └── plotting/      # Plotting helpers for comparison figures
```

## Environment Setup

Dependencies are managed with conda. To create and activate the environment:

```bash
conda env create -f environment.yml
conda activate 10XvParse
```

Install the local `XvP_utils` package in editable mode so that the scripts and notebooks can import it:

```bash
pip install -e .
```

## Running an Analysis

Each analysis has a dedicated script under `Scripts/` and a config file under `Configs/`. The scripts download raw reads from SRA/ERA, build or reuse a kallisto index, pseudoalign, and write `.h5ad` files to `Data/`.

```bash
python Scripts/analysis2.py   # Thymocyte dataset (mouse)
python Scripts/analysis3.py   # PBMC dataset (human)
python Scripts/analysis4.py   # K562/mESC barnyard dataset
python Scripts/analysis5.py   # Frozen PBMC dataset (human)
```

Downstream exploration is done in the corresponding Jupyter notebooks under `Notebooks/`.

## Datasets

### Analysis 2
From: [Comparative transcriptomic analyses of thymocytes using 10x Genomics and Parse scRNA-seq technologies - BMC Genomics](https://link.springer.com/article/10.1186/s12864-024-10976-x)

Sample: *Mus musculus* thymocytes from two female mice aged 6 months

Datasets (find [here](https://www.ncbi.nlm.nih.gov/sra?term=SRP484103)):
- GSM8020231-9: Parse Evercode WT v2, sublibrary 1–9
- GSM8020240 and GSM8020241: Parse Evercode Mini v2, sublibrary 1 and 2
- GSM8020242: 10x Genomics Next Gem v3, Gene Expression
- GSM8020243: 10x Genomics Next Gem v3, TotalSeq™-B hashtag antibodies (BioLegend)

### Analysis 3
From: [Comparative Analysis of Single-Cell RNA Sequencing Methods with and without Sample Multiplexing](https://www.mdpi.com/1422-0067/25/7/3828)

Sample: *Homo sapiens* PBMCs from two healthy individuals

Datasets (find [here](https://www.ncbi.nlm.nih.gov/sra?term=SRP469371)):
- GSM7873659,61,63,65,67,69,71,73: Parse Evercode WT v2 H1 replicates 1–8
- GSM7873660,62,64,66,68,70,72,74: Parse Evercode WT v2 H2 replicates 1–8
- GSM7873657: 10x Genomics Next Gem v3 H1
- GSM7873658: 10x Genomics Next Gem v3 H2

### Analysis 4
From: [Comparison of Single Cell Transcriptome Sequencing Methods: Of Mice and Men](https://www.mdpi.com/2073-4425/14/12/2226)

Sample: Mixture of K562 human multiple myeloma and mESC mouse embryonic stem cell lines (barnyard)

Datasets uploaded to EBI under PRJEB67544 (10x) and PRJEB67549 (Parse):
- ERR12398015: 10x Next Gem v3
- ERR12167397 and ERR12167398: Parse Evercode Mini v2 sublibrary 1 and 2

### Analysis 5
From: [Comparative Analysis of Commercial Single-Cell RNA Sequencing Technologies](https://www.biorxiv.org/content/10.1101/2024.06.18.599579v1.full#sec-15)

Sample: Frozen *Homo sapiens* PBMCs

Datasets (find [here](https://www.ncbi.nlm.nih.gov/sra?term=SRP484103)):
- SRR28867558: Parse Evercode WT v2
- SRR28867563 and SRR28867562: 10x Next Gem v3, technical replicates 1 and 2

## Adding a New Dataset

Follow these steps to add a new analysis:

### 1. Create a config file

Copy an existing config from `Configs/` and edit it for your dataset. The config has four required sections:

```yaml
name: Analysis_N

# SRA or ERA accession numbers, grouped by assay type
SRA:
  10x:
    - SRR_XXXXXXX
  parse:
    - SRR_XXXXXXX

# R1/R2 file suffixes produced by fasterq-dump (these vary across datasets — check with getReadNum.out)
read_num:
  10x:
    R1: 1
    R2: 2
  parse:
    R1: 1
    R2: 2

# Ensembl reference files for building the kallisto index
reference:
  species: human   # or "mouse"; controls index reuse from Index/
  fasta: https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  gtf: https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

# kb-python technology strings
tech:
  10x: 10XV3
  parse: "1,10,18,1,48,56,1,78,86:1,0,10:0,0,0"
```

Use `ERA` instead of `SRA` for European Nucleotide Archive accessions. For barnyard (dual-species) experiments, provide both a human and mouse `fasta`/`gtf` — see `analysis4.yaml` for the pattern.

### 2. Add Parse barcode config files

For each Parse assay in the config, create a subdirectory `Configs/Analysis_N/<assay_name>/` and place the following five static files there. These files are specific to the Parse kit version (e.g. v2 vs. v3) and the demultiplexing strategy (e.g all vs. subset of wells used for sample) — copy them from an existing analysis that used the same kit version:

| File | Contents |
|------|----------|
| `r1_R.txt` | Round 1 Random primer barcode sequences used by splitcode to demultiplex by well |
| `r1_T.txt` | Round 1 PolyT primer barcode sequences used by splitcode to demultiplex by well |
| `onlist.txt` | Barcode whitelist with allowed sequences and error-correction targets |
| `replace.txt` | Barcode replacement table mapping raw barcodes to corrected sequences |
| `bcs_to_wells.txt` | Mapping from barcode sequences to well positions |

The pipeline also writes generated files (`config_RT_parse.txt`, `parse_keep.txt`, `randOpolyT_keep.txt`) into the same subdirectory at run time — these do not need to be created manually.

### 3. Write a script

Create `Scripts/analysisN.py` based on an existing script. At minimum, call `load_10x` and `load_parse` (and `subsample_*` variants if you want depth-matched comparisons):

```python
from XvP_utils.preprocessing import load_10x, load_parse, subsample_parse, subsample_10x, get_subsample_num, setup_logger
from XvP_utils.preprocessing import RunSettings
from pathlib import Path
import os

if __name__ == "__main__":
    settings = RunSettings(
        root_dir=Path(__file__).parent.parent,
        config_name="analysisN.yaml",
        overwrite=False,
        threads=16,
        max_workers=4,
    )

    config_file = settings.root_dir / "Configs" / settings.config_name
    os.makedirs(settings.root_dir / "Logs", exist_ok=True)
    logger = setup_logger(settings.root_dir / "Logs" / "analysisN.txt")

    load_10x(settings, config_file, "10x", logger)
    load_parse(settings, config_file, "parse", logger)

    subsample_num = get_subsample_num(
        settings, config_file,
        ten_x_assays=["10x"],
        parse_assays=["parse"],
        logger=logger,
    )
    subsample_10x(settings, config_file, "10x", subsample_num, logger)
    subsample_parse(settings, config_file, "parse", subsample_num, logger)
```

### 4. Run the script

```bash
python Scripts/analysisN.py
```

kb python outputs with `.h5ad` files will be written to `Data/Analysis_N/{assat}/kb_python`. Logs go to `Logs/analysisN.txt`.

### 5. Add a notebook

Create a directory `Notebooks/Analysis_N/` and add a Jupyter notebook for downstream analysis and figures, following the pattern in `Notebooks/Analysis_2/` or `Notebooks/Analysis_3/`.

### 6. Update this README

Add an entry for the new dataset under the **Datasets** section above, including the paper link, sample description, and accession numbers.
