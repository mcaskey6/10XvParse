# A Comparison of Parse Biosciences Evercode WT to 10x Genomics Chromium Single Cell 3'

This repository reproduces and compares published single-cell RNA-seq datasets generated with Parse Biosciences Evercode and 10x Genomics Chromium Single Cell 3', using a unified [Snakemake](https://snakemake.github.io/) preprocessing workflow built on `kb-python` (kallisto|bustools). From raw reads it produces gene-count matrices (`.h5ad`), subsampling the two technologies to a matched sequencing depth so they can be compared fairly, plus optional STARsolo + RSeQC gene-body coverage plots.

## Repository Structure

```
10XvParse/
├── Config/                  # Snakemake configuration
│   ├── config.yaml          # thread count + the analyses to build (name -> config file)
│   ├── indexes.yaml         # per-species reference URLs (+ HRT Atlas URL; barnyard keys)
│   └── analysisN.yaml       # one file per analysis (assays, read sources, comparisons)
├── Envs/                    # conda environment specs you create by hand
│   ├── environment.yml      # the 10XvParse env (Snakemake + tools + XvP_utils via pip)
│   └── goseq.yaml           # small separate R env for the goseq GO-enrichment step
├── workflow/                # the Snakemake workflow (self-contained; lowercase by Snakemake convention)
│   ├── Snakefile            # all rules
│   ├── paths.py             # every input/output path (the naming conventions)
│   ├── config_helpers.py    # reading + interpreting the per-analysis configs
│   ├── scripts/             # the few Python steps (Parse config gen, batch file, STAR params)
│   ├── envs/                # per-rule conda envs for `--use-conda` (distinct from Envs/)
│   └── profiles/default/    # default run settings (cores)
├── Resources/               # committed static inputs
│   ├── parse_info/          # Parse kit barcode reference (kits_info.txt + barcodes/*.csv)
│   └── Analysis_N/<assay>/  # per-assay committed inputs (hashtags.tsv, hto_demultiplexed.csv)
├── Generated/               # workflow-generated files: Parse/STAR configs, 10x whitelists, read_counts (gitignored)
├── Data/                    # FASTQs, count matrices (.h5ad), BAMs, and plots (gitignored)
├── Index/                   # kallisto + STAR indices and genome references, by species
├── Tests/                   # synthetic end-to-end test + CI environment
├── Notebooks/               # Jupyter notebooks for downstream analysis and figures
└── XvP_utils/               # plotting/analysis helpers imported by the notebooks
```

> The workflow separates **committed static inputs** (`Config/`, `Resources/`) from files it **generates** (`Generated/`, `Data/`). Everything under `Generated/` and `Data/` is git-ignored and rebuilt on demand, so only the true inputs are version-controlled. `XvP_utils` holds the plotting/analysis helpers the notebooks import; it is not used by the workflow, which depends on nothing outside `workflow/`.

### Parse barcode config files are auto-generated

For each Parse assay the workflow generates the following files in `Generated/Analysis_N/<assay>/` at run time, derived from the kit name and optional well list in the analysis config:

| File | Contents |
|------|----------|
| `r1_R.txt` | Round 1 randO primer barcode sequences for the selected wells (used by splitcode) |
| `r1_T.txt` | Round 1 polyT primer barcode sequences for the selected wells (used by splitcode) |
| `onlist.txt` | Per-round barcode whitelist for kb-python error correction, including the sublibrary barcode |
| `replace.txt` | Maps randO bc1 sequences to their polyT counterparts for kb-python |
| `bcs_to_wells.txt` | Mapping from all bc1 sequences to well positions (used in analysis notebooks) |
| `lib_bc.txt`, `star_bc1/2/3.txt` | Per-round barcode whitelists for STARsolo (gene-body coverage) |
| `sublibraries.txt` | Which sublibrary barcode was assigned to which sublibrary name |
| `x_string.txt` | The kb-python technology x-string (shifted for the sublibrary barcode) |

The source of truth for all barcode sequences is the committed `Resources/parse_info/`: `kits_info.txt` lists the barcode files for each kit version, and `barcodes/*.csv` holds the sequences. Alongside the files above, the workflow writes a few small splitcode/STAR helper files into `Generated/Analysis_N/<assay>/` (`config_RT_parse.txt`, `parse_keep.txt`, `randOpolyT_keep.txt`, and the `star_*` position/whitelist files).

## Environment Setup

Dependencies are managed with conda:

```bash
conda env create -f Envs/environment.yml
conda activate 10XvParse
```

This provides Snakemake and the command-line tools the workflow calls — `kb-python`, `splitcode`, `sra-tools`, `seqtk`, `pigz`, `STAR`, `samtools`, `gffread`, and RSeQC's `geneBody_coverage.py` — and (via `pip install -e .`) the local `XvP_utils` package the `Notebooks/` import, along with its Python dependencies from `pyproject.toml`. The `workflow/` package itself imports nothing beyond the standard library.

The GO-enrichment step (`XvP_utils.cross_comparison`) runs **goseq** through R, which lives in a small separate environment — kept apart so `XvP_utils`' resolver doesn't pick this env's goseq-less `Rscript`:

```bash
conda env create -f Envs/goseq.yaml
# point XvP at it (or let the sibling-env search find it):
export XVP_RSCRIPT=$CONDA_PREFIX/../goseq/bin/Rscript
```

## Running the Pipeline

Run from the repository root, in the `10XvParse` environment:

```bash
snakemake --workflow-profile workflow/profiles/default all
```

`all` builds the gene-count matrices for every analysis registered in `Config/config.yaml`. To build only part of the workflow, name specific targets:

```bash
# one assay's count matrix
snakemake --workflow-profile workflow/profiles/default \
  Data/Analysis_5/10x/kb_python/10x_out/counts_unfiltered/adata.h5ad

# a gene-body coverage plot for a comparison group (empty tag -> ".geneBodyCoverage")
snakemake --workflow-profile workflow/profiles/default \
  Data/Analysis_5/Plots/.geneBodyCoverage.txt
```

Add `-n` for a dry run (prints the jobs without running them). For a given target the workflow will, as needed: fetch the reads (SRA via `prefetch`/`fasterq-dump`, ENA via FTP, or pre-downloaded local files), build or reuse the kallisto and STAR indices, remultiplex libraries with splitcode, filter/split Parse reads, subsample each comparison group to its shared minimum read count, pseudoalign with `kb count`, and write `.h5ad` matrices under `Data/<analysis>/<assay>/kb_python/`. Gene-body plots align the subsampled reads with STARsolo and run RSeQC.

### Per-rule conda environments

Every rule declares a conda environment (`workflow/envs/*.yaml`, grouped by tool — `kb`, `sra`, `splitcode`, `seqtk`, `star`, `coverage`, and a small `base` for the download/unzip/Python-script rules). Pass `--use-conda` (or `--software-deployment-method conda`) to run each rule in its own pinned environment, which Snakemake creates on first use:

```bash
snakemake --workflow-profile workflow/profiles/default --use-conda all
```

Without `--use-conda`, rules run in whatever environment is active (e.g. the `10XvParse` env from `Envs/environment.yml`), which must then provide the tools itself.

### Logs

Each job redirects its stderr to a per-job log under `Logs/` (git-ignored): analysis work under `Logs/<analysis>/<rule>/…` and shared reference/index builds under `Logs/reference/<rule>/…`. Snakemake's own run log stays under `.snakemake/log/`; STAR additionally writes its `Log.*` files next to each BAM.

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

Sample: Frozen *Homo sapiens* PBMCs from two healthy individuals

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
Now: [A comprehensive analysis framework for evaluating commercial single-cell RNA sequencing technologies] https://academic.oup.com/nar/article/53/2/gkae1186/7924191#501285411

Sample: Frozen *Homo sapiens* PBMCs

Datasets (find [here](https://www.ncbi.nlm.nih.gov/sra?term=SRP505235)):
- SRR28867558: Parse Evercode WT v2
- SRR28867563 and SRR28867562: 10x Next Gem v3.1, technical replicates 1 and 2

### Analysis 6
From: [Comparison of high-throughput single-cell RNA-seq methods for ex vivo drug screening](https://academic.oup.com/nargab/article/6/1/lqae001/7591100?login=true#460158720)

Sample: *Homo sapiens* glucocorticoid-resistant E/R+ ALL Reh cell line multiplexed according to drug treatment (6 total perturbation experiments)

Datasets (find [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229617))
- SRR24154339: 10X Next Gem v3
- SRR24154340: Multi-seq barcodes
- SRR24154341-2: Parse Evercode mini v1 replicate 1 and 2

### Analysis 7
From: Comparative analysis of multiplex single-cell mRNA sequencing of resting and activated 
PBMCs using droplet-based and split-pool methods (Yet to be published)

Sample: Frozen *Homo sapiens* PBMCs

Datasets:
10x: AZ13332/AZ_cDNA_S1_L002; 10X Chromium v4
10x_hashtags: LMO/RPI7_S0_L001; 10X Multi-Seq
parse: AZ12601/AZ_PS_5k_S5_L002 and AZ12601/AZ_PS_10k_S6_L002; Parse Evercode mini v3


## Paper Figure Labels

The comparison figures in `Notebooks/Comparisons/` use the following experiment labels, ordered by species/tissue. Labels are colored by tissue type in the figures.

| Paper Label | Repository | Sample | Species |
|---|---|---|---|
| Exp 1 | Analysis_6 | Cancer cell perturbation (Reh) | Human |
| Exp 2a | Analysis_3 (H1) | Frozen PBMC, Donor 1 | Human |
| Exp 2b | Analysis_3 (H2) | Frozen PBMC, Donor 2 | Human |
| Exp 3 | Analysis_5 | Frozen PBMC | Human |
| Exp 4 | Analysis_7 | Frozen PBMC | Human |
| Exp 5a | Analysis_2 (standard) | Thymocytes | Mouse |
| Exp 5b | Analysis_2 (mini) | Thymocytes | Mouse |
| Exp 6 | Analysis_4 | K562/mESC barnyard | Human/Mouse |

## Adding a New Dataset

Adding an analysis is entirely config-driven — no workflow code changes are needed.

### 1. Add the reference (only if it's a new species)

Reference URLs live under the `indexes:` key of `Config/indexes.yaml`. `human`, `mouse`, and the `human_mouse` barnyard reference already exist. To add a species or a different Ensembl release, add an entry:

```yaml
indexes:
  human:
    fasta: https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    gtf: https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
    hrt_atlas_url: https://raw.githubusercontent.com/Bidossessih/HRT_Atlas/master/www/Housekeeping_GenesHuman.csv
```

The index directory (`Index/{species}/`) is shared across all analyses of the same species. A barnyard (dual-species) reference for species `"{sp1}_{sp2}"` sets `barnyard: true` and adds `fasta_2`/`gtf_2`/`hrt_atlas_url_2` for the second species — see the `human_mouse` entry.

### 2. Write the analysis config

Create `Config/analysisN.yaml`. Each assay is named as a key under `tech:`; its reads come from an `SRA:`, `ERA:`, or `local:` block, and `comparisons:` says which 10x/Parse assays are subsampled together.

```yaml
name: Analysis_N

species: human            # or "mouse" / "human_mouse"; picks the Config/indexes.yaml entry

# Read source. Use one of SRA / ERA / local per assay. Parse assays group their
# accessions by sublibrary (see "Parse sublibraries" below); 10x assays are flat.
SRA:
  10x:
    - SRR_XXXXXXX
  parse:
    sub1:
      - SRR_XXXXXXX
      - SRR_XXXXXXX   # a second run of the same sublibrary
    sub2:
      - SRR_XXXXXXX

# R1/R2 file suffixes (SRA: the fasterq-dump split index; ERA: the ENA file suffix).
# This varies across datasets. Add a per-sublibrary override by nesting an {R1, R2}
# dict under the sublibrary name.
read_num:
  10x:   {R1: 1, R2: 2}
  parse: {R1: 1, R2: 2}

# kb-python technology string for 10x; kit name (<kit>_v<chem>) for Parse.
tech:
  10x: 10XV3
  parse: WT_v2

  # Optional: restrict Parse to a subset of wells (all wells if omitted).
  wells:
    parse: [A1, A2, A3]

# Optional: trim a feature-barcode (hashtag) assay's R2 to N bases.
# trim:
#   10x_hashtags: 8

# Subsample comparison groups. Each group lists the 10x and Parse assays that are
# downsampled together to their shared minimum read count. Give each group a `tag`
# when the same assay is subsampled at more than one depth (e.g. one 10x compared
# against both a standard and a mini Parse kit); use "" for a single depth. Each
# group also defines a gene-body coverage comparison (its one 10x vs its Parse).
comparisons:
  - tag: ""
    tenx: [10x]
    parse: [parse]
```

For Parse assays the `tech` value must be a kit name of the form `<kit>_v<chem>` (e.g. `WT_v2`, `WT_mini_v3`); it is looked up in `Resources/parse_info/kits_info.txt`. Feature-barcode ("hashtag") assays are recognised by `hashtag` in the assay name and use a kite index built from a committed `Resources/Analysis_N/<assay>/hashtags.tsv`.

**Pre-downloaded (local) reads.** For files you already have, use a `local:` block instead of `SRA`/`ERA`. Each library is an `[R1, R2]` filename pair (in the assay's read order) living in `Data/Analysis_N/<assay>/FASTA/Dumped/`; filenames are used as given, so arbitrary Illumina names work without renaming, and each pair is one Parse sublibrary in the order listed:

```yaml
local:
  10x:
    - [AZ_cDNA_S1_L002_R1_001.fastq.gz, AZ_cDNA_S1_L002_R2_001.fastq.gz]
  parse:
    - [AZ_PS_5k_S5_L002_R1_001.fastq.gz, AZ_PS_5k_S5_L002_R2_001.fastq.gz]
    - [AZ_PS_10k_S6_L002_R1_001.fastq.gz, AZ_PS_10k_S6_L002_R2_001.fastq.gz]
```

### 3. Register the analysis

Add it under `analyses:` in `Config/config.yaml`:

```yaml
analyses:
  Analysis_N: Config/analysisN.yaml
```

### 4. Run it

```bash
snakemake --workflow-profile workflow/profiles/default \
  Data/Analysis_N/10x/kb_python/10x_out/counts_unfiltered/adata.h5ad
# ...or just `all` to (re)build everything.
```

### 5. Add a notebook and update this README

Add a `Notebooks/Analysis_N/` notebook for downstream analysis, and add an entry under **Datasets** above with the paper link, sample description, and accession numbers.

### Parse sublibraries

A Parse sublibrary is distinguished by its sequencing index, not by the combinatorial barcodes themselves, so the same bc1/bc2/bc3 combination in two sublibraries belongs to two different cells. Sublibraries therefore cannot simply be concatenated. Group each Parse assay's accessions under a sublibrary name; `splitcode --remultiplex --bclen=4` then prepends a distinct 4 bp sequence to read 2 of each sublibrary, and that sequence is carried downstream as a fourth cell barcode. Accessions listed under the same name are runs of one sublibrary and share a barcode; for a `local:` block, each `[R1, R2]` pair is its own sublibrary.

Two consequences worth knowing:

- Barcodes are assigned in the order sublibraries first appear in the config, so **reordering them changes the barcode assignment of already-processed data**. The assignment for a given run is recorded in `Generated/Analysis_N/<assay>/sublibraries.txt`.
- Parse cell barcodes in the resulting H5ADs are consequently 4 bp longer, with the sublibrary barcode leading. bc1 remains the trailing 8 bp, so notebook logic keyed on the end of the barcode is unaffected.
