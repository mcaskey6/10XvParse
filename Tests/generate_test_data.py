"""Generate a tiny synthetic dataset for the 10XvParse end-to-end test.

Everything here is committable *code* rather than committed binary blobs: the
reference, the 10x reads, and the Parse reads are all synthesized deterministically
(fixed RNG seed) at test time, so the repo stays free of large FASTQ/FASTA files
(which .gitignore excludes anyway) and no network download is needed.

It writes, under the repo root (all git-ignored, staged fresh each run):

  Index/testsp/ref.fa.gz, ref.gtf.gz
      A ~3 kb single-chromosome reference with two 2-exon genes. Writing it as the
      *download_reference outputs* means Snakemake skips the download and `kb ref`
      builds a tiny index from these directly.

  Data/Analysis_test/10x/FASTA/Dumped/Lib0_{1,2}.fastq.gz
      Synthetic 10XV3 reads: R1 = 16 bp cell barcode (drawn from real 10XV3
      whitelist entries, so kb's onlist correction keeps them) + 12 bp UMI;
      R2 = a 90 bp window of a transcript.

  Data/Analysis_test/parse/FASTA/Dumped/Lib{0,1}_{1,2}.fastq.gz
      Synthetic WT_mini_v3 Parse reads for two sublibraries. R1 = cDNA; R2 = the
      barcode read laid out for the kit's x_string (UMI + bc3 + bc2 + bc1, with the
      linker gaps the kit expects), bc1/bc2/bc3 drawn from the committed
      Resources/parse_info barcode tables. Both polyT (T-type bc1) and randO (R-type
      bc1) reads are emitted so the barcode-split branch is exercised.

The read layout deliberately matches what the workflow expects *before*
remultiplexing: splitcode later prepends the 4 bp sublibrary barcode to the Parse
barcode read, which is why the raw barcode read starts at the un-shifted positions.

Run: python Tests/generate_test_data.py
"""
from __future__ import annotations

import csv
import gzip
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PARSE_INFO = ROOT / "Resources" / "parse_info"

RNG = random.Random(1234)

# 16 real 10XV3 whitelist barcodes, sampled far apart across the whitelist so kb's
# 1-mismatch correction can't collapse two of them into one cell.
TENX_V3_BARCODES = [
    "AAACCCAAGAAACACT", "ACATTTCGTCCAGGTC", "AGCGATTAGGCACTCC", "ATCGATGTCAAATGGG",
    "CACGGGTCAAGGCTTT", "CCCGGAATCGCCAATA", "CGTAGTACATCCGCAG", "CTGTGGGAGAGCGCGT",
    "GATCAGTGTGAACGGT", "GGATCTAAGGTGCTCT", "GTCATGATCAGGCGAA", "TACACCCCACCCAAGC",
]


def randseq(n: int) -> str:
    return "".join(RNG.choice("ACGT") for _ in range(n))


def write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    """Write (name, seq) records as a gzipped FASTQ with all-'I' (Q40) qualities."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


# --------------------------------------------------------------------------- #
# reference: one chromosome, two genes, two exons each
# --------------------------------------------------------------------------- #

# (gene_id, [(exon_start, exon_end), ...]) — 1-based inclusive, all on the + strand.
GENES = {
    "geneA": [(101, 300), (501, 700)],
    "geneB": [(1101, 1300), (1601, 1800)],
}
CHROM = "testchr1"
CHROM_LEN = 3000


def build_reference() -> dict[str, str]:
    """Write ref.fa.gz + ref.gtf.gz and return {gene_id: spliced_cDNA}."""
    genome = randseq(CHROM_LEN)
    index_dir = ROOT / "Index" / "testsp"
    index_dir.mkdir(parents=True, exist_ok=True)

    with gzip.open(index_dir / "ref.fa.gz", "wt") as fa:
        fa.write(f">{CHROM}\n")
        for i in range(0, len(genome), 60):
            fa.write(genome[i:i + 60] + "\n")

    cdnas: dict[str, str] = {}
    with gzip.open(index_dir / "ref.gtf.gz", "wt") as gtf:
        for gene_id, exons in GENES.items():
            tx_id = gene_id.replace("gene", "tx")
            gene_start, gene_end = exons[0][0], exons[-1][1]
            attr_g = f'gene_id "{gene_id}"; gene_name "{gene_id}";'
            attr_t = f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_name "{gene_id}";'
            gtf.write(f"{CHROM}\tsynthetic\tgene\t{gene_start}\t{gene_end}\t.\t+\t.\t{attr_g}\n")
            gtf.write(f"{CHROM}\tsynthetic\ttranscript\t{gene_start}\t{gene_end}\t.\t+\t.\t{attr_t}\n")
            cdna = ""
            for start, end in exons:
                gtf.write(f"{CHROM}\tsynthetic\texon\t{start}\t{end}\t.\t+\t.\t{attr_t}\n")
                cdna += genome[start - 1:end]  # 1-based inclusive -> 0-based slice
            cdnas[gene_id] = cdna
    return cdnas


def cdna_window(cdna: str, length: int = 90) -> str:
    start = RNG.randint(0, len(cdna) - length)
    return cdna[start:start + length]


# --------------------------------------------------------------------------- #
# 10x reads
# --------------------------------------------------------------------------- #

def generate_tenx(cdnas: dict[str, str]) -> None:
    r1_records, r2_records = [], []
    n = 0
    for cb in TENX_V3_BARCODES:
        for gene_id, cdna in cdnas.items():
            for _ in range(4):  # a few UMIs per cell/gene
                umi = randseq(12)
                name = f"tenx_read{n}"
                r1_records.append((name, cb + umi))          # 16 bp CB + 12 bp UMI
                r2_records.append((name, cdna_window(cdna)))  # cDNA
                n += 1
    dumped = ROOT / "Data" / "Analysis_test" / "10x" / "FASTA" / "Dumped"
    write_fastq(dumped / "Lib0_1.fastq.gz", r1_records)
    write_fastq(dumped / "Lib0_2.fastq.gz", r2_records)


# --------------------------------------------------------------------------- #
# Parse reads (WT_mini_v3: bc1=n38_R1_v3_8, bc2=v1, bc3=R3_v3)
# --------------------------------------------------------------------------- #

def _load_bc(name: str) -> list[dict[str, str]]:
    with open(PARSE_INFO / "barcodes" / f"bc_data_{name}.csv") as fh:
        return list(csv.DictReader(fh))


def _parse_barcode_picks():
    bc1_rows = _load_bc("n38_R1_v3_8")
    wells = ["A1", "A2", "A3"]
    bc1_T = {r["well"]: r["sequence"] for r in bc1_rows if r["stype"] == "T"}
    bc1_R = {r["well"]: r["sequence"] for r in bc1_rows if r["stype"] == "R"}
    # T-type bc1 -> polyT branch, R-type -> randO branch (paired by well so the
    # randO->polyT replace table merges them).
    bc1_picks = [(bc1_T[w], "T") for w in wells] + [(bc1_R[w], "R") for w in wells]
    bc2_picks = [r["sequence"] for r in _load_bc("v1") if r["stype"] == "L"][:2]
    bc3_picks = [r["sequence"] for r in _load_bc("R3_v3") if r["stype"] == "L"][:2]
    return bc1_picks, bc2_picks, bc3_picks


def _barcode_read(umi: str, bc3: str, bc2: str, bc1: str) -> str:
    # WT_mini_v3 x_string 1,10,18,1,30,38,1,50,58 : 1,0,10 (before the 4 bp
    # sublibrary barcode splitcode prepends): UMI[0:10] bc3[10:18] bc2[30:38]
    # bc1[50:58]; the [18:30] and [38:50] gaps are kit linker sequence.
    linker = "GATCGATCGATC"  # 12 bp filler, position-only (never read as a barcode)
    return umi + bc3 + linker + bc2 + linker + bc1  # 10+8+12+8+12+8 = 58 bp


def generate_parse(cdnas: dict[str, str]) -> None:
    bc1_picks, bc2_picks, bc3_picks = _parse_barcode_picks()
    for lib in (0, 1):
        r1_records, r2_records = [], []  # R1 = cDNA, R2 = barcode read
        n = 0
        for bc1, _bc1type in bc1_picks:
            for bc2 in bc2_picks:
                for bc3 in bc3_picks:
                    for gene_id, cdna in cdnas.items():
                        umi = randseq(10)
                        name = f"parse_lib{lib}_read{n}"
                        r1_records.append((name, cdna_window(cdna)))
                        r2_records.append((name, _barcode_read(umi, bc3, bc2, bc1)))
                        n += 1
        dumped = ROOT / "Data" / "Analysis_test" / "parse" / "FASTA" / "Dumped"
        write_fastq(dumped / f"Lib{lib}_1.fastq.gz", r1_records)
        write_fastq(dumped / f"Lib{lib}_2.fastq.gz", r2_records)


def main() -> None:
    cdnas = build_reference()
    generate_tenx(cdnas)
    generate_parse(cdnas)
    print("Wrote synthetic reference to Index/testsp/ and reads to Data/Analysis_test/")


if __name__ == "__main__":
    main()
