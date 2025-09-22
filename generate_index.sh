#! /bin/bash

outdir="/mnt/data1/10XvParse/Indexes/"
mkdir -p $outdir

if [[ $1 == m || $1 == mouse ]]; then
    echo "Building Mouse Index"

    echo "Retrieving Mouse Genome"
    wget -nc https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz -O "${outdir}"mmus.fa.gz

    echo "Retrieving Mouse GTF"
    wget -nc https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz -O "${outdir}"mmus.gtf.gz

    echo "Unzipping GTF"
    gunzip -k ${outdir}mmus.gtf.gz

    echo "Generating BED file from GTF"
    gffread "${outdir}"mmus.gtf -o "${outdir}"mmus.bed --bed

    echo "Building Kallisto Index for Mouse Genome"
    kb ref -i ${outdir}mouse.idx -g ${outdir}mouse_t2g.txt -f1 ${outdir}mouse_f1.fasta ${outdir}mmus.fa.gz ${outdir}mmus.gtf.gz -t 16

elif [[ $1 == h || $1 == human ]]; then
    echo "Building Human Index"

    echo "Retrieving Human Genome"
    wget -nc https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O "${outdir}"hsapiens.fa.gz

    echo "Retrieving Human GTF"
    wget -nc https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz -O "${outdir}"hsapiens.gtf.gz

    echo "Unzipping GTF"
    gunzip -k ${outdir}hsapiens.gtf.gz

    echo "Generating BED file from GTF"
    gffread "${outdir}"hsapiens.gtf -o "${outdir}"hsapiens.bed --bed

    echo "Building Kallisto Index for Human Genome"
    kb ref -i ${outdir}human.idx -g ${outdir}human_t2g.txt -f1 ${outdir}human_f1.fasta ${outdir}hsapiens.fa.gz ${outdir}hsapiens.gtf.gz -t 16

else 
    echo "Usage: $0 [mouse|m|human|h]"
    exit 1
fi