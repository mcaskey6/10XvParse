#!/bin/bash

dir="/home/mcaskey/10XvParse"

A2_paths=(
    "$dir/Data/Analysis_2/10x/FASTA/Processed/10x_0.fastq.gz"
    "$dir/Data/Analysis_2/parse/FASTA/Processed/parse_0.fastq.gz"
    "$dir/Data/Analysis_2/parse/FASTA/Processed/polyT_0.fastq.gz"
    "$dir/Data/Analysis_2/parse/FASTA/Processed/randO_0.fastq.gz"
    "$dir/Data/Analysis_2/parse_mini/FASTA/Processed/parse_mini_0.fastq.gz"
    "$dir/Data/Analysis_2/parse_mini/FASTA/Processed/polyT_0.fastq.gz"
    "$dir/Data/Analysis_2/parse_mini/FASTA/Processed/randO_0.fastq.gz"
)

A3_paths=(
    "$dir/Data/Analysis_3/10x_H1/FASTA/Processed/10x_H1_0.fastq.gz"
    "$dir/Data/Analysis_3/parse_H1/FASTA/Processed/parse_H1_0.fastq.gz"
    "$dir/Data/Analysis_3/parse_H1/FASTA/Processed/polyT_0.fastq.gz"
    "$dir/Data/Analysis_3/parse_H1/FASTA/Processed/randO_0.fastq.gz"
    "$dir/Data/Analysis_3/10x_H2/FASTA/Processed/10x_H2_0.fastq.gz"
    "$dir/Data/Analysis_3/parse_H2/FASTA/Processed/parse_H2_0.fastq.gz"
    "$dir/Data/Analysis_3/parse_H2/FASTA/Processed/polyT_0.fastq.gz"
    "$dir/Data/Analysis_3/parse_H2/FASTA/Processed/randO_0.fastq.gz"
)

for path in "${A2_paths[@]}"; do
    if [ ! -f "$path" ]; then
        echo "File not found: $path"
        exit 1
    fi
done

echo "Analysis 2 num reads:"
for path in "${A2_paths[@]}"; do
    num=$(( $(zcat "$path" | wc -l) / 4 ))
    echo "$path: $num"
done

for path in "${A3_paths[@]}"; do
    if [ ! -f "$path" ]; then
        echo "File not found: $path"
        exit 1
    fi
done

echo "Analysis 3 num reads:"
for path in "${A3_paths[@]}"; do
    num=$(( $(zcat "$path" | wc -l) / 4 ))
    echo "$path: $num"
done