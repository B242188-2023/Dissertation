#!/bin/bash
samples=("LD1" "LD2" "LL1" "LL2" "SD1" "SD2" "SL1" "SL2")
adapter_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

for sample in "${samples[@]}"
do
    echo "Processing $sample..."
    fastp -i ${sample}_R1_001.fastq -I ${sample}_R2_001.fastq \
          -o ${sample}_R1_processed.fastq -O ${sample}_R2_processed.fastq \
          --adapter_sequence $adapter_R1 --adapter_sequence_r2 $adapter_R2 \
          --trim_front1 10 --trim_front2 10 \
          --length_required 30 \
          --dedup \
          -j ${sample}_fastp_report.json -h ${sample}_fastp_report.html
done