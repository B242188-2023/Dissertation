#!/bin/bash
adapter_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

samples=("LD1" "LD2" "LL1" "LL2" "SD1" "SD2" "SL1" "SL2")
for sample in "${samples[@]}"
do
  echo "Processing $sample..."
  cutadapt -a $adapter_R1 -o ${sample}_R1_trimmed.fastq ${sample}_R1_001.fastq
  cutadapt -a $adapter_R2 -o ${sample}_R2_trimmed.fastq ${sample}_R2_001.fastq
done
echo "Adapter trimming completed for all samples."