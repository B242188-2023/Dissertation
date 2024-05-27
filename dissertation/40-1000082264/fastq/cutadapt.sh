#!/bin/bash
adapter_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

samples=("LD1" "LD2" "LL1" "LL2" "SD1" "SD2" "SL1" "SL2")
for sample in "${samples[@]}"
do
  echo "Processing $sample..."
  # Step 1: Cutadapt Remove adapter
  cutadapt -a $adapter_R1 -o ${sample}_R1_trimmed.fastq ${sample}_R1_001.fastq
  cutadapt -a $adapter_R2 -o ${sample}_R2_trimmed.fastq ${sample}_R2_001.fastq

  # Step 2: Use seqtk filter empty sequence
  seqtk seq -L 1 ${sample}_R1_trimmed.fastq > ${sample}_R1_filtered.fastq
  seqtk seq -L 1 ${sample}_R2_trimmed.fastq > ${sample}_R2_filtered.fastq

  # Step 3: Use fastx_trimmer trims first 10 bases
  fastx_trimmer -f 11 -i ${sample}_R1_filtered.fastq -o ${sample}_R1_headtrimmed.fastq
  fastx_trimmer -f 11 -i ${sample}_R2_filtered.fastq -o ${sample}_R2_headtrimmed.fastq

  # Step 4: Use clumpify remove repeated sequences 
  clumpify.sh in=${sample}_R1_headtrimmed.fastq out=${sample}_R1_dedup.fastq dedupe
  clumpify.sh in=${sample}_R2_headtrimmed.fastq out=${sample}_R2_dedup.fastq dedupe
done
echo "Adapter trimming completed for all samples."