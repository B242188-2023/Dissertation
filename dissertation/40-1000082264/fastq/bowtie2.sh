#!/bin/bash
#bowtie2-build AeUmbellulata_TA1851_v1-cds.fasta Index-cds
#bowtie2-build AeUmbellulata_TA1851_v1.fasta Index

samples=("LD1" "LD2" "LL1" "LL2" "SD1" "SD2" "SL1" "SL2")
#samples=("LD1")

for sample in "${samples[@]}"
do
    echo "Processing $sample..."
    bowtie2 -x Index-cds -1 ${sample}_R1_processed.fastq -2 ${sample}_R2_processed.fastq \
    -p 32 --very-sensitive-local | \
    samtools sort -@ 160 -o ../alignment/${sample}_cds_sorted.bam
    echo "Completed aligning $sample."
done