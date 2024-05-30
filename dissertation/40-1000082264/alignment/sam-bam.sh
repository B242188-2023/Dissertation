#!/bin/bash
samples=("LD1" "LD2" "LL1" "LL2" "SD1" "SD2" "SL1" "SL2")
for sample in "${samples[@]}"
do
    echo "Processing $sample..."
    # format conversion sam->bam
    samtools view -b ${sample}_aligned.sam -o ${sample}_aligned.bam
    # sort
    samtools sort -o ${sample}_sorted.bam ${sample}_aligned.bam
    
    echo "Completed processing $sample."
done