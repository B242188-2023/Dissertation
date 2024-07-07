#!/bin/bash
samples=("LD1" "LD2" "LL1" "LL2" "SD1" "SD2" "SL1" "SL2")
#samples=("LD1")
for sample in "${samples[@]}"
do
    echo "Processing $sample..."

    # extracting unmatched reads 
    #samtools view -@ 160 -b ${sample}_cds_sorted.bam > ${sample}_cds_mapped_sorted.bam
    # bam -> fastq
    #samtools fastq ${sample}_cds_mapped_sorted.bam > ${sample}_cds_mapped.fastq
    # de novo assembly
    #spades.py -s ${sample}_cds_unmapped.fastq -o spades_output_${sample}
    # matching
    #bowtie2 -x ../fastq/Index-cds -f spades_output_${sample}/contigs.fasta\
    # | samtools view -bS - | samtools sort -o ${sample}_cds_denovo_sorted.bam
    
    
    # sort
    # samtools sort -@ 160 ${sample}_cds_unmapped.bam -o ${sample}_cds_unmapped_sorted.bam
    # statistical analyze
    #samtools flagstat ${sample}_cds_sorted.bam > ${sample}_cds_sorted.stats
    
    # index
    # samtools index -@ 160 ${sample}_cds_sorted.bam

    featureCounts -a AeUmbellulata_TA1851_v1.gtf -o ${sample}_counts.txt \
    ${sample}_genome_sorted.bam > ${sample}_featureCounts.log 2>&1

    # identify variations
    bcftools mpileup --threads 160 -Ou -f ../fastq/AeUmbellulata_TA1851_v1-cds.fasta\
     ${sample}_cds_sorted.bam | bcftools call -mv -Oz -o ${sample}_variants.vcf.gz

    # filtering low quality variants
    bcftools view ${sample}_variants.vcf.gz | bcftools filter --exclude 'QUAL < 20' \
    | bgzip -c > ${sample}_filtered_variants.vcf.gz

    bcftools view -v snps -m2 -M2 --min-ac 1:minor ${sample}_filtered_variants.vcf.gz > ${sample}_snp.vcf
    echo "Completed processing $sample."
done