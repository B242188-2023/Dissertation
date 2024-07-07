#!/bin/bash
#snpEff build -gff3 -v AeUmbellulata -c snpEff.config > output.log 2>&1
#python3 id-transform.py 

samples=("LD1" "LD2" "LL1" "LL2" "SD1" "SD2" "SL1" "SL2")
#samples=("LD1")
for sample in "${samples[@]}"
do
    echo "Processing $sample..."
    #bcftools annotate --rename-chrs rename.txt ${sample}_snp.vcf -o ${sample}_snp_renamed.vcf

    snpEff ann -htmlStats ../annotation/${sample}_snpEff_report.html \
     AeUmbellulata ${sample}_output.vcf -c snpEff.config > \
     ../annotation/${sample}_annotated_variants.vcf
    echo "Completed processing $sample."
done