#!/bin/bash
snpEff build -gff3 -v AeUmbellulata -c snpEff.config > output.log 2>&1
python3 id-transform.py 

samples=("LD1" "LD2" "LL1" "LL2" "SD1" "SD2" "SL1" "SL2")
#samples=("LD1")
for sample in "${samples[@]}"
do
    echo "Processing $sample..."
    bcftools annotate --rename-chrs rename.txt ${sample}_variants.vcf\
     | snpEff ann AeUmbellulata - -c snpEff.config > \
     ../annotation/${sample}_annotated_variants.vcf
    echo "Completed processing $sample."
done
