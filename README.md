# xQTL_genomic_architecture_analysis

# DGRP
## Generate Sample Files through:
````
#!/bin/bash

VCF="filtered-all.vcf.gz"

# Build header from sample names
samples=($(bcftools query -l $VCF))
header="CHROM\tPOS"
for s in "${samples[@]}"; do
    header+="\tREF_${s}\tALT_${s}"
done
echo -e "$header"

# Extract AD (allele depth) field: REF_count,ALT_count per sample
bcftools query -f '%CHROM\t%POS[\t%AD]\n' $VCF | awk '{
    printf "%s\t%s", $1, $2
    for(i=3; i<=NF; i++) {
        split($i, a, ",")
        if($i == "." || a[1] == ".") printf "\t0\t0"
        else printf "\t%s\t%s", a[1], a[2]
    }
    printf "\n"
}'

````

