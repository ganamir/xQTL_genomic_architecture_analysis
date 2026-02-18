# xQTL_genomic_architecture_analysis

# DGRP
## 1. Generate Sample Files through:
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

## 2. Trim files
````
#!/bin/bash
DIR="/mnt/d/xQTL_2025_Data/DGRP_XQTL/trimmed_fastqs"

for R1 in "$DIR"/*_merged_R1_001.fastq.gz; do
    BASE=$(basename "$R1" _merged_R1_001.fastq.gz)
    R1_FILE="$DIR/${BASE}_merged_R1_001.fastq.gz"
    R2_FILE="$DIR/${BASE}_merged_R2_001.fastq.gz"

    if [[ -f "$R1_FILE" && -f "$R2_FILE" ]]; then
        echo "Trimming: $BASE"
        trim_galore --paired --length 40 --max_n 1 -q 20 -j 8 \
            -o "/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/RawProcessing" \
            "$R1_FILE" "$R2_FILE"
    else
        echo "Missing file(s) for $BASE"
    fi
done

````

## 3. BWA mem merge reads
````
# Directory containing trimmed .gz files
input_dir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/RawProcessing"
output_dir="$input_dir/bams"
mkdir -p "$output_dir"

# Reference genome
reference_genome="/mnt/d/xQTL_2025_Data/ref/dm6.fa"

# Number of threads for BWA
threads=50

# Iterate over all R1 files in the specified directory
for file1 in "$input_dir"/*_merged_R1_*_val_1.fq.gz; do
    file2="${file1/_R1_/_R2_}"
    file2="${file2/_val_1/_val_2}"
    echo "$file1"
    echo "$file2"

    if [ -f "$file2" ]; then
        output_file1=$(basename "${file1%_R1_*.fq.gz}_temp_aligned.bam")
        output_file2=$(basename "${file1%_R1_*.fq.gz}_aligned.bam")
        bwa mem "$reference_genome" -M -t "$threads" -v 3 "$file1" "$file2" | samtools view -bS - > "$output_dir/$output_file1"
        samtools sort "$output_dir/$output_file1" -o "$output_dir/$output_file2"
        samtools index "$output_dir/$output_file2"
    else
        echo "Corresponding R2 file for $file1 not found."
    fi
done

````
