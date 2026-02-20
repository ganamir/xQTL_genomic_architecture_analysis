# xQTL_genomic_architecture_analysis

# DGRP & DSPR Processing Pipeline
## 1. Merge samples from different Lanes <<< DGRP ONLY!! >>> Skip to step 2 for DSPR
````
#!/bin/bash

# Run directly with: bash merge_lanes_simple.sh

SOURCE_DIR="/mnt/d/xQTL_2025_Data/DGRP_XQTL"
OUTPUT_DIR="/mnt/d/xQTL_2025_Data/DGRP_XQTL/trimmed_fastqs"

echo "=== Starting Lane Merge ==="
echo "Source: $SOURCE_DIR"
echo "Output: $OUTPUT_DIR"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Counter
TOTAL=0
SUCCESS=0
FAILED=0

# Loop through all unique samples
for L001_R1 in "$SOURCE_DIR"/2_[A-E]*_L001_R1_001.fastq.gz; do
    
    # Extract base name
    BASE=$(basename "$L001_R1")
    BASE=${BASE%%_L001_R1_001.fastq.gz}
    
    ((TOTAL++))
    
    echo "[$TOTAL] Processing: $BASE"
    
    # Define file paths
    L001_R1="$SOURCE_DIR/${BASE}_L001_R1_001.fastq.gz"
    L001_R2="$SOURCE_DIR/${BASE}_L001_R2_001.fastq.gz"
    L002_R1="$SOURCE_DIR/${BASE}_L002_R1_001.fastq.gz"
    L002_R2="$SOURCE_DIR/${BASE}_L002_R2_001.fastq.gz"
    
    # Check all files exist
    if [[ -f "$L001_R1" && -f "$L001_R2" && -f "$L002_R1" && -f "$L002_R2" ]]; then
        
        # Merge R1
        cat "$L001_R1" "$L002_R1" > "$OUTPUT_DIR/${BASE}_merged_R1_001.fastq.gz"
        
        # Merge R2
        cat "$L001_R2" "$L002_R2" > "$OUTPUT_DIR/${BASE}_merged_R2_001.fastq.gz"
        
        echo "  ✓ Complete"
        ((SUCCESS++))
    else
        echo "  ✗ ERROR: Missing files"
        ((FAILED++))
    fi
done

echo ""
echo "=== Summary ==="
echo "Total samples: $TOTAL"
echo "Successful: $SUCCESS"
echo "Failed: $FAILED"
echo ""
echo "Output location: $OUTPUT_DIR"
ls -lh "$OUTPUT_DIR"/*_merged_*.fastq.gz 2>/dev/null | wc -l | xargs echo "Merged files created:"


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

## 3. BWA mem merge reads, sort, and output bams
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

## 4. Remove Duplicates via dedup
````
#!/bin/bash
INPUT_DIR="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/RawProcessing/bams"
OUTPUT_DIR="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/RawProcessing/bams/DuplicateClear"

mkdir -p "$OUTPUT_DIR"

for BAM in "$INPUT_DIR"/*.bam; do
    BASENAME=$(basename "$BAM" .bam)
    FINAL_OUTPUT="$OUTPUT_DIR/${BASENAME}_dedup.bam"

    if [[ -f "$FINAL_OUTPUT" ]]; then
        echo "✅ Already processed: $FINAL_OUTPUT"
        continue
    fi

    echo "🔄 Processing $BAM..."

    samtools sort -n -o "$OUTPUT_DIR/${BASENAME}_namesorted.bam" "$BAM"
    samtools fixmate -m "$OUTPUT_DIR/${BASENAME}_namesorted.bam" "$OUTPUT_DIR/${BASENAME}_fixmate.bam"
    samtools sort -o "$OUTPUT_DIR/${BASENAME}_positionsorted.bam" "$OUTPUT_DIR/${BASENAME}_fixmate.bam"
    samtools markdup -r -s "$OUTPUT_DIR/${BASENAME}_positionsorted.bam" "$FINAL_OUTPUT" 2> "$OUTPUT_DIR/${BASENAME}_dedup_stats.txt"

    echo "✅ Output written to $FINAL_OUTPUT"

    rm "$OUTPUT_DIR/${BASENAME}_namesorted.bam"
    rm "$OUTPUT_DIR/${BASENAME}_fixmate.bam"
    rm "$OUTPUT_DIR/${BASENAME}_positionsorted.bam"
done
````

## 5. Add ReadGroups & Rename Files:
````
#!/bin/bash
indir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/Spino3"
outRG="${indir}/RG_groups"
mkdir -p "$outRG"

# Treatment samples A1-A8 (A6-A12 + B1)
declare -A t_samples
t_samples=([A6]=S120 [A7]=S124 [A8]=S128 [A9]=S132 [A10]=S136 [A11]=S140 [A12]=S144 [B1]=S98)
rg_num=1
for sample in A6 A7 A8 A9 A10 A11 A12 B1; do
    s=${t_samples[$sample]}
    inbam="2_${sample}_${s}_merged_aligned_dedup.bam"
    shortname=$(basename "$inbam" .bam)
    rgid="A${rg_num}"
    rgsm="A${rg_num}"
    echo "Processing $inbam -> RGID=$rgid"
    picard AddOrReplaceReadGroups \
        I="${indir}/${inbam}" \
        O="${outRG}/${shortname}.RGfixed.bam" \
        SORT_ORDER=coordinate \
        RGPL=illumina RGPU=D109LACXX RGLB=Lib1 \
        RGID=${rgid} RGSM=${rgsm} \
        VALIDATION_STRINGENCY=LENIENT
    mv "${outRG}/${shortname}.RGfixed.bam" "${outRG}/A${rg_num}.bam"
    ((rg_num++))
done

# Control samples C1-C8 (B2-B9)
declare -A c_samples
c_samples=([B2]=S103 [B3]=S108 [B4]=S113 [B5]=S117 [B6]=S121 [B7]=S125 [B8]=S129 [B9]=S133)
rg_num=1
for sample in B2 B3 B4 B5 B6 B7 B8 B9; do
    s=${c_samples[$sample]}
    inbam="2_${sample}_${s}_merged_aligned_dedup.bam"
    shortname=$(basename "$inbam" .bam)
    rgid="C${rg_num}"
    rgsm="C${rg_num}"
    echo "Processing $inbam -> RGID=$rgid"
    picard AddOrReplaceReadGroups \
        I="${indir}/${inbam}" \
        O="${outRG}/${shortname}.RGfixed.bam" \
        SORT_ORDER=coordinate \
        RGPL=illumina RGPU=D109LACXX RGLB=Lib1 \
        RGID=${rgid} RGSM=${rgsm} \
        VALIDATION_STRINGENCY=LENIENT
    mv "${outRG}/${shortname}.RGfixed.bam" "${outRG}/C${rg_num}.bam"
    ((rg_num++))
done


````

## 6. Merge Founders & Experimental Samples via mpileup
````
````



