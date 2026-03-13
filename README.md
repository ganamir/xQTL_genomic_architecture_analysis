# xQTL_genomic_architecture_analysis
# DGRP Founder Prep

## 1. Download DGRP LINES data (drgp_lines.tsv) from https://dgrpool.epfl.ch/

## 2. Filter DGRP Line data to contain Founder lines of interest (file: dgrp_mapping_clean)

## 3. Download Huang 2014 DGRP data from https://genome.cshlp.org/content/24/7/1193/suppl/DC1

## 4. Align out lines with 2014 resequenced lines & remove 454 sequencing data and only keep illumina (file: my_100_lines_illumina.csv)

## 5. Create sra_accessions file to download the DGRP_# SRA files (file: sra_accessions_final.txt)

## 6. Download the SRA files from NCBI using step 5 output:
<details>
<summary>Click to expand code</summary>

````
#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --output=logs/sra_%A_%a.out
#SBATCH --error=logs/sra_%A_%a.err
#SBATCH --array=1-126%14
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

ACCESSIONS="sra_accessions_final.txt"
SRA_DIR="sra_downloads"
FASTQ_DIR="fastq_folder"

mkdir -p "$SRA_DIR" "$FASTQ_DIR" logs

# Get this task's line
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCESSIONS")
DGRP=$(echo "$LINE" | cut -d',' -f1)
ACC=$(echo "$LINE"  | cut -d',' -f2 | tr -d '[:space:]')

echo "[INFO] Task $SLURM_ARRAY_TASK_ID: $DGRP - $ACC"

# Skip if already successfully downloaded
if ls "$FASTQ_DIR"/${DGRP}_*.fastq 1>/dev/null 2>&1; then
    echo "[SKIP] $DGRP - $ACC already exists"
    exit 0
fi

# Download into its own subdirectory to avoid race conditions
prefetch "$ACC" --output-directory "$SRA_DIR/$ACC"

# Find the SRA file ONLY within this accession's own directory
SRR_DIR=$(find "$SRA_DIR/$ACC" -name "*.sra" -o -name "*.sralite" | grep -v vdbcache | head -1 | xargs dirname 2>/dev/null)

if [ -z "$SRR_DIR" ]; then
    echo "[FAILED] prefetch did not produce an .sra file for $DGRP - $ACC"
    exit 1
fi

SRR=$(basename "$SRR_DIR")
echo "[INFO] Resolved $ACC -> $SRR"

# Extract fastq into a temp per-accession folder to avoid collisions
TMPFASTQ="$FASTQ_DIR/tmp_${ACC}"
mkdir -p "$TMPFASTQ"

fasterq-dump "$SRR_DIR" \
    --outdir "$TMPFASTQ" \
    --threads 4 \
    --progress

# Check fastq files were actually created
FASTQS=("$TMPFASTQ"/*.fastq)
if [ ${#FASTQS[@]} -eq 0 ] || [ ! -f "${FASTQS[0]}" ]; then
    echo "[FAILED] fasterq-dump produced no fastq files for $DGRP - $ACC"
    rm -rf "$TMPFASTQ"
    exit 1
fi

# Verify this SRR is not identical to an already-existing one
for f in "${FASTQS[@]}"; do
    FSIZE=$(stat -c%s "$f")
    FNAME=$(basename "$f")
    # Check if any existing fastq has the same size AND different DGRP name (potential duplicate)
    DUPE=$(find "$FASTQ_DIR" -maxdepth 1 -name "DGRP_*_${SRR}*.fastq" ! -path "$TMPFASTQ/*" 2>/dev/null | head -1)
    if [ -n "$DUPE" ]; then
        echo "[WARNING] $DGRP - $ACC ($SRR) appears to duplicate existing file: $DUPE"
        echo "[WARNING] This may indicate multiple SRX accessions map to the same SRR"
    fi
done

# Rename and move to final fastq folder: DGRP_21_SRR933566_1.fastq
for f in "${FASTQS[@]}"; do
    base=$(basename "$f")
    mv "$f" "$FASTQ_DIR/${DGRP}_${base}"
done

rm -rf "$TMPFASTQ"

# Final check
if ! ls "$FASTQ_DIR"/${DGRP}_*.fastq 1>/dev/null 2>&1; then
    echo "[FAILED] No fastq found after rename for $DGRP - $ACC"
    exit 1
fi

echo "[DONE] $DGRP - $ACC -> $SRR"

````

</details>


## 2. Trim files

<details>
<summary>Click to expand code</summary>

````
#!/bin/bash
#SBATCH --job-name=trim_galore
#SBATCH --output=logs/trim_%A_%a.out
#SBATCH --error=logs/trim_%A_%a.err
#SBATCH --array=1-100%5
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

FASTQ_DIR="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/fastq_folder"
TRIM_DIR="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/trimmed_fastqs"

mkdir -p "$TRIM_DIR" logs

DGRP=$(ls "$FASTQ_DIR"/DGRP_*.fastq | grep -oP 'DGRP_\d+' | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p")

echo "[INFO] Task $SLURM_ARRAY_TASK_ID: $DGRP"

# Skip if already trimmed
if ls "$TRIM_DIR"/${DGRP}_*trimmed* 1>/dev/null 2>&1; then
    echo "[SKIP] $DGRP already trimmed"
    exit 0
fi

# Paired-end runs
for R1 in "$FASTQ_DIR"/${DGRP}_*_1.fastq; do
    [ -f "$R1" ] || continue
    R2="${R1/_1.fastq/_2.fastq}"
    if [[ -f "$R2" ]]; then
        echo "[PE] Trimming $R1 + $R2"
        trim_galore --paired --length 30 --max_n 1 -q 20 -j 8 \
            -o "$TRIM_DIR" "$R1" "$R2"
    fi
done

# Single-end runs
for SE in "$FASTQ_DIR"/${DGRP}_*.fastq; do
    [[ "$SE" == *_1.fastq ]] && continue
    [[ "$SE" == *_2.fastq ]] && continue
    echo "[SE] Trimming $SE"
    trim_galore --length 30 --max_n 1 -q 20 -j 8 \
        -o "$TRIM_DIR" "$SE"
done

echo "[DONE] $DGRP"


````

</details>

## 3. BWA mem merge reads, sort, and output bams

<details>
<summary>Click to expand code</summary>

````
input_dir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/trimmed_fastqs"
output_dir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/bams"
mkdir -p "$output_dir"
reference_genome="/mnt/d/xQTL_2025_Data/ref/dm6.fa"
threads=50

# --- Paired-end ---
for file1 in "$input_dir"/DGRP_*_1_val_1.fq; do
    [ -f "$file1" ] || continue
    file2="${file1/_1_val_1.fq/_2_val_2.fq}"
    if [ ! -f "$file2" ]; then
        echo "R2 not found for $file1, skipping."
        continue
    fi

    base=$(basename "$file1" _1_val_1.fq)
    line_num=$(echo "$base" | cut -d_ -f2)
    rg="@RG\tID:${line_num}\tSM:${base}\tPL:ILLUMINA"

    temp_bam="$output_dir/${base}_temp_aligned.bam"
    out_bam="$output_dir/${base}_aligned.bam"

    echo "PE: $file1 + $file2"
    bwa mem -M -t "$threads" -v 3 -R "$rg" "$reference_genome" "$file1" "$file2" \
        | samtools view -bS - > "$temp_bam"
    samtools sort "$temp_bam" -o "$out_bam"
    samtools index "$out_bam"
    rm "$temp_bam"
done

# --- Single-end ---
for file1 in "$input_dir"/DGRP_*_trimmed.fq; do
    [ -f "$file1" ] || continue

    base=$(basename "$file1" _trimmed.fq)
    line_num=$(echo "$base" | cut -d_ -f2)
    rg="@RG\tID:${line_num}\tSM:${base}\tPL:ILLUMINA"

    temp_bam="$output_dir/${base}_temp_aligned.bam"
    out_bam="$output_dir/${base}_aligned.bam"

    echo "SE: $file1"
    bwa mem -M -t "$threads" -v 3 -R "$rg" "$reference_genome" "$file1" \
        | samtools view -bS - > "$temp_bam"
    samtools sort "$temp_bam" -o "$out_bam"
    samtools index "$out_bam"
    rm "$temp_bam"
done

````

</details>

## 4. Dedup & Merge
<details>
<summary>Click to expand code</summary>

````
#!/bin/bash
bam_dir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/bams"
out_dir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/final_bams"
mkdir -p "$out_dir"

# Get unique line numbers
ls "$bam_dir"/*_aligned.bam | grep -oP 'DGRP_\d+' | sort -u | while read line; do
    bams=("$bam_dir"/${line}_*_aligned.bam)
    deduped_bams=()

    # Dedup each SRR BAM separately
    for bam in "${bams[@]}"; do
        base=$(basename "$bam" .bam)
        deduped="$out_dir/${base}_dedup.bam"
        deduped_bams+=("$deduped")

        samtools sort -n -o "$out_dir/${base}_namesorted.bam" "$bam"
        samtools fixmate -m "$out_dir/${base}_namesorted.bam" "$out_dir/${base}_fixmate.bam"
        samtools sort -o "$out_dir/${base}_positionsorted.bam" "$out_dir/${base}_fixmate.bam"
        samtools markdup -r -s "$out_dir/${base}_positionsorted.bam" "$deduped" 2> "$out_dir/${base}_dedup_stats.txt"

        rm "$out_dir/${base}_namesorted.bam" "$out_dir/${base}_fixmate.bam" "$out_dir/${base}_positionsorted.bam"
    done

    # Merge deduped BAMs per line
    final="$out_dir/${line}.bam"
    if [ "${#deduped_bams[@]}" -gt 1 ]; then
        echo "Merging ${#deduped_bams[@]} deduped BAMs for $line"
        samtools merge -f "$final" "${deduped_bams[@]}"
        rm "${deduped_bams[@]}"
    else
        mv "${deduped_bams[0]}" "$final"
    fi

    samtools index "$final"
    echo "Done: $final"
done

````

</details>

## 5. Add RG
<details>
<summary>Click to expand code</summary>

````
#!/bin/bash
#SBATCH --job-name=add_RG
#SBATCH --output=/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/logs/RG_%a.out
#SBATCH --array=1-100%20
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

in_dir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/final_bams"
out_dir="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/final_bams_RG"
mkdir -p "$out_dir"

# Get the Nth BAM for this array task
bam=$(ls "$in_dir"/DGRP_*.bam | sed -n "${SLURM_ARRAY_TASK_ID}p")
line=$(basename "$bam" .bam)
line_num=$(echo "$line" | cut -d_ -f2)

picard AddOrReplaceReadGroups \
    I="$bam" \
    O="$out_dir/${line}.bam" \
    RGID=${line_num} \
    RGSM=DGRP_${line_num} \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGLB=Lib1 \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=LENIENT

samtools index "$out_dir/${line}.bam"
echo "Done: $line"

````

</details>


# DGRP & DSPR Processing Pipeline
## 1. Merge samples from different Lanes <<< DGRP ONLY!! >>> Skip to step 2 for DSPR

<details>
<summary>Click to expand code</summary>

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

</details>

## 2. Trim files

<details>
<summary>Click to expand code</summary>

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

</details>

## 3. BWA mem merge reads, sort, and output bams

<details>
<summary>Click to expand code</summary>

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

</details>

## 4. Remove Duplicates via dedup

<details>
<summary>Click to expand code</summary>

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

</details>

## 5. Add ReadGroups & Rename Files:

<details>
<summary>Click to expand code</summary>

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

</details>

# Processing xQTL samples through Haplotype Inferencing pipeline:

## 1. Combine Founders & Samples:

### First find RG names of the founders and input them into the haplotype.parameters.R file

<details>
<summary>Click to expand code</summary>

````
for bam in /mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/final_bams_RG/DGRP_*.bam; do
    samtools view -H "$bam" | grep '^@RG' | grep -oP 'SM:\S+' | sed 's/SM://'
done | sort -V | awk 'BEGIN{printf "c("} {printf sep"\""$0"\""; sep=","} END{print ")"}' \
> /mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/founder_RG.txt
````

# Then combine and create RefAlt Tables:

````
#!/bin/bash
#SBATCH --job-name=refalt
#SBATCH --output=/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/Spino3/RG_groups/haplotype_inf_final/process/logs/refalt_%a.out
#SBATCH --array=1-5
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G

ref="/mnt/d/xQTL_2025_Data/ref/dm6.fa"
bam_dir_founders="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/DGRP_Founder_100/final_bams_RG"
bam_dir_pool="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/Spino3/RG_groups"
output="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/Spino3/RG_groups/haplotype_inf_final/process"
mkdir -p "$output/logs"

ls "$bam_dir_founders"/DGRP_*.bam "$bam_dir_pool"/{A,C}*.bam > /tmp/all_bams.txt

chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "Processing $mychr"

bcf_out="${output}/calls.${mychr}.bcf"

bcftools mpileup -I -d 1000 -r "$mychr" -a "FORMAT/AD,FORMAT/DP" -f "$ref" -b /tmp/all_bams.txt \
    --threads 10 | \
    bcftools call -mv -Ob -o "$bcf_out"

bcftools index "$bcf_out"

echo -ne "CHROM\tPOS" > "${output}/RefAlt.${mychr}.txt"
bcftools query -l "$bcf_out" | awk '{printf("\tREF_%s\tALT_%s",$1,$1)}' >> "${output}/RefAlt.${mychr}.txt"
echo -ne "\n" >> "${output}/RefAlt.${mychr}.txt"

bcftools view -m2 -M2 -v snps -i 'QUAL>59' "$bcf_out" | \
    bcftools query -f '%CHROM\t%POS[\t%AD{0}\t%AD{1}]\n' | \
    grep -v '\.' >> "${output}/RefAlt.${mychr}.txt"

echo "Finished $mychr"

````


</details>

## 1.5. Debugging code for assesing h_cutoff, and various other problems associated with selecting options, such as wald test cut_off. 

<details>
<summary>Click to expand code</summary>

````
library(tidyverse)
library(limSolve)
library(data.table)
library(dtplyr)

mychr <- "chr2L"
myparfile <- "haplotype.parameters.R"
mydir <- "process/"

source(myparfile)

filein  <- paste0(mydir, "/RefAlt.",     mychr, ".txt")
rdsfile <- paste0(mydir, "/R.haps.",     mychr, ".rds")
fileout <- paste0(mydir, "/R.haps.",     mychr, ".out.rds")

cat("=== Chromosome:", mychr, "===\n")
cat("Input file:", filein, "\n")

# =============================================================================
# FUNCTIONS
# =============================================================================

# Estimate haplotype composition for all samples in a window
# Works on pre-pivoted wide-format data (no per-window pivot_wider needed)
est_hap <- function(spotsdf, df_wide) {
  chunk <- df_wide %>%
    filter(CHROM == spotsdf$CHROM &
             POS > spotsdf$start &
             POS < spotsdf$end)
  
  founder_cols <- chunk %>% select(all_of(founders))
  pool_names   <- intersect(names_in_bam, names(chunk))
  
  # Cluster once using founder data only (shared across all samples)
  f_mat <- as.matrix(founder_cols[complete.cases(founder_cols), ])
  if (nrow(f_mat) < 10 | ncol(f_mat) < 2) {
    return(tibble(
      sample = pool_names,
      Groups = list(NA), Haps = list(NA), 
      Err = list(NA), Names = list(NA)
    ))
  }
  Groups <- cutree(hclust(dist(t(f_mat))), h = h_cutoff)
  cat(sprintf("SNPs after complete.cases: %d, n_groups: %d, max_dist: %.2f\n", 
              nrow(f_mat), length(unique(Groups)), max(dist(t(f_mat)))))
  
  results <- lapply(pool_names, function(samp) {
    sampdf <- data.frame(freq = chunk[[samp]], founder_cols)
    est_hap2(sampdf, Groups)
  })
  
  tibble(
    sample = pool_names,
    Groups = lapply(results, `[[`, "Groups"),
    Haps   = lapply(results, `[[`, "Haps"),
    Err    = lapply(results, `[[`, "Err"),
    Names  = lapply(results, `[[`, "Names")
  )
}

# Estimate haplotype proportions for a single sample via lsei
# sampdf has columns: freq (pool), then one column per founder
est_hap2 <- function(sampdf, Groups) {
  founder_mat   <- sampdf %>% select(-freq)
  Y             <- sampdf$freq
  
  good          <- !is.na(Y) & complete.cases(founder_mat)
  Y             <- Y[good]
  founder_mat   <- founder_mat[good, ]
  m_founder_mat <- as.matrix(founder_mat)
  
  if (nrow(m_founder_mat) < 10 | ncol(m_founder_mat) < 2) {
    return(list(Groups = NA, Haps = NA, Err = NA, Names = NA))
  }
  
  d <- ncol(m_founder_mat)
  
  out <- lsei(
    A = m_founder_mat, B = Y,
    E = t(matrix(rep(1, d))), F = 1,
    G = diag(d), H = matrix(rep(0.0003, d)),
    verbose = FALSE, fulloutput = TRUE
  )
  
  # Collapse using the shared clustering
  unique_groups <- sort(unique(Groups))
  n_groups      <- length(unique_groups)
  
  T_mat <- matrix(0, nrow = n_groups, ncol = d)
  for (i in seq_along(unique_groups)) {
    T_mat[i, which(Groups == unique_groups[i])] <- 1
  }
  
  collapsed_haps <- as.vector(T_mat %*% out$X)
  collapsed_err  <- T_mat %*% out$cov %*% t(T_mat)
  group_names    <- paste0("G", unique_groups)
  names(collapsed_haps) <- group_names
  
  list(
    Groups = Groups,
    Haps   = collapsed_haps,
    Err    = collapsed_err,
    Names  = group_names
  )
}

# =============================================================================
# DATA LOADING
# =============================================================================

cat("Loading RefAlt file...\n")

df <- lazy_dt(read.table(filein, header = TRUE))

df2 <- df %>%
  pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
  mutate(RefAlt = str_sub(lab, 1, 3),
         name   = str_sub(lab, 5)) %>%
  select(-lab) %>%
  pivot_wider(names_from = RefAlt, values_from = count) %>%
  mutate(freq = REF / (REF + ALT),
         N    = REF + ALT) %>%
  select(-c("REF", "ALT")) %>%
  as_tibble()

rm(df)
cat("df2 loaded:", nrow(df2), "rows\n")

# =============================================================================
# FOUNDER FREQUENCY ROUNDING
# =============================================================================

df3 <- df2 %>%
  mutate(freq = case_when(
    name %in% founders & freq >  0.55 ~  1,
    name %in% founders & freq <  0.45 ~  0,
    name %in% founders & freq >= 0.45 & freq <= 0.55 ~ NA_real_,
    TRUE ~ freq
  )) %>%
  filter(!(name %in% founders & is.na(freq)))

rm(df2)
cat("df3 after founder rounding:", nrow(df3), "rows\n")

# =============================================================================
# SAVE PROCESSED DATA (long format, for downstream use)
# =============================================================================

saveRDS(df3, file = rdsfile)

# =============================================================================
# PIVOT TO WIDE FORMAT ONCE
# =============================================================================

cat("Pivoting to wide format (one-time operation)...\n")

df3_wide <- df3 %>%
  select(-N) %>%
  pivot_wider(names_from = name, values_from = freq)

cat("df3_wide:", nrow(df3_wide), "rows x", ncol(df3_wide), "cols\n")

# --- DIAGNOSTIC ---
test_chunk <- df3_wide %>% 
  filter(CHROM == mychr) %>% 
  slice(1:500)

founder_cols <- test_chunk %>% select(all_of(founders))
cat("Total SNPs in chunk:", nrow(founder_cols), "\n")
cat("SNPs after complete.cases:", sum(complete.cases(founder_cols)), "\n")

f_mat <- as.matrix(founder_cols[complete.cases(founder_cols), ])
if (nrow(f_mat) > 1) {
  d <- max(dist(t(f_mat)))
  cat("Max pairwise distance:", d, "\n")
  cat("Groups at h=7.5:", length(unique(cutree(hclust(dist(t(f_mat))), h = 7.5))), "\n")
}


# Test with an actual window, not just first 500 rows
mid_pos <- median(df3_wide$POS[df3_wide$CHROM == mychr])
test_chunk <- df3_wide %>% 
  filter(CHROM == mychr & POS > (mid_pos - size) & POS < (mid_pos + size))

founder_cols <- test_chunk %>% select(all_of(founders))
cat("Total SNPs in window:", nrow(founder_cols), "\n")
cat("SNPs after complete.cases:", sum(complete.cases(founder_cols)), "\n")

f_mat <- as.matrix(founder_cols[complete.cases(founder_cols), ])
d <- max(dist(t(f_mat)))
cat("Max pairwise distance:", d, "\n")

# Check the distribution of distances
all_d <- as.vector(dist(t(f_mat)))
cat("Distance quantiles:\n")
print(quantile(all_d, c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 1.0)))

# Try a few cutoffs
for (h in c(0.5, 1.0, 1.5, 2.0)) {
  g <- length(unique(cutree(hclust(dist(t(f_mat))), h = h)))
  cat(sprintf("h=%.1f: %d groups\n", h, g))
}

# --- END DIAGNOSTIC ---

rm(df3)
gc()

# =============================================================================
# DEFINE SCAN WINDOWS
# =============================================================================

minpos <- min(df3_wide$POS)
maxpos <- max(df3_wide$POS)
myseq  <- seq(0, maxpos, step)
myseq  <- myseq[myseq > minpos + size & myseq < maxpos - size]

spots <- data.frame(
  CHROM = rep(mychr, length(myseq)),
  pos   = myseq,
  start = myseq - size,
  end   = myseq + size
)

UU    <- unique(df3_wide$POS)
spots <- spots %>%
  rowwise() %>%
  mutate(NN = sum(start < UU) - sum(end < UU)) %>%
  filter(NN >= 50) %>%
  select(-NN)

cat(sprintf("Scan windows: %d\n", nrow(spots)))

# =============================================================================
# HAPLOTYPE SCAN (sequential)
# =============================================================================

total_windows <- nrow(spots)
cat(sprintf("Starting scan: %d windows, h_cutoff = %.1f\n", total_windows, h_cutoff))

t_start <- Sys.time()

spots2 <- spots %>%
  group_nest(row_number()) %>%
  mutate(out = imap(data, function(x, i) {
    if (i %% 100 == 0) {
      cat(sprintf("Progress: %d/%d windows (%.1f%%)\n",
                  i, total_windows, 100 * i / total_windows))
      flush.console()
    }
    est_hap(x, df3_wide)
  })) %>%
  unnest(data) %>%
  select(-c(start, end)) %>%
  select(-`row_number()`) %>%
  unnest_wider(out)

t_end <- Sys.time()
cat(sprintf("Scan finished in %.1f minutes\n", difftime(t_end, t_start, units = "mins")))

saveRDS(spots2, file = fileout)
cat("Output saved:", fileout, "\n")





# =============================================================================
# DEBUG HAPLOTYPEINF
# =============================================================================

chr2L <- read.table("process/Spino3/Spino3.pseudoscan.chr2R.txt")
filt_chr2L <- na.omit(chr2L)

design.df <- read.table("input_table.txt")

library(tidyverse)
library(limSolve)
library(abind)

filt_chr2L %>%
  ggplot(aes(x = pos, y = Pseu_log10p)) + geom_point()



#########
# Functions
#########

average_variance <- function(cov_matrix, tolerance = 1e-10) {
  n <- nrow(cov_matrix)  
  # Calculate eigenvalues
  eigenvalues <- eigen(cov_matrix, only.values = TRUE)$values  
  # Filter out eigenvalues that are effectively zero or negative
  positive_eigenvalues <- eigenvalues[eigenvalues > tolerance]  
  # Calculate the product of positive eigenvalues
  log_det <- sum(log(positive_eigenvalues))  
  # Use the number of positive eigenvalues for the root
  n_positive <- length(positive_eigenvalues)  
  # Calculate log of average variance
  log_avg_var <- log_det / n_positive  
  # Convert back to original scale
  avg_var <- exp(log_avg_var)  
  return(list(avg_var = avg_var, n_positive = n_positive, n_total = n))
}

wald.test3 = function(p1,p2,covar1,covar2,nrepl=1,N1=NA,N2=NA){
  
  # Wald test for multinomial frequencies
  # if nrepl = 1: (one replicate, analogous to chi square):
  #  p1 and p2 are vectors of relative frequencies to be compared
  # covar1 and covar2 are the reconstruction error 
  # covariance matrices from limSolve
  # the sampling covariance matrices are generated within limSolve
  # if nrepl > 1 (multiple replicates, analogous to CMH):
  #   p1 and p2 are matrices, each row is frequency vector for one replicate
  # covar1 and covar2 are tensors (3-dimensional arrays, third dimension 
  #  denotes replicate) for the linSolve covariance matrices
  # N1 (initial) and N2 (after treatment) 
  # are sample sizes, they are vectors when there is more than one replicate
  # N1[i], N2[i] are then for replicate i
  if (nrepl>1){
    N1.eff=rep(NA,nrepl)
    N2.eff=rep(NA,nrepl)
    lp1 = length(p1[1,])
    cv1=array(NA,c(lp1,lp1,nrepl))
    cv2=array(NA,c(lp1,lp1,nrepl))
    for (i in 1:nrepl){
      
      covmat1  = mn.covmat((N1[i]*p1[i,]+N2[i]*p2[i,])/(N1[i]+N2[i]),2*N1[i])
      covmat2  = mn.covmat((N1[i]*p1[i,]+N2[i]*p2[i,])/(N1[i]+N2[i]),2*N2[i])
      
      N1.eff[i] = sum(diag(covmat1))*4*N1[i]^2/(sum(diag(covmat1))*2*N1[i]+2*N1[i]*sum(diag(covar1[,,i])) )
      N2.eff[i] = sum(diag(covmat2))*4*N2[i]^2/(sum(diag(covmat2))*2*N2[i]+2*N2[i]*sum(diag(covar2[,,i])) )
      cv1[,,i]= (covmat1 + covar1[ , ,i])  * (N1.eff[i])^2
      cv2[,,i]= (covmat2 + covar2[ , ,i])  * (N2.eff[i])^2
      
    }
    
    p1 = N1.eff %*% p1 / sum(N1.eff)
    p2 = N2.eff %*% p2 / sum(N2.eff)
    
    covar1= rowSums(cv1, dims = 2) / sum(N1.eff)^2
    covar2= rowSums(cv2, dims = 2) / sum(N2.eff)^2
    # browser()
  }
  else {
    covmat1  = mn.covmat((N1*p1+N2*p2)/(N1+N2),2*N1)
    covmat2  = mn.covmat((N1*p1+N2*p2)/(N1+N2),2*N2)
    covar1 = covar1 + covmat1
    covar2 = covar2 + covmat2
  }
  
  df = length(p1)-1
  covar=covar1+covar2
  eg <- eigen(covar)
  # remove last eigenvector which corresponds to eigenvalue zero
  keep <- eg$values[1:df] > 1e-4  # threshold for meaningful eigenvalues
  ev <- eg$vectors[,1:df][,keep]
  eval <- eg$values[1:df][keep]
  df <- sum(keep)  # adjust degrees of freedom
  trafo<-diag(1/sqrt(eval)) %*% t(ev) 
  # set extremely small values to zero
  #new.covar[new.covar < 10^-9]=0
  p1= as.vector(p1); p2=as.vector(p2)
  tstat <- sum((trafo %*% (p1 - p2))^2)
  pval<- exp(pchisq(tstat,df,lower.tail=FALSE,log.p=TRUE))
  list(wald.test=tstat, p.value=pval, avg.var=average_variance(covar)$avg_var)
}

mn.covmat= function(p,n,min.p=0){
  # generate multinomial covariance matrix
  # p is vector of multinomial relative frequencies
  # n is sample size
  # compute covariance matrix for relative frequencies, for absolute frequencies multiply by n^2
  # if min.p >0, then values of p smaller than min.p are set to min.p and the resulting vector is rescaled.
  p[p<min.p] = min.p; p=p/sum(p)
  mat = - tcrossprod(p)
  diag(mat) = p*(1-p)
  mat = mat/n
  mat
}

pseudoN.test = function(p1,p2,covar1,covar2,nrepl,N1,N2){
  pseudoN_C = rep(NA,nrepl)
  pseudoN_Z = rep(NA,nrepl)
  for(i in 1:nrepl){
    pseudoN_C[i] = (2 * N1[i] * sum(p1[i,] * (1-p1[i]))) / (2 * N1[i] * sum(diag(covar1[,,i])) + sum(p1[i,] * (1-p1[i])))
    pseudoN_Z[i] = (2 * N2[i] * sum(p2[i,] * (1-p2[i]))) / (2 * N2[i] * sum(diag(covar2[,,i])) + sum(p2[i,] * (1-p2[i])))
  }
  Count1 = round(p1*pseudoN_C,0)
  Count2 = round(p2*pseudoN_Z,0)
  lowCountFounder = apply(rbind(Count1,Count2),2,sum)
  if(sum(lowCountFounder>=5)<2){
    log10p = NA
  }else{
    Count1 = Count1[,lowCountFounder >= 5]		
    Count2 = Count2[,lowCountFounder >= 5]		
    if(nrepl==1){
      out=chisq.test(rbind(Count1,Count2),correct=TRUE)
    }else{
      nF = ncol(Count1)
      tdf = data.frame(Count=c(as.numeric(t(Count1)),as.numeric(t(Count2))),
                       founder=rep(1:nF,2*nrepl),
                       TRT = c(rep(1,nF*nrepl),rep(2,nF*nrepl)),
                       REP = c(rep(1:nrepl,each=nF),rep(1:nrepl,each=nF)))
      D.x = xtabs(Count ~ founder + TRT + REP, data = tdf)
      out = mantelhaen.test(D.x,correct=TRUE)
    }
    log10p = -log10(out$p.value)
  }
  log10p
}

add_genetic = function(df){
  df$cM = rep(NA,nrow(df))
  fm=read.table("flymap.r6.txt",header=FALSE)
  colnames(fm)=c("chr","pos","cM")
  library(splines)
  for(chrs in c("chrX","chr2L","chr2R","chr3L","chr3R")){
    fmX = fm %>% filter(chr==chrs)
    out = ksmooth(fmX$pos,fmX$cM,kernel="normal",bandwidth=3e6)
    f_of_x = splinefun(out$x,out$y)
    temp = f_of_x(df$pos[df$chr==chrs])
    df$cM[df$chr==chrs] = temp
  }
  df
}

Heritability = function(p1, p2, nrepl, ProportionSelect, af_cutoff){
  nF = ncol(p1)
  tdf = data.frame(freq=c(as.numeric(t(p1)),as.numeric(t(p2))),
                   founder=rep(1:nF,2*nrepl),
                   TRT = c(rep("C",nF*nrepl),rep("Z",nF*nrepl)),
                   REP = c(rep(1:nrepl,each=nF),rep(1:nrepl,each=nF)))
  
  Falconer_H2 = tdf %>%
    pivot_wider(names_from = TRT, values_from = freq) %>%
    mutate(mean_diff_sq = (Z-C)^2) %>%
    mutate(mean_af_C = case_when(C <= af_cutoff ~ af_cutoff, .default = C)) %>%
    mutate(H2temp = mean_diff_sq/mean_af_C) %>%
    group_by(REP) %>%
    summarize(H2temp_sum = sum(H2temp)) %>%
    ungroup() %>%
    left_join(ProportionSelect,by="REP") %>%
    filter(!is.na(Proportion)) %>%
    mutate(Falcon_i = dnorm(qnorm(1-Proportion))/Proportion) %>%
    group_by(REP) %>%
    summarize(H2 = 200 * H2temp_sum / Falcon_i^2) %>%
    ungroup() %>%
    summarize(mH2 = mean(H2)) %>%
    pull(mH2)
  
  Cutler_H2 = tdf %>%
    pivot_wider(names_from = TRT, values_from = freq) %>%
    left_join(ProportionSelect,by="REP") %>%
    filter(!is.na(Proportion)) %>%
    mutate(Penetrance = (Z * Proportion)/C) %>%
    mutate(Penetrance = case_when(Penetrance <= Proportion/2 ~ Proportion/2,
                                  Penetrance >= 2*Proportion ~ 2*Proportion,
                                  .default = Penetrance)) %>% 
    mutate(Affect = qnorm(1-Proportion) - qnorm(1-Penetrance)) %>%
    mutate(marg_Va = Affect^2 * C) %>%
    group_by(REP) %>%
    mutate(H2 = 200*sum(marg_Va)) %>%
    ungroup() %>%
    summarize(mH2 = mean(H2)) %>%
    pull(mH2)
  
  list(Falconer_H2=Falconer_H2, Cutler_H2=Cutler_H2)
}

doscan = function(df,chr,Nfounders){
  sexlink = 1
  if(chr=="chrX"){ sexlink=0.75 }
  
  # I tested with xx2$data[[1]]
  df2 = df %>%
    unnest(cols = c(sample, Groups, Haps, Err, Names)) %>%
    left_join(recodeTable) %>%
    select(-sample) %>% mutate(sample=pool) %>% select(-pool) %>%
    left_join(Numflies, join_by(sample==pool)) %>%
    separate(sample,into=c("longTRT","REP","REPrep"),remove=FALSE) %>%
    left_join(TreatmentMapping)
  
  # only analyze data for which all founders are discernable..
  allFounders = as.numeric(df2 %>% mutate(mm = max(unlist(Groups))) %>% summarize(max(mm)))	
  
  ll = list(Wald_log10p = NA, Pseu_log10p = NA, Falc_H2 = NA, Cutl_H2 = NA, avg.var = NA)
  if(allFounders!=Nfounders){ return(ll) }
  
  ##  now cases where all founders are OK
  ##  now collapse any pure replicates.  This is tidy ugly.  But I feel there 
  ##  is value in keeping dataframe columns as lists...
  df3 = df2 %>%
    select(-Groups) %>%
    group_by(TRT,REP) %>%	
    summarise(Err_mean = list(reduce(map(Err, ~as.matrix(.x)), `+`)/length(Err)),
              Haps_mean = list(reduce(map(Haps, ~as.vector(.x)), `+`)/length(Haps)),
              Names = list(first(Names)),
              Num_mean = sexlink*mean(Num)) %>%
    rename(Haps=Haps_mean,Num=Num_mean,Err=Err_mean)
  
  ## these summaries of the data are pretty useful for tests
  p1 = df3 %>% filter(TRT=="C") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t() 
  row.names(p1) <- NULL
  p2 = df3 %>% filter(TRT=="Z") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
  row.names(p2) <- NULL
  covar1 = do.call(abind, c(df3 %>% filter(TRT=="C") %>% pull(Err), along = 3))
  covar2 = do.call(abind, c(df3 %>% filter(TRT=="Z") %>% pull(Err), along = 3))
  nrepl = df3 %>% filter(TRT=="C") %>% nrow()
  nrepl == df3 %>% filter(TRT=="Z") %>% nrow()
  N1 = df3 %>% filter(TRT=="C") %>% pull(Num)
  N2 = df3 %>% filter(TRT=="Z") %>% pull(Num)
  
  wt = tryCatch(
    wald.test3(p1,p2,covar1,covar2,nrepl,N1,N2),
    error = function(e) list(wald.test=NA, p.value=NA, avg.var=NA)
  )
  Wald_log10p = -log10(wt$p.value)
  Pseu_log10p = tryCatch(
    pseudoN.test(p1,p2,covar1,covar2,nrepl,N1,N2),
    error = function(e) NA
  )
  
  af_cutoff = 0.01     # 1% --- heritability estimators can be off for really low allele frequencies
  temp = Heritability(p1, p2, nrepl, ProportionSelect, af_cutoff)
  Falc_H2 = temp$Falconer_H2
  Cutl_H2 = temp$Cutler_H2
  
  ll = list(Wald_log10p = Wald_log10p, Pseu_log10p = Pseu_log10p,
            Falc_H2 = Falc_H2, Cutl_H2 = Cutl_H2, avg.var = wt$avg.var)
  ll
}

doscan2 = function(df,chr,Nfounders){
  sexlink = 1
  if(chr=="chrX"){ sexlink=0.75 }
  
  # I tested with xx2$data[[1]]
  df2 = df %>%
    unnest(cols = c(sample, Groups, Haps, Err, Names)) %>%
    left_join(design.df, join_by(sample==bam)) %>%
    filter(!is.na(TRT))
  
  # only analyze data for which all founders are discernable..
  allFounders = as.numeric(df2 %>% mutate(mm = max(unlist(Groups))) %>% summarize(max(mm)))	
  
  ll = list(Wald_log10p = NA, Pseu_log10p = NA, Falc_H2 = NA, Cutl_H2 = NA, avg.var = NA)
  if (is.na(allFounders) || allFounders < 2) { return(ll) }
  
  ##  now cases where all founders are OK
  ##  now collapse any pure replicates.  This is tidy ugly.  But I feel there 
  ##  is value in keeping dataframe columns as lists...
  df3 = df2 %>%
    select(-Groups) %>%
    group_by(TRT,REP) %>%	
    summarise(Err_mean = list(reduce(map(Err, ~as.matrix(.x)), `+`)/length(Err)),
              Haps_mean = list(reduce(map(Haps, ~as.vector(.x)), `+`)/length(Haps)),
              Names = list(first(Names)),
              Num_mean = sexlink*mean(Num)) %>%
    rename(Haps=Haps_mean,Num=Num_mean,Err=Err_mean)
  
  ## these summaries of the data are pretty useful for tests
  p1 = df3 %>% filter(TRT=="C") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t() 
  row.names(p1) <- NULL
  p2 = df3 %>% filter(TRT=="Z") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
  row.names(p2) <- NULL
  covar1 = do.call(abind, c(df3 %>% filter(TRT=="C") %>% pull(Err), along = 3))
  covar2 = do.call(abind, c(df3 %>% filter(TRT=="Z") %>% pull(Err), along = 3))
  nrepl = df3 %>% filter(TRT=="C") %>% nrow()
  nrepl == df3 %>% filter(TRT=="Z") %>% nrow()
  N1 = df3 %>% filter(TRT=="C") %>% pull(Num)
  N2 = df3 %>% filter(TRT=="Z") %>% pull(Num)
  
  wt = tryCatch(
    wald.test3(p1,p2,covar1,covar2,nrepl,N1,N2),
    error = function(e) list(wald.test=NA, p.value=NA, avg.var=NA)
  )
  Wald_log10p = -log10(wt$p.value)
  Pseu_log10p = tryCatch(
    pseudoN.test(p1,p2,covar1,covar2,nrepl,N1,N2),
    error = function(e) NA
  )
  
  af_cutoff = 0.01     # 1% --- heritability estimators can be off for really low allele frequencies
  temp = Heritability(p1, p2, nrepl, ProportionSelect, af_cutoff)
  Falc_H2 = temp$Falconer_H2
  Cutl_H2 = temp$Cutler_H2
  
  ll = list(Wald_log10p = Wald_log10p, Pseu_log10p = Pseu_log10p,
            Falc_H2 = Falc_H2, Cutl_H2 = Cutl_H2, avg.var = wt$avg.var)
  ll
}

xx1 = readRDS("process/R.haps.chr2L.out.rds")
Nfounders=length(xx1$Groups[[1]][[1]])
ProportionSelect = design.df %>% filter(TRT=="Z") %>% select(REP,Proportion) %>% arrange(REP)

bb1 = xx1 %>%
  head(n=100) %>%
  group_by(CHROM,pos) %>%
  nest() %>%
  mutate(out = map2(data, CHROM, doscan2, Nfounders=Nfounders)) %>%
  unnest_wider(out)
bb2 = bb1 %>% select(-data) %>% rename(chr=CHROM)
bb3 = add_genetic(bb2)

# ---- DIAGNOSTIC: inspect eigenvalues on a bad window ----
if (any(bb3$Wald_log10p > 200, na.rm = TRUE)) {
  bad_row <- bb3 %>% filter(Wald_log10p > 200) %>% slice(1)
  bad_pos <- bad_row$pos
  bad_chr <- bad_row$chr
  cat("Inspecting bad window at", bad_chr, "pos =", bad_pos, "\n")
  
  bad_data <- xx1 %>% filter(pos == bad_pos)
  
  df2_bad <- bad_data %>%
    unnest(cols = c(sample, Groups, Haps, Err, Names)) %>%
    left_join(design.df, join_by(sample==bam)) %>%
    filter(!is.na(TRT))
  
  df3_bad <- df2_bad %>%
    select(-Groups) %>%
    group_by(TRT,REP) %>%
    summarise(Err_mean = list(reduce(map(Err, ~as.matrix(.x)), `+`)/length(Err)),
              Haps_mean = list(reduce(map(Haps, ~as.vector(.x)), `+`)/length(Haps)),
              Names = list(first(Names)),
              Num_mean = mean(Num)) %>%
    rename(Haps=Haps_mean, Num=Num_mean, Err=Err_mean)
  
  p1_bad = df3_bad %>% filter(TRT=="C") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
  p2_bad = df3_bad %>% filter(TRT=="Z") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
  N1_bad = df3_bad %>% filter(TRT=="C") %>% pull(Num)
  N2_bad = df3_bad %>% filter(TRT=="Z") %>% pull(Num)
  covar1_bad = do.call(abind, c(df3_bad %>% filter(TRT=="C") %>% pull(Err), along = 3))
  covar2_bad = do.call(abind, c(df3_bad %>% filter(TRT=="Z") %>% pull(Err), along = 3))
  nrepl_bad = df3_bad %>% filter(TRT=="C") %>% nrow()
  
  N1.eff = N2.eff = rep(NA, nrepl_bad)
  lp1 = ncol(p1_bad)
  cv1 = cv2 = array(NA, c(lp1, lp1, nrepl_bad))
  for (i in 1:nrepl_bad) {
    covmat1 = mn.covmat((N1_bad[i]*p1_bad[i,]+N2_bad[i]*p2_bad[i,])/(N1_bad[i]+N2_bad[i]), 2*N1_bad[i])
    covmat2 = mn.covmat((N1_bad[i]*p1_bad[i,]+N2_bad[i]*p2_bad[i,])/(N1_bad[i]+N2_bad[i]), 2*N2_bad[i])
    N1.eff[i] = sum(diag(covmat1))*4*N1_bad[i]^2/(sum(diag(covmat1))*2*N1_bad[i]+2*N1_bad[i]*sum(diag(covar1_bad[,,i])))
    N2.eff[i] = sum(diag(covmat2))*4*N2_bad[i]^2/(sum(diag(covmat2))*2*N2_bad[i]+2*N2_bad[i]*sum(diag(covar2_bad[,,i])))
    cv1[,,i] = (covmat1 + covar1_bad[,,i]) * (N1.eff[i])^2
    cv2[,,i] = (covmat2 + covar2_bad[,,i]) * (N2.eff[i])^2
  }
  covar_total = rowSums(cv1, dims = 2)/sum(N1.eff)^2 + rowSums(cv2, dims = 2)/sum(N2.eff)^2
  
  evals = eigen(covar_total)$values
  cat("Number of clusters:", lp1, "\n")
  cat("All eigenvalues:\n")
  print(signif(sort(evals, decreasing = TRUE), 4))
  cat("\nEigenvalue range:", signif(max(evals), 4), "to", signif(min(evals), 4), "\n")
  stop("Diagnostic done - remove this block and set proper threshold")
} else {
  cat("No windows with Wald_log10p > 200 found\n")
}

# Add to the diagnostic, after covar_total:
cat("\nDiag of lsei covariance (reconstruction error):\n")
print(signif(diag(covar1_bad[,,1]), 4))
cat("\nDiag of multinomial covariance (sampling error):\n")
pooled_p <- (N1_bad[1]*p1_bad[1,]+N2_bad[1]*p2_bad[1,])/(N1_bad[1]+N2_bad[1])
print(signif(diag(mn.covmat(pooled_p, 2*N1_bad[1])), 4))



cat("\nHaps for Control rep 1:\n")
print(signif(p1_bad[1,], 4))
cat("\nHaps for Treatment rep 1:\n")
print(signif(p2_bad[1,], 4))
cat("\nDifference:\n")
print(signif(p2_bad[1,] - p1_bad[1,], 4))



# Check how many windows have out-of-range proportions
bad_count <- 0
for (i in 1:nrow(xx1)) {
  haps <- xx1$Haps[[i]]
  if (!is.list(haps)) next
  for (h in haps) {
    if (any(!is.na(h) & (h < -0.01 | h > 1.01))) {
      bad_count <- bad_count + 1
      break
    }
  }
}
cat("Windows with invalid proportions:", bad_count, "out of", nrow(xx1), "\n")



# I drop the loci for which the scan gives and NA
#bb4 = bb1 %>%
#	filter(!is.na(Pseu_log10p)) %>%
#	select(-c(Wald_log10p, Pseu_log10p, Falc_H2, Cutl_H2, avg.var, data)) %>%
#	left_join(xx1) %>%
#	select(-c(Err,Groups)) %>%
#	unnest(c(sample,Haps,Names)) %>%
#	unnest(c(Haps,Names)) %>%
#	rename(chr=CHROM,pool=sample,freq=Haps,founder=Names) %>%
#	left_join(design.df, by=c("pool"="bam")) %>%
#	select(c(chr,pos,TRT,REP,REPrep,freq,founder)) %>%
#	filter(!is.na(TRT)) %>%
#	group_by(chr,pos,TRT,REP,founder) %>%
#	summarize(freq=mean(freq,na.rm=TRUE))

write.table(bb3, fileout)
#write.table(bb4, fileout_meansBySample)

````

</details>

Because we have ~100 DGRP founders, we need to make sure that hierachical clustering works properly and creates a proper amount of haplotype blocks (x>1).

Number of haplotypes has to stay biologically relevant, while not being overdiscriminatory (<100)

## 2. Modify your input_table and haplotype.parameters files:
haplotype.parameters

<details>
<summary>Click to expand code</summary>

````
# JUICE Haplotype Parameters
# This file defines the parameters for haplotype estimation in the JUICE dataset

# Founder samples (these should match the column names in your REFALT data)
founders <- c("DGRP_21","DGRP_26","DGRP_31","DGRP_48","DGRP_69","DGRP_73","DGRP_83","DGRP_85","DGRP_88","DGRP_91","DGRP_100","DGRP_101","DGRP_105","DGRP_109","DGRP_129","DGRP_136","DGRP_142","DGRP_149","DGRP_153","DGRP_161","DGRP_176","DGRP_177","DGRP_181","DGRP_189","DGRP_217","DGRP_227","DGRP_228","DGRP_239","DGRP_287","DGRP_301","DGRP_306","DGRP_309","DGRP_310","DGRP_313","DGRP_315","DGRP_318","DGRP_321","DGRP_324","DGRP_332","DGRP_335","DGRP_338","DGRP_340","DGRP_350","DGRP_356","DGRP_359","DGRP_365","DGRP_370","DGRP_371","DGRP_373","DGRP_374","DGRP_375","DGRP_380","DGRP_381","DGRP_382","DGRP_385","DGRP_386","DGRP_390","DGRP_391","DGRP_399","DGRP_405","DGRP_406","DGRP_427","DGRP_437","DGRP_440","DGRP_441","DGRP_443","DGRP_502","DGRP_509","DGRP_513","DGRP_517","DGRP_535","DGRP_551","DGRP_559","DGRP_563","DGRP_566","DGRP_589","DGRP_634","DGRP_703","DGRP_705","DGRP_712","DGRP_732","DGRP_737","DGRP_748","DGRP_757","DGRP_776","DGRP_783","DGRP_786","DGRP_801","DGRP_804","DGRP_812","DGRP_819","DGRP_820","DGRP_821","DGRP_852","DGRP_857","DGRP_859","DGRP_861","DGRP_890","DGRP_900","DGRP_908")
#founders <- readLines("founder_names.txt")
# Step size for scanning positions along chromosome
step <- 1000  # 10kb steps (increased from 5kb for speed)

#try 1kb, 50kb, 100kb, 200kb

#maybe 10kb

# Window size - COMMENTED OUT, comes from command line
size <- 50000  # Window size comes from command line parameter

#try 100kb, 200kb, 300kb, 400kb, 500kb

# Clustering parameters
h_cutoff <- 2.5  # Default h_cutoff for fixed window hierarchical clustering

# Samples to process (allows selecting subset from large REFALT files)
names_in_bam=c("C1","C2","C3","C4","C5","C6","C7","C8","A1","A2","A3","A4","A5","A6","A7","A8")

````

</details>

input_table

<details>
<summary>Click to expand code</summary>

````
"bam" "REP" "REPrep" "Num" "Proportion" "TRT"
"1" "C1" 1 1 500 NA "C"
"2" "C2" 2 1 500 NA "C"
"3" "C3" 3 1 500 NA "C"
"4" "C4" 4 1 500 NA "C"
"5" "C5" 5 1 500 NA "C"
"6" "C6" 6 1 500 NA "C"
"7" "C7" 7 1 500 NA "C"
"8" "C8" 8 1 500 NA "C"
"9" "A1" 1 1 500 0.046 "Z"
"10" "A2" 2 1 500 0.056 "Z"
"11" "A3" 3 1 500 0.086 "Z"
"12" "A4" 4 1 500 0.076 "Z"
"13" "A5" 5 1 500 0.080 "Z"
"14" "A6" 6 1 500 0.063 "Z"
"15" "A7" 7 1 500 0.060 "Z"
"16" "A8" 8 1 500 0.033 "Z"
````

</details>

## 3. Run REFALT2haps.Andreas.sh (for DGRP, running a modified REFALT2haps.Andreas.code.r, reference below!)
Make sure that REFALT.chr . txt are located in the process folder
````
sbatch REFALT2haps.Andreas.sh haplotype.parameters.R "process/
````
These scripts have to be located in the same directory as the shell exectuble
````
REFALT2haps.Andreas.sh <<< You run this
REFALT2haps.Andreas.r <<< This gets piped in .sh & executed code.r
REFALT2haps.Andreas.code.r <<< This is the statistic side of the script
````

CONCERNS!!!!! <<<< DGRP ONLY >>>>

Founders data: really low coverage ~5-50x which creates a multitude of problems in REFALT2haps step.
1. Code is expecting that all founders have at least 1 read for that specific SNP. Not the biggest contributor to a problem.

2. Code is expecting homozygosity at every site, ALL founder samples have to be <0.03 or >0.97. Which is problematic in founder samples that are low coverage (DGRP RILs are ~11x coverage). Which leads to a lot of samples actually being non-homozygotic, i.e. somewhere between 0.03 and 0.97. Which gets filtered out (this is where most of the sites are getting filtered out)
<img width="852" height="414" alt="image" src="https://github.com/user-attachments/assets/2851764b-e950-47bc-962c-54f9ad737f08" />

3. Some sites are monomorphic between ALL founder samples (AF = 1, or AF = 0), which is bad for haplotype cluster estimation, how do you determine which haplotype it belongs to? So that get's filtered out.

In the end, a ~500,000 SNP chr2L ends up with ~1621 SNPs across the entire genome post filtering.


### Modified DGRP Scripts:
REFALT2haps.Andreas.code.r
Key modification Made:

1. Binerization of Founder AF (frequencies >0.55 → 1, <0.45 → 0, ambiguous 0.45–0.55 → NA. This reflects the biological expectation that inbred lines should be homozygous.)

2. Hierarchical clustering Cut off Testing & finding appropriate clustering distance (h_cutoff <- 1.5 DGRP | h_cutoff <- 2.5 DSPR)

3. Added clustering of similar haplotypes (was not done in original xQTL Anthony Long pipeline), now clusters are calculated and grouped via h_cutoff (lower h_cutoff means smaller cut off to define dissimilarity (Less number of SNP different between founders))

4. Clustering occurs before lsei runs, removes calculations errors associated with complex matrix with large number of founder samples (that are very similar to one another across the entire window)

<details>
<summary>Click to expand code</summary>

````
# =============================================================================
# FUNCTIONS
# =============================================================================

# Estimate haplotype composition for all samples in a window
# Works on pre-pivoted wide-format data (no per-window pivot_wider needed)
est_hap <- function(spotsdf, df_wide) {
  chunk <- df_wide %>%
    filter(CHROM == spotsdf$CHROM &
             POS > spotsdf$start &
             POS < spotsdf$end)
  
  founder_cols <- chunk %>% select(all_of(founders))
  pool_names   <- intersect(names_in_bam, names(chunk))
  
  # Cluster once using founder data only (shared across all samples)
  f_mat <- as.matrix(founder_cols[complete.cases(founder_cols), ])
  if (nrow(f_mat) < 10 | ncol(f_mat) < 2) {
    return(tibble(
      sample = pool_names,
      Groups = list(NA), Haps = list(NA), 
      Err = list(NA), Names = list(NA)
    ))
  }
  Groups <- cutree(hclust(dist(t(f_mat))), h = h_cutoff)
  
  results <- lapply(pool_names, function(samp) {
    sampdf <- data.frame(freq = chunk[[samp]], founder_cols)
    est_hap2(sampdf, Groups)
  })
  
  tibble(
    sample = pool_names,
    Groups = lapply(results, `[[`, "Groups"),
    Haps   = lapply(results, `[[`, "Haps"),
    Err    = lapply(results, `[[`, "Err"),
    Names  = lapply(results, `[[`, "Names")
  )
}

# Estimate haplotype proportions for a single sample via lsei
# sampdf has columns: freq (pool), then one column per founder
est_hap2 <- function(sampdf, Groups) {
  founder_mat   <- sampdf %>% select(-freq)
  Y             <- sampdf$freq
  
  good          <- !is.na(Y) & complete.cases(founder_mat)
  Y             <- Y[good]
  founder_mat   <- founder_mat[good, ]
  m_founder_mat <- as.matrix(founder_mat)
  
  if (nrow(m_founder_mat) < 10 | ncol(m_founder_mat) < 2) {
    return(list(Groups = NA, Haps = NA, Err = NA, Names = NA))
  }
  
  # Collapse founder matrix into cluster-level matrix BEFORE lsei
  unique_groups <- sort(unique(Groups))
  n_groups      <- length(unique_groups)
  
  cluster_mat <- sapply(unique_groups, function(cl) {
    members <- which(Groups == cl)
    if (length(members) == 1) {
      m_founder_mat[, members]
    } else {
      rowMeans(m_founder_mat[, members])
    }
  })
  
  d <- ncol(cluster_mat)
  
  if (d < 2) {
    return(list(Groups = NA, Haps = NA, Err = NA, Names = NA))
  }
  
  out <- lsei(
    A = cluster_mat, B = Y,
    E = t(matrix(rep(1, d))), F = 1,
    G = diag(d), H = matrix(rep(0.0003, d)),
    verbose = FALSE, fulloutput = TRUE
  )
  
  group_names <- paste0("G", unique_groups)
  names(out$X) <- group_names
  
  list(
    Groups = Groups,
    Haps   = out$X,
    Err    = out$cov,
    Names  = group_names
  )
}

# =============================================================================
# DATA LOADING
# =============================================================================

cat("Loading RefAlt file...\n")

df <- lazy_dt(read.table(filein, header = TRUE))

df2 <- df %>%
  pivot_longer(c(-CHROM, -POS), names_to = "lab", values_to = "count") %>%
  mutate(RefAlt = str_sub(lab, 1, 3),
         name   = str_sub(lab, 5)) %>%
  select(-lab) %>%
  pivot_wider(names_from = RefAlt, values_from = count) %>%
  mutate(freq = REF / (REF + ALT),
         N    = REF + ALT) %>%
  select(-c("REF", "ALT")) %>%
  as_tibble()

rm(df)
cat("df2 loaded:", nrow(df2), "rows\n")

# =============================================================================
# FOUNDER FREQUENCY ROUNDING
# =============================================================================

df3 <- df2 %>%
  mutate(freq = case_when(
    name %in% founders & freq >  0.55 ~  1,
    name %in% founders & freq <  0.45 ~  0,
    name %in% founders & freq >= 0.45 & freq <= 0.55 ~ NA_real_,
    TRUE ~ freq
  )) %>%
  filter(!(name %in% founders & is.na(freq)))

rm(df2)
cat("df3 after founder rounding:", nrow(df3), "rows\n")

# =============================================================================
# SAVE PROCESSED DATA (long format, for downstream use)
# =============================================================================

saveRDS(df3, file = rdsfile)

# =============================================================================
# PIVOT TO WIDE FORMAT ONCE
# =============================================================================

cat("Pivoting to wide format (one-time operation)...\n")

df3_wide <- df3 %>%
  select(-N) %>%
  pivot_wider(names_from = name, values_from = freq)

cat("df3_wide:", nrow(df3_wide), "rows x", ncol(df3_wide), "cols\n")

rm(df3)
gc()

# =============================================================================
# DEFINE SCAN WINDOWS
# =============================================================================

minpos <- min(df3_wide$POS)
maxpos <- max(df3_wide$POS)
myseq  <- seq(0, maxpos, step)
myseq  <- myseq[myseq > minpos + size & myseq < maxpos - size]

spots <- data.frame(
  CHROM = rep(mychr, length(myseq)),
  pos   = myseq,
  start = myseq - size,
  end   = myseq + size
)

UU    <- unique(df3_wide$POS)
spots <- spots %>%
  rowwise() %>%
  mutate(NN = sum(start < UU) - sum(end < UU)) %>%
  filter(NN >= 50) %>%
  select(-NN)

cat(sprintf("Scan windows: %d\n", nrow(spots)))

# =============================================================================
# HAPLOTYPE SCAN (sequential)
# =============================================================================

total_windows <- nrow(spots)
cat(sprintf("Starting scan: %d windows, h_cutoff = %.1f\n", total_windows, h_cutoff))

t_start <- Sys.time()

spots2 <- spots %>%
  group_nest(row_number()) %>%
  mutate(out = imap(data, function(x, i) {
    if (i %% 100 == 0) {
      cat(sprintf("Progress: %d/%d windows (%.1f%%)\n",
                  i, total_windows, 100 * i / total_windows))
      flush.console()
    }
    est_hap(x, df3_wide)
  })) %>%
  unnest(data) %>%
  select(-c(start, end)) %>%
  select(-`row_number()`) %>%
  unnest_wider(out)

t_end <- Sys.time()
cat(sprintf("Scan finished in %.1f minutes\n", difftime(t_end, t_start, units = "mins")))

saveRDS(spots2, file = fileout)
cat("Output saved:", fileout, "\n")

````

</details>

## 4. Run haps2scan.Apr2025
````
sbatch haps2scan.Apr2025.sh input_table.txt "process/" "Spino3"

````
These scripts have to be located in the same directroy as the shell executable
````
haps2scan.Apr2025.sh <<< You run this
haps2scan.Apr2025.r <<< This gets piped in .sh & executed code.r
haps2scan.Apr2025.code.r <<< Haplotype Inf. stats
scan_functions.r <<< More Haplotype Inf. stats
````

### Modified DGRP Scripts:
haps2scan.Apr2025.code
Key modification Made:

1. doscan2 works with NAs that are introduced in complex matrix calculations

2. Commented out founder contribution dataset. Due to clustering of haplotypes, we no longer are able to identify which founder set contributed to the window, which is fine for DGRP data.

<details>
<summary>Click to expand code</summary>

````
doscan2 = function(df,chr,Nfounders){
  sexlink = 1
  if(chr=="chrX"){ sexlink=0.75 }
  
  # I tested with xx2$data[[1]]
  df2 = df %>%
    unnest(cols = c(sample, Groups, Haps, Err, Names)) %>%
    left_join(design.df, join_by(sample==bam)) %>%
    filter(!is.na(TRT))
  
  # only analyze data for which all founders are discernable..
  allFounders = as.numeric(df2 %>% mutate(mm = max(unlist(Groups))) %>% summarize(max(mm)))	
  
  ll = list(Wald_log10p = NA, Pseu_log10p = NA, Falc_H2 = NA, Cutl_H2 = NA, avg.var = NA)
  if (is.na(allFounders) || allFounders < 2) { return(ll) }
  
  ##  now cases where all founders are OK
  ##  now collapse any pure replicates.  This is tidy ugly.  But I feel there 
  ##  is value in keeping dataframe columns as lists...
  df3 = df2 %>%
    select(-Groups) %>%
    group_by(TRT,REP) %>%	
    summarise(Err_mean = list(reduce(map(Err, ~as.matrix(.x)), `+`)/length(Err)),
              Haps_mean = list(reduce(map(Haps, ~as.vector(.x)), `+`)/length(Haps)),
              Names = list(first(Names)),
              Num_mean = sexlink*mean(Num)) %>%
    rename(Haps=Haps_mean,Num=Num_mean,Err=Err_mean)
  
  ## these summaries of the data are pretty useful for tests
  p1 = df3 %>% filter(TRT=="C") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t() 
  row.names(p1) <- NULL
  p2 = df3 %>% filter(TRT=="Z") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
  row.names(p2) <- NULL
  covar1 = do.call(abind, c(df3 %>% filter(TRT=="C") %>% pull(Err), along = 3))
  covar2 = do.call(abind, c(df3 %>% filter(TRT=="Z") %>% pull(Err), along = 3))
  nrepl = df3 %>% filter(TRT=="C") %>% nrow()
  nrepl == df3 %>% filter(TRT=="Z") %>% nrow()
  N1 = df3 %>% filter(TRT=="C") %>% pull(Num)
  N2 = df3 %>% filter(TRT=="Z") %>% pull(Num)
  
  wt = tryCatch(
    wald.test3(p1,p2,covar1,covar2,nrepl,N1,N2),
    error = function(e) list(wald.test=NA, p.value=NA, avg.var=NA)
  )
  Wald_log10p = -log10(wt$p.value)
  Pseu_log10p = tryCatch(
    pseudoN.test(p1,p2,covar1,covar2,nrepl,N1,N2),
    error = function(e) NA
  )
  
  af_cutoff = 0.01     # 1% --- heritability estimators can be off for really low allele frequencies
  temp = Heritability(p1, p2, nrepl, ProportionSelect, af_cutoff)
  Falc_H2 = temp$Falconer_H2
  Cutl_H2 = temp$Cutler_H2
  
  ll = list(Wald_log10p = Wald_log10p, Pseu_log10p = Pseu_log10p,
            Falc_H2 = Falc_H2, Cutl_H2 = Cutl_H2, avg.var = wt$avg.var)
  ll
}

xx1 = readRDS(filein)
Nfounders=length(xx1$Groups[[1]][[1]])
ProportionSelect = design.df %>% filter(TRT=="Z") %>% select(REP,Proportion) %>% arrange(REP)

bb1 = xx1 %>%
#	head(n=100) %>%
	group_by(CHROM,pos) %>%
	nest() %>%
	mutate(out = map2(data, CHROM, doscan2, Nfounders=Nfounders)) %>%
	unnest_wider(out)
bb2 = bb1 %>% select(-data) %>% rename(chr=CHROM)
bb3 = add_genetic(bb2)

# I drop the loci for which the scan gives and NA
#bb4 = bb1 %>%
#	filter(!is.na(Pseu_log10p)) %>%
#	select(-c(Wald_log10p, Pseu_log10p, Falc_H2, Cutl_H2, avg.var, data)) %>%
#	left_join(xx1) %>%
#	select(-c(Err,Groups)) %>%
#	unnest(c(sample,Haps,Names)) %>%
#	unnest(c(Haps,Names)) %>%
#	rename(chr=CHROM,pool=sample,freq=Haps,founder=Names) %>%
#	left_join(design.df, by=c("pool"="bam")) %>%
#	select(c(chr,pos,TRT,REP,REPrep,freq,founder)) %>%
#	filter(!is.na(TRT)) %>%
#	group_by(chr,pos,TRT,REP,founder) %>%
#	summarize(freq=mean(freq,na.rm=TRUE))

write.table(bb3, fileout)
#write.table(bb4, fileout_meansBySample)

````

</details>

### Modified DGRP Scripts:
scan_functions
Key modification Made:

1. Wald.test3 now has eigenvalue filtering, with many founders being so similar to one another, haplotype contributions coming from all of them are small, therefore we filter values that are smaller than e-4. (Helps to minimize NAs calculated from overly complex matrix & removes uneeded data that wouldn't have helped with inferences. Some signal might be lost, but without filtering, no signal would even be generated.)

NOTE: REVISIT AFTER FOUNDER REWORK!!!!!!!!!!!


<details>
<summary>Click to expand code</summary>

````
#########
# Functions
#########

average_variance <- function(cov_matrix, tolerance = 1e-10) {
  n <- nrow(cov_matrix)  
  # Calculate eigenvalues
  eigenvalues <- eigen(cov_matrix, only.values = TRUE)$values  
  # Filter out eigenvalues that are effectively zero or negative
  positive_eigenvalues <- eigenvalues[eigenvalues > tolerance]  
  # Calculate the product of positive eigenvalues
  log_det <- sum(log(positive_eigenvalues))  
  # Use the number of positive eigenvalues for the root
  n_positive <- length(positive_eigenvalues)  
  # Calculate log of average variance
  log_avg_var <- log_det / n_positive  
  # Convert back to original scale
  avg_var <- exp(log_avg_var)  
  return(list(avg_var = avg_var, n_positive = n_positive, n_total = n))
}

wald.test3 = function(p1,p2,covar1,covar2,nrepl=1,N1=NA,N2=NA){
    
    # Wald test for multinomial frequencies
    # if nrepl = 1: (one replicate, analogous to chi square):
    #  p1 and p2 are vectors of relative frequencies to be compared
    # covar1 and covar2 are the reconstruction error 
    # covariance matrices from limSolve
    # the sampling covariance matrices are generated within limSolve
    # if nrepl > 1 (multiple replicates, analogous to CMH):
    #   p1 and p2 are matrices, each row is frequency vector for one replicate
    # covar1 and covar2 are tensors (3-dimensional arrays, third dimension 
    #  denotes replicate) for the linSolve covariance matrices
    # N1 (initial) and N2 (after treatment) 
    # are sample sizes, they are vectors when there is more than one replicate
    # N1[i], N2[i] are then for replicate i
    if (nrepl>1){
      N1.eff=rep(NA,nrepl)
      N2.eff=rep(NA,nrepl)
      lp1 = length(p1[1,])
      cv1=array(NA,c(lp1,lp1,nrepl))
      cv2=array(NA,c(lp1,lp1,nrepl))
      for (i in 1:nrepl){
          
        covmat1  = mn.covmat((N1[i]*p1[i,]+N2[i]*p2[i,])/(N1[i]+N2[i]),2*N1[i])
        covmat2  = mn.covmat((N1[i]*p1[i,]+N2[i]*p2[i,])/(N1[i]+N2[i]),2*N2[i])
    
        N1.eff[i] = sum(diag(covmat1))*4*N1[i]^2/(sum(diag(covmat1))*2*N1[i]+2*N1[i]*sum(diag(covar1[,,i])) )
        N2.eff[i] = sum(diag(covmat2))*4*N2[i]^2/(sum(diag(covmat2))*2*N2[i]+2*N2[i]*sum(diag(covar2[,,i])) )
        cv1[,,i]= (covmat1 + covar1[ , ,i])  * (N1.eff[i])^2
        cv2[,,i]= (covmat2 + covar2[ , ,i])  * (N2.eff[i])^2
        
      }
        
      p1 = N1.eff %*% p1 / sum(N1.eff)
      p2 = N2.eff %*% p2 / sum(N2.eff)
     
      covar1= rowSums(cv1, dims = 2) / sum(N1.eff)^2
      covar2= rowSums(cv2, dims = 2) / sum(N2.eff)^2
     # browser()
    }
    else {
      covmat1  = mn.covmat((N1*p1+N2*p2)/(N1+N2),2*N1)
      covmat2  = mn.covmat((N1*p1+N2*p2)/(N1+N2),2*N2)
      covar1 = covar1 + covmat1
      covar2 = covar2 + covmat2
    }
  
  df = length(p1)-1
  covar=covar1+covar2
  eg <- eigen(covar)
  # remove last eigenvector which corresponds to eigenvalue zero
  keep <- eg$values[1:df] > 1e-4  # threshold for meaningful eigenvalues
  ev <- eg$vectors[,1:df][,keep]
  eval <- eg$values[1:df][keep]
  df <- sum(keep)  # adjust degrees of freedom
  trafo<-diag(1/sqrt(eval)) %*% t(ev) 
  # set extremely small values to zero
  #new.covar[new.covar < 10^-9]=0
  p1= as.vector(p1); p2=as.vector(p2)
  tstat <- sum((trafo %*% (p1 - p2))^2)
  pval<- exp(pchisq(tstat,df,lower.tail=FALSE,log.p=TRUE))
  list(wald.test=tstat, p.value=pval, avg.var=average_variance(covar)$avg_var)
}

mn.covmat= function(p,n,min.p=0){
  # generate multinomial covariance matrix
  # p is vector of multinomial relative frequencies
  # n is sample size
  # compute covariance matrix for relative frequencies, for absolute frequencies multiply by n^2
  # if min.p >0, then values of p smaller than min.p are set to min.p and the resulting vector is rescaled.
  p[p<min.p] = min.p; p=p/sum(p)
  mat = - tcrossprod(p)
  diag(mat) = p*(1-p)
  mat = mat/n
  mat
}
    
pseudoN.test = function(p1,p2,covar1,covar2,nrepl,N1,N2){
	pseudoN_C = rep(NA,nrepl)
	pseudoN_Z = rep(NA,nrepl)
	for(i in 1:nrepl){
		pseudoN_C[i] = (2 * N1[i] * sum(p1[i,] * (1-p1[i]))) / (2 * N1[i] * sum(diag(covar1[,,i])) + sum(p1[i,] * (1-p1[i])))
		pseudoN_Z[i] = (2 * N2[i] * sum(p2[i,] * (1-p2[i]))) / (2 * N2[i] * sum(diag(covar2[,,i])) + sum(p2[i,] * (1-p2[i])))
		}
	Count1 = round(p1*pseudoN_C,0)
	Count2 = round(p2*pseudoN_Z,0)
	lowCountFounder = apply(rbind(Count1,Count2),2,sum)
	if(sum(lowCountFounder>=5)<2){
		log10p = NA
		}else{
		Count1 = Count1[,lowCountFounder >= 5]		
		Count2 = Count2[,lowCountFounder >= 5]		
		if(nrepl==1){
			out=chisq.test(rbind(Count1,Count2),correct=TRUE)
			}else{
			nF = ncol(Count1)
			tdf = data.frame(Count=c(as.numeric(t(Count1)),as.numeric(t(Count2))),
				founder=rep(1:nF,2*nrepl),
				TRT = c(rep(1,nF*nrepl),rep(2,nF*nrepl)),
				REP = c(rep(1:nrepl,each=nF),rep(1:nrepl,each=nF)))
			D.x = xtabs(Count ~ founder + TRT + REP, data = tdf)
			out = mantelhaen.test(D.x,correct=TRUE)
			}
		log10p = -log10(out$p.value)
		}
	log10p
	}
        
add_genetic = function(df){
	df$cM = rep(NA,nrow(df))
	fm=read.table("flymap.r6.txt",header=FALSE)
	colnames(fm)=c("chr","pos","cM")
	library(splines)
	for(chrs in c("chrX","chr2L","chr2R","chr3L","chr3R")){
		fmX = fm %>% filter(chr==chrs)
		out = ksmooth(fmX$pos,fmX$cM,kernel="normal",bandwidth=3e6)
		f_of_x = splinefun(out$x,out$y)
		temp = f_of_x(df$pos[df$chr==chrs])
		df$cM[df$chr==chrs] = temp
		}
	df
	}

Heritability = function(p1, p2, nrepl, ProportionSelect, af_cutoff){
	nF = ncol(p1)
	tdf = data.frame(freq=c(as.numeric(t(p1)),as.numeric(t(p2))),
		founder=rep(1:nF,2*nrepl),
		TRT = c(rep("C",nF*nrepl),rep("Z",nF*nrepl)),
		REP = c(rep(1:nrepl,each=nF),rep(1:nrepl,each=nF)))

	Falconer_H2 = tdf %>%
		pivot_wider(names_from = TRT, values_from = freq) %>%
		mutate(mean_diff_sq = (Z-C)^2) %>%
		mutate(mean_af_C = case_when(C <= af_cutoff ~ af_cutoff, .default = C)) %>%
		mutate(H2temp = mean_diff_sq/mean_af_C) %>%
		group_by(REP) %>%
		summarize(H2temp_sum = sum(H2temp)) %>%
		ungroup() %>%
		left_join(ProportionSelect,by="REP") %>%
		filter(!is.na(Proportion)) %>%
		mutate(Falcon_i = dnorm(qnorm(1-Proportion))/Proportion) %>%
		group_by(REP) %>%
		summarize(H2 = 200 * H2temp_sum / Falcon_i^2) %>%
		ungroup() %>%
		summarize(mH2 = mean(H2)) %>%
		pull(mH2)
			
	Cutler_H2 = tdf %>%
		pivot_wider(names_from = TRT, values_from = freq) %>%
		left_join(ProportionSelect,by="REP") %>%
		filter(!is.na(Proportion)) %>%
		mutate(Penetrance = (Z * Proportion)/C) %>%
		mutate(Penetrance = case_when(Penetrance <= Proportion/2 ~ Proportion/2,
						  Penetrance >= 2*Proportion ~ 2*Proportion,
						  .default = Penetrance)) %>% 
		mutate(Affect = qnorm(1-Proportion) - qnorm(1-Penetrance)) %>%
		mutate(marg_Va = Affect^2 * C) %>%
		group_by(REP) %>%
		mutate(H2 = 200*sum(marg_Va)) %>%
		ungroup() %>%
		summarize(mH2 = mean(H2)) %>%
		pull(mH2)

	list(Falconer_H2=Falconer_H2, Cutler_H2=Cutler_H2)
	}

doscan = function(df,chr,Nfounders){
  sexlink = 1
  if(chr=="chrX"){ sexlink=0.75 }
  
  # I tested with xx2$data[[1]]
  df2 = df %>%
    unnest(cols = c(sample, Groups, Haps, Err, Names)) %>%
    left_join(recodeTable) %>%
    select(-sample) %>% mutate(sample=pool) %>% select(-pool) %>%
    left_join(Numflies, join_by(sample==pool)) %>%
    separate(sample,into=c("longTRT","REP","REPrep"),remove=FALSE) %>%
    left_join(TreatmentMapping)
  
  # only analyze data for which all founders are discernable..
  allFounders = as.numeric(df2 %>% mutate(mm = max(unlist(Groups))) %>% summarize(max(mm)))	
  
  ll = list(Wald_log10p = NA, Pseu_log10p = NA, Falc_H2 = NA, Cutl_H2 = NA, avg.var = NA)
  if(allFounders!=Nfounders){ return(ll) }
  
  ##  now cases where all founders are OK
  ##  now collapse any pure replicates.  This is tidy ugly.  But I feel there 
  ##  is value in keeping dataframe columns as lists...
  df3 = df2 %>%
    select(-Groups) %>%
    group_by(TRT,REP) %>%	
    summarise(Err_mean = list(reduce(map(Err, ~as.matrix(.x)), `+`)/length(Err)),
              Haps_mean = list(reduce(map(Haps, ~as.vector(.x)), `+`)/length(Haps)),
              Names = list(first(Names)),
              Num_mean = sexlink*mean(Num)) %>%
    rename(Haps=Haps_mean,Num=Num_mean,Err=Err_mean)
  
  ## these summaries of the data are pretty useful for tests
  p1 = df3 %>% filter(TRT=="C") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t() 
  row.names(p1) <- NULL
  p2 = df3 %>% filter(TRT=="Z") %>% pull(Haps) %>% as.data.frame() %>% as.matrix() %>% t()
  row.names(p2) <- NULL
  covar1 = do.call(abind, c(df3 %>% filter(TRT=="C") %>% pull(Err), along = 3))
  covar2 = do.call(abind, c(df3 %>% filter(TRT=="Z") %>% pull(Err), along = 3))
  nrepl = df3 %>% filter(TRT=="C") %>% nrow()
  nrepl == df3 %>% filter(TRT=="Z") %>% nrow()
  N1 = df3 %>% filter(TRT=="C") %>% pull(Num)
  N2 = df3 %>% filter(TRT=="Z") %>% pull(Num)
  
  wt = tryCatch(
    wald.test3(p1,p2,covar1,covar2,nrepl,N1,N2),
    error = function(e) list(wald.test=NA, p.value=NA, avg.var=NA)
  )
  Wald_log10p = -log10(wt$p.value)
  Pseu_log10p = tryCatch(
    pseudoN.test(p1,p2,covar1,covar2,nrepl,N1,N2),
    error = function(e) NA
  )
  
  af_cutoff = 0.01     # 1% --- heritability estimators can be off for really low allele frequencies
  temp = Heritability(p1, p2, nrepl, ProportionSelect, af_cutoff)
  Falc_H2 = temp$Falconer_H2
  Cutl_H2 = temp$Cutler_H2
  
  ll = list(Wald_log10p = Wald_log10p, Pseu_log10p = Pseu_log10p,
            Falc_H2 = Falc_H2, Cutl_H2 = Cutl_H2, avg.var = wt$avg.var)
  ll
}


````

</details>


# Outdoor samples CVTK pipeline

## 1. Create CVTK Conda env & Set-up the jupyter notebook env:

<details>
<summary>Click to expand code</summary>

````
conda create -n cvtk python=3.7
conda activate cvtk

git clone https://github.com/vsbuffalo/cvtkpy.git
cd cvtkpy
python setup.py install

````

</details>

## 2. Activate Jupyter Notebook:

<details>
<summary>Click to expand code</summary>

### In terminal

````
pip install jupyter
pip install notebook
jupyter notebook --no-browser

# In VSC; insert for existing kernel: http://localhost:8888/?token=84e2443fdd99d5c85071427b03de889ad49834cbcc11fef0 
# Now you should be connected to the conda env python.
````

</details>

## 3. Install grenedalf & filter 2023 BCF file for E & S samples:

<details>
<summary>Click to expand code</summary>

````
# === Setting up vcf data for .sync output === #

cd /mnt/d/xQTL_2025_Data/Final_Window_Analysis/OutdoorSample_CVTK/input_files

sed -i 's/\r$//' metadata_2023.tsv

bcftools query -l bcftools_filtered-all.vcf.gz > vcf_samples.txt

awk -F'\t' 'NR>1 && $1<=195 && ($8=="E"||$8=="S") && $4>=1 && $4<=4 { key=$6"\t"$4; if(!seen[key]++) cage_tpt[$6]++ } END { for(c in cage_tpt) if(cage_tpt[c]==4) print c }' metadata_2023.tsv | sort -n > complete_cages_full.tx

cat complete_cages_full.txt

awk -F'\t' 'NR==FNR{cages[$1];next} FNR==1{next} $1<=195&&($8=="E"||$8=="S")&&$4>=1&&$4<=4&&$6 in cages&&$12=="spino"{key=$6 SUBSEP $4;if(!seen[key]||($11=="single"&&prev_rep[key]=="dupl")){seen[key]=$1;prev_rep[key]=$11}} $1<=195&&$8=="E"&&$4==0{founders[$1]=1} END{for(k in seen){split(k,a,SUBSEP);cage_count[a[1]]++;samps[a[1],a[2]]=seen[k]} for(c in cage_count){if(cage_count[c]==4){for(t=1;t<=4;t++) print samps[c,t]}} for(f in founders) print f}' complete_cages_full.txt metadata_2023.tsv | sort -n > samples_all.txt

grep -Fxf vcf_samples.txt samples_all.txt > samples_final.txt

wc -l samples_final.txt

awk -F'\t' 'NR==FNR{samps[$1];next} FNR>1&&$1 in samps{print "samp="$1,"tpt="$4,"cage="$6,"treat="$8,"rep="$11}' samples_final.txt metadata_2023.tsv | sort -t= -k3,3n -k2,2n

bcftools view -S samples_final.txt bcftools_filtered-all.vcf.gz -Oz -o filtered_E_S.vcf.gz
bcftools index filtered_E_S.vcf.gz

````


````

# === Grenedalf === #
git clone --recursive https://github.com/lczech/grenedalf.git
cd grenedalf
make

ls grenedalf/grenedalf/bin/
./grenedalf/grenedalf/bin/grenedalf --help

````

</details>

## 4. Fix chromosome data and run grenedalf sync:

<details>
<summary>Click to expand code</summary>

````
echo -e "2L\tchr2L\n2R\tchr2R\n3L\tchr3L\n3R\tchr3R\n4\tchr4\nX\tchrX\nY\tchrY" > chr_rename.tx

bcftools annotate --rename-chrs chr_rename.txt filtered_E_S.vcf.gz -Oz -o filtered_E_S_chr.vcf.g

bcftools index filtered_E_S_chr.vcf.gz

./grenedalf/grenedalf/bin/grenedalf sync --vcf-path filtered_E_S_chr.vcf.gz --reference-genome-fasta /mnt/d/xQTL_2025_Data/ref/dm6.fa --filter-sample-min-read-depth 2 --filter-sample-max-read-depth 1000 --filter-sample-only-biallelic-snps --filter-region [chr2L,chr2R,chr3L,chr3R,chrX] --out-dir sync_file/ --threads 32

# Make sure that sync.sync file does not exist, and the folder is deleted if you got an error running grenedalf and need to rerun it

````


</details>


## 5. Jupyter notebook analysis:

````
!sed -i 's/sklearn/scikit-learn/' /mnt/d/xQTL_2025_Data/Final_Window_Analysis/OutdoorSample_CVTK/cvtkpy/setup.py
!sed -i 's/panda$/pandas/' /mnt/d/xQTL_2025_Data/Final_Window_Analysis/OutdoorSample_CVTK/cvtkpy/setup.py
!pip install /mnt/d/xQTL_2025_Data/Final_Window_Analysis/OutdoorSample_CVTK/cvtkpy/

cd "/mnt/d/xQTL_2025_Data/Final_Window_Analysis/OutdoorSample_CVTK/input_files"

````

bash commands to clean up the sync file
````
gzip /mnt/d/xQTL_2025_Data/Final_Window_Analysis/OutdoorSample_CVTK/input_files/sync_file/sync.sync
zcat /mnt/d/xQTL_2025_Data/Final_Window_Analysis/OutdoorSample_CVTK/input_files/sync_file/sync.sync.gz | tail -n +2 | gzip > /mnt/d/xQTL_2025_Data/Final_Window_Analysis/OutdoorSample_CVTK/input_files/sync_file/sync_noheader.sync.gz
zcat /mnt/d/xQTL_2025_Data/Final_Window_Analysis/OutdoorSample_CVTK/input_files/sync_file/sync_noheader.sync.gz | sed 's/\.\:\.\:\.\:\.\:\.\:\./0:0:0:0:0:0/g' | gzip > /mnt/d/xQTL_2025_Data/Final_Window_Analysis/OutdoorSample_CVTK/input_files/sync_file/sync_clean.sync.gz
````

## Notes for changes to the CVTK pipeline for outdoor samples;

1. We only have 3 founders; CVTK is expecting equal design across timepoints (9 E samples and 5 S samples)
	- So we either sum/average depth values & average frequency values for the SNPs and create 1 compound founder that we align as a reference T0 to all the samples.
    - So I settled on sum the depth values & average the frequency. Reason; given that the 3 founder pools are coming from the same population, fractional coverages generated from averaging make no biological sense. Rounding to whole numbers could be a solution, but seems wrong to do as its direct modification of the data.
  
2. TiledTemporalFreqs function; using the same tile range of 100kBP window size, the only thing changing is the input of data, now seperating data into E and S treatments and pairing all the timepoints together T0 > T1 > T2 > T3.
	- Changing diploids from 1000 to 240 because we have pools of 120 individuals.

WIP RETURN LATER & FINISH

WORK ON 2017 DATA


<img width="793" height="456" alt="image" src="https://github.com/user-attachments/assets/bbbd0a3b-2b0a-4c5b-bffc-9fddac4acdae" />




