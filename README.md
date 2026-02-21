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

## 6. Create sample table info with bams
````
#!/bin/bash
# Usage: bash make_bam_list.sh /path/to/output_list.txt

out_list=$1

if [ -z "$out_list" ]; then
    echo "Usage: $0 <output_list.txt>"
    exit 1
fi

# Define your two BAM directories
dir1="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/Spino3/RG_groups"
dir2="/mnt/d/grenepipe/gp_analysis/rudflies_2024/dedup"

> "$out_list"

for bam_dir in "$dir1" "$dir2"; do
    if [ ! -d "$bam_dir" ]; then
        echo "Warning: directory $bam_dir does not exist, skipping."
        continue
    fi
    for bam in "$bam_dir"/*.bam; do
        [ -e "$bam" ] || continue
        realpath "$bam" >> "$out_list"
    done
    echo "Added BAMs from $bam_dir"
done

echo "BAM list written to $out_list"
````

## 7. Make Bam table:
````
#!/usr/bin/env bash
# Usage: ./make_bam_table.sh output.tsv

OUT_TSV=${1:-bam_list.tsv}

# Define your two BAM directories
dir1="/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/Spino3/RG_groups"
dir2="/mnt/d/grenepipe/gp_analysis/rudflies_2024/dedup"

# Header
echo -e "sample\tbam" > "$OUT_TSV"

for BAM_DIR in "$dir1" "$dir2"; do
    if [ ! -d "$BAM_DIR" ]; then
        echo "Warning: directory $BAM_DIR does not exist, skipping."
        continue
    fi
    find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" | sort | while read -r bamfile; do
        sample=$(basename "$bamfile" .bam)
        echo -e "${sample}\t${bamfile}" >> "$OUT_TSV"
    done
    echo "Added BAMs from $BAM_DIR"
done

echo "✅ TSV file created: $OUT_TSV"
````

## 8. Combine Founders & Samples:
````
#!/bin/bash
# Usage: bash combined_local_simple.sh bam_list.txt output_dir

ref="/mnt/d/xQTL_2025_Data/ref/dm6.fa"
bams=$1
output=$2
mkdir -p "$output"

# --- Check and index BAMs if needed ---
echo "=== Checking BAM indexes ==="
while read -r bam; do
    if [[ ! -f "${bam}.bai" ]]; then
        echo "Indexing: $bam"
        samtools index "$bam"
    else
        echo "Already indexed: $(basename $bam)"
    fi
done < "$bams"
echo "=== Indexing complete ==="

# --- Run mpileup per chromosome in parallel ---
declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")

run_chrom() {
    mychr=$1
    echo "Processing chromosome $mychr"

    bcf_out="${output}/calls.${mychr}.bcf"

    bcftools mpileup -I -d 1000 -r "$mychr" -a "FORMAT/AD,FORMAT/DP" \
        -f "$ref" -b "$bams" --threads 12 | \
        bcftools call -mv --threads 12 -Ob -o "$bcf_out"

    bcftools index "$bcf_out"

    echo -ne "CHROM\tPOS" > "${output}/RefAlt.${mychr}.txt"
    bcftools query -l "$bcf_out" | awk '{printf("\tREF_%s\tALT_%s",$1,$1)}' >> "${output}/RefAlt.${mychr}.txt"
    echo -ne "\n" >> "${output}/RefAlt.${mychr}.txt"

    bcftools view -m2 -M2 -v snps -i 'QUAL>59' "$bcf_out" |
        bcftools query -f '%CHROM %POS [ %AD{0} %AD{1}] [%GT]\n' |
        grep -v '\.' | awk 'NF-=1' >> "${output}/RefAlt.${mychr}.txt"

    echo "Finished chromosome $mychr"
}

export -f run_chrom
export ref bams output

# Run all chromosomes in parallel
for mychr in "${chrs[@]}"; do
    run_chrom "$mychr" &
done
wait

echo "=== All chromosomes complete ==="
````
## 8.5. Do H cut off assessment (hierarchical clustering): <<< DGRP ONLY >>>
````
library(tidyverse)

# --- inputs ---
setwd("/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/Spino3/RG_groups/haplotype_inf")
filein <- "/mnt/d/xQTL_2025_Data/Final_Window_Analysis/DGRP_xQTL/Spino3/RG_groups/haplotype_inf/RefAlt.chrX.txt"
source("haplotype.parameters.r")  # loads founders, size, h_cutoff

# --- read just what we need ---
df <- read.table("RefAlt.chr2L.txt", header=TRUE)

df2 <- df %>%
  pivot_longer(c(-CHROM,-POS), names_to="lab", values_to="count") %>%
  mutate(RefAlt = str_sub(lab,1,3),
         name   = str_sub(lab,5)) %>%
  select(-lab) %>%
  pivot_wider(names_from=RefAlt, values_from=count) %>%
  mutate(freq = REF/(REF+ALT), N = REF+ALT) %>%
  select(-c("REF","ALT"))

# --- pick median window ---
mid <- median(df2$POS[df2$CHROM=="chr2L"])
test_window <- df2 %>%
  filter(CHROM=="chr2L" & 
           POS > (mid - size/2) & 
           POS < (mid + size/2) & 
           name %in% founders) %>%
  select(-c(CHROM,N)) %>%
  pivot_wider(names_from=name, values_from=freq)

m_test <- as.matrix(test_window %>% select(-POS))

# replace NAs with 0
m_test[is.na(m_test)] <- 0

d <- dist(t(m_test))

# --- diagnostics ---
cat("Distance summary:\n")
print(summary(as.numeric(d)))
cat("\nCurrent h_cutoff:", h_cutoff, "\n")

# test a range of cutoffs to find one that gives ~8-20 groups
for(h in c(2.5, 5, 10, 20, 30, 40, 50)){
  n_groups <- max(cutree(hclust(d), h=h))
  cat("h_cutoff =", h, "->", n_groups, "groups\n")
}

for(h in seq(7.0, 9.5, by=0.1)){
  n_groups <- max(cutree(hclust(d), h=h))
  cat("h_cutoff =", h, "->", n_groups, "groups\n")
}

for(test_pos in c(2000000, 6000000, 10000000, 14000000, 18000000)){
  tw <- df2 %>%
    filter(CHROM=="chr2L" & POS > (test_pos - size/2) & POS < (test_pos + size/2) & name %in% founders) %>%
    select(-c(CHROM,N)) %>%
    pivot_wider(names_from=name, values_from=freq)
  m <- as.matrix(tw %>% select(-POS))
  m[is.na(m)] <- 0
  d <- dist(t(m))
  cat("pos", test_pos, "->", max(cutree(hclust(d), h=7.5)), "groups\n")
}

# visualize
hist(as.numeric(d), breaks=100,
     main="Pairwise founder distances in test window",
     xlab="Euclidean distance")
abline(v=h_cutoff, col="red", lwd=2, lty=2)

````
<img width="860" height="416" alt="image" src="https://github.com/user-attachments/assets/c3227001-80fa-4859-9f2f-35d108896eb8" />

Because we have ~226 DGRP founders, we need to make sure that hierachical clustering works properly and creates a proper amount of haplotype blocks (x>1).

Number of haplotypes has to stay biologically relevant, while not being overdiscriminatory (>100)

## 9. Modify your input_table and haplotype.parameters files:
haplotype.parameters
````
# JUICE Haplotype Parameters
# This file defines the parameters for haplotype estimation in the JUICE dataset

# Founder samples (these should match the column names in your REFALT data)
founders <- c("196","197","198","199","200","201","202","203","204","205","206","207","208","209","210","211","212","213","214","215","216","217","218","219","220","221","222","223","224","225","226","227","228","229","230","231","232","233","234","235","236","237","238","239","240","241","242","243","244","245","246","247","248","249","250","251","252","253","254","255","256","257","258","259","260","261","262","263","264","265","266","267","268","269","270","271","272","273","274","275","276","277","278","279","280","281","282","283","284","285","286","287","288","289","290","291","292","293","294","295","296","297","298","299","300","301","302","303","304","305","306","307","308","309","310","311","312","313","314","315","316","317","318","319","320","321","322","323","324","325","326","327","328","329","330","331","332","333","334","335","336","337","338","339","340","341","342","343","344","345","346","347","348","349","350","351","352","353","354","355","356","357","358","359","360","361","362","363","364","365","366","368","369","370","371","372","373","374","375","376","377","378","379","380","381","382","383","384","385","386","387","388","389","390","391","392","393","394","395","396","397","398","399","400","401","402","404","405","406","407","408","409","410","411","412","413","414","415","416","417","418","419","420","421","422","423")

# Step size for scanning positions along chromosome
step <- 1000  # 10kb steps (increased from 5kb for speed)

#try 1kb, 50kb, 100kb, 200kb

#maybe 10kb

# Window size - COMMENTED OUT, comes from command line
size <- 50000  # Window size comes from command line parameter

#try 100kb, 200kb, 300kb, 400kb, 500kb

# Clustering parameters
h_cutoff <- 7.5  # Default h_cutoff for fixed window hierarchical clustering , 2.5 for DSPR

# Samples to process (allows selecting subset from large REFALT files)
names_in_bam=c("C1","C2","C3","C4","C5","C6","C7","C8","A1","A2","A3","A4","A5","A6","A7","A8")

````

input_table
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

## 10. Run REFALT2haps.Andreas.sh (for DGRP, running a modified REFALT2haps.Andreas.code.r, reference below!)
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

Founders data: really low coverage ~11x which creates a multitude of problems in REFALT2haps step.
1. Code is expecting that all founders have at least 1 read for that specific SNP. When you have 8 founders that is a non-issue, but when its 226, 99.6% of the sites get filtered out.

2. Code is expecting homozygosity at every site, ALL founder samples have to be <0.03 or >0.97. Which is problematic in founder samples that are low coverage (DGRP RILs are ~11x coverage). Which leads to a lot of samples actually being non-homozygotic, i.e. somewhere between 0.03 and 0.97. Which gets filtered out (the rest up to 99.9% gets filtered here)
<img width="852" height="414" alt="image" src="https://github.com/user-attachments/assets/2851764b-e950-47bc-962c-54f9ad737f08" />

3. Some sites are monomorphic between ALL founder samples (AF = 1, or AF = 0), which is bad for haplotype cluster estimation, how do you determine which haplotype it belongs to? So that get's filtered out.

In the end, a ~500,000 SNP chr2L ends up with ~1621 SNPs across the entire genome post filtering.


### Late Night Solution (add filtering that imitates hafpipe)
````
df3 <- df2 %>%
  mutate(freq = case_when(
    name %in% founders & freq > 0.55 ~ 1,
    name %in% founders & freq < 0.45 ~ 0,
    name %in% founders & freq >= 0.45 & freq <= 0.55 ~ NA_real_,
    TRUE ~ freq
  )) %>%
  filter(!(name %in% founders & is.na(freq)))
````



## 11. Run haps2scan.Apr2025
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

## Modified DGRP Scripts:
REFALT2haps.Andreas.code.r <<< Modified est_hap2 function to work with NAs in the data, and modified df2 filtering, no longer performs zero, notfixed, and informative filtering. Plus added a progress bar during window analysis.
````
nrow_subset = function(spotsdf, df3){
        df3 %>%
                filter(CHROM==spotsdf$CHROM &
                        POS > spotsdf$start &
                        POS < spotsdf$end &
                        (name %in% founders | name %in% names_in_bam)) %>%
                nrow()
        }




est_hap = function(spotsdf, df3){
        # spotsdf = spots$data[[1]]  testing
        temp_mat = df3 %>%
                filter(CHROM==spotsdf$CHROM &
                        POS > spotsdf$start &
                        POS < spotsdf$end &
                        (name %in% founders | name %in% names_in_bam)) %>%
                select(-c(CHROM,N)) %>%
                pivot_wider(names_from=name, values_from=freq) %>%
                pivot_longer(!c("POS",matches(founders)),names_to = "sample", values_to = "freq") %>%
                select(-POS)

        sample_mat = temp_mat %>%
                group_by(sample) %>%
                nest() %>%
                mutate(haps=map(data,est_hap2)) %>%
                select(-data) %>%
                unnest_wider(haps)

        sample_mat
        }
        
# function for estimating the haplotype for a given sample, called from est_hap
# returns Groups from cuttree, hap freq estimates, errors by founder
est_hap2 = function(sampdf){
  founder_mat = sampdf %>% select(matches(founders))
  Y = sampdf$freq
  
  # Original: only filters NA in Y (pool)
  # Fix: also filter rows where any founder is NA
  good = !is.na(Y) & complete.cases(founder_mat)
  
  Y = Y[good]
  founder_mat = founder_mat[good,] 
  m_founder_mat = as.matrix(founder_mat)
  
  # Guard against empty matrix after filtering
  if(nrow(m_founder_mat) < 10 | ncol(m_founder_mat) < 2){
    return(list(Groups=NA, Haps=NA, Err=NA, Names=NA))
  }
  
  Groups = cutree(hclust(dist(t(m_founder_mat))), h=h_cutoff)
  d = ncol(m_founder_mat)         
  out = lsei(A=m_founder_mat, B=Y, E=t(matrix(rep(1,d))), F=1,
             G=diag(d), H=matrix(rep(0.0003,d)), verbose=TRUE, fulloutput=TRUE)
  Haps = out$X
  Err = out$cov
  list(Groups=Groups, Haps=Haps, Err=Err, Names=names(Haps))
}

df = lazy_dt(read.table(filein,header=TRUE))
df2 = df %>%
	pivot_longer(c(-CHROM,-POS), names_to = "lab", values_to = "count") %>%
	mutate(RefAlt = str_sub(lab,1,3)) %>%
	mutate(name = str_sub(lab,5)) %>%
	select(-lab) %>%
#	separate(lab, c("RefAlt", "name"), "_", extra = "merge") %>%
	pivot_wider(names_from = RefAlt, values_from = count) %>%
	mutate(freq = REF/(REF+ALT), N = REF+ALT) %>%
	select(-c("REF","ALT")) %>%
	as_tibble()

rm(df)
cat("df2 is now made\n")

# --- debug block ---
#good_SNPs_debug = df2 %>%
#  filter(name %in% founders) %>%
#  group_by(CHROM, POS) %>%
#  summarize(
#    zeros       = sum(N==0),
#    notfixed    = sum(N!=0 & freq > 0.03 & freq < 0.97),
#    informative = (sum(freq, na.rm=TRUE) > 0.05 | sum(freq, na.rm=TRUE) < 0.95)
#  ) %>%
#  ungroup()
#cat("Total SNPs before filtering:", nrow(good_SNPs_debug), "\n")
#cat("Fail zeros==0:", sum(good_SNPs_debug$zeros > 0), "\n")
#cat("Fail notfixed==0:", sum(good_SNPs_debug$notfixed > 0), "\n")
#cat("Fail informative:", sum(!good_SNPs_debug$informative), "\n")
#cat("Pass all filters:", sum(good_SNPs_debug$zeros==0 & good_SNPs_debug$notfixed==0 & good_SNPs_debug$informative), "\n")
#rm(good_SNPs_debug)
# --- end debug ---

# identify SNPs that are NOT problematic in the set of founders
#good_SNPs = df2 %>%
#  filter(name %in% founders) %>%
#  group_by(CHROM,POS) %>%
#  summarize(
#    zeros      = sum(N==0),
#    notfixed   = sum(N!=0 & freq > 0.03 & freq < 0.97),
#    informative = (sum(freq)>0.05 | sum(freq) < 0.95)
#  ) %>%
#  ungroup() %>%
#  filter(zeros < 20 & notfixed < 5 & informative=="TRUE") %>%
#  select(c(CHROM,POS))

# now subset the entire dataset for the good SNPs only
#df3 = good_SNPs %>% left_join(df2, multiple = "all")
df3 <- df2

rm(df2)  # <-- df2 deleted here, after debug block is done
cat("df3 is now made (no SNP filtering)\n")

# --- h_cutoff diagnostic: run once then comment out ---
# pick midpoint of the chromosome's data range for the test window
#mid <- median(df3$POS[df3$CHROM == mychr])
#test_window <- df3 %>%
#  filter(CHROM == mychr & POS > (mid - size/2) & POS < (mid + size/2) & name %in% founders) %>%
#  select(-c(CHROM, N)) %>%
#  pivot_wider(names_from = name, values_from = freq)

#m_test <- as.matrix(test_window %>% select(-POS))
#m_test[is.na(m_test)] <- 0  # ADD THIS LINE
#d <- dist(t(m_test))
#cat("Distance summary:\n")
#print(summary(as.numeric(d)))
#cat("Current h_cutoff:", h_cutoff, "\n")
#cat("Founders collapsed at h_cutoff:", 
#    max(cutree(hclust(d), h=h_cutoff)), "groups from", length(founders), "founders\n")

# --- end diagnostic ---
saveRDS(df3, file = rdsfile)

# df3 = readRDS(rdsfile)
# spots are the locations at which we will estimate haplotypes
# every <step> bp (i.e., 10kb) on the step (i.e, 0, 10, 20, ... kb)
# I define a window +/- size on those steps, and fix the ends
minpos = min(df3$POS)
maxpos = max(df3$POS)
myseq = seq(0,maxpos,step)
myseq = myseq[myseq > minpos + size & myseq < maxpos - size]
spots = data.frame(CHROM=rep(mychr,length(myseq)), pos=myseq, start=myseq-size, end=myseq+size)
# get rid of windows with fewer than 50 SNPs
# i.e., <50 SNPs in 100kb is pretty strange
UU = unique(df3$POS)
spots = spots %>%
  rowwise() %>%
  mutate(NN = sum(start < UU) - sum(end < UU)) %>%
  filter(NN >= 50) %>%
  select(-NN)

# this is the actual scan
#library(furrr)
#options(future.globals.maxSize = 8 * 1024^3)  # 8GB limit
#plan(multisession, workers=4)

total_windows <- nrow(spots)
cat(sprintf("Starting scan: %d windows to process\n", total_windows))

spots2 = spots %>%
  group_nest(row_number()) %>%
  mutate(out = imap(data, function(x, i) {
    if(i %% 100 == 0) cat(sprintf("Progress: %d/%d windows (%.1f%%)\n",
                                  i, total_windows, 100*i/total_windows))
    flush.console()
    est_hap(x, df3)
  })) %>%
  unnest(data) %>%
  select(-c(start,end)) %>%
  select(-`row_number()`) %>%
  unnest_wider(out)

cat("Scan complete\n")
saveRDS(spots2, file=fileout)


````




