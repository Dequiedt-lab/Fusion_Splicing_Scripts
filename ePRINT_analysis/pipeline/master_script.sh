#!/bin/bash
#
#SBATCH --job-name=master_eprint
#SBATCH --partition=all_5hrs
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --output=/path/to/your/eprint/Logs/master_script.%j.log

printf "#########################################\n\nStarting master script\n\n#########################################\n\n"

# ============================================================
# Per-sample variables — update these for each new sample
# ============================================================

# Directory containing raw FASTQ files for this sample
export path2fastq=/path/to/your/rawdata/fastq/

# Sample prefix: base filename shared by R1 and R2 FASTQ files
# Expected file names: ${prefix}_R1_001.fastq.gz and ${prefix}_R2_001.fastq.gz
export prefix=SAMPLENAME_S1_L001

# ============================================================
# Constant variables — update once for your environment
# ============================================================

export path2wd=/path/to/your/eprint/run_directory
export path2genomes=/path/to/your/genomes
export path2sif=/path/to/your/SIF_files
export path2scripts=/path/to/your/eprint/scripts   # directory containing umi2dedup.sh etc.

# ============================================================
# Step 1: UMI extraction → trimming → alignment → deduplication
# ============================================================

jid1=$(sbatch --export=ALL ${path2scripts}/umi2dedup.sh | cut -d ' ' -f4)
printf "Launching job number %s for umi2dedup\n\n\n\n" "$jid1"

# ============================================================
# Step 2: Quality controls (FastQC, PreSeq, RSeQC)
# ============================================================

jid2=$(sbatch --export=ALL --dependency=afterok:$jid1 ${path2scripts}/qc_umi2dedup.sh | cut -d ' ' -f4)
printf "Launching job number %s for QCs\n\n\n\n" "$jid2"

# ============================================================
# Step 3: Peak calling with CLIPper
# ============================================================

jid3=$(sbatch --export=ALL --dependency=afterok:$jid2 ${path2scripts}/clipper_eprint2.sh | cut -d ' ' -f4)
printf "Launching job number %s for CLIPper\n\n\n\n" "$jid3"

