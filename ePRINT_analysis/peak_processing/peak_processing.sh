#!/usr/bin/env bash
# =============================================================================
# ePRINT peak processing pipeline for a single experiment (two conditions,
# two biological replicates each).
#
# Execution model: run interactively inside a Singularity bedtools shell.
# Launch the shell first, then source or paste this script.
#
# Pipeline overview (bedtools steps interleaved with R processing):
#   1. Sort & merge per-replicate peak BED files
#   2. Compute IP and input read coverage over merged peaks
#   3. [R step 1] Normalise coverage, filter peaks → <cond>.norm.bed
#   4. Intersect replicates to identify reproducible peaks
#   5. [R step 2] Filter intersected peaks for reproducibility → <cond>.reprod.bed
#   6. Sort & merge reproducible peaks
#   7. [R step 3] Identify common peaks across conditions → common_peaks.unprocessed.bed
#   8. Sort & merge common peaks
#   9. Generate random background regions (bedtools shuffle)
# =============================================================================

# ---------------------------------------------------------------------------
# 0. User-configurable paths and experiment metadata
#    Replace each placeholder with the actual absolute path / value.
# ---------------------------------------------------------------------------

path2sif=/path/to/your/singularity_images/      # Directory containing .sif files
path2bam=/path/to/your/experiment/dedup/        # Deduplicated BAM files
path2bed=/path/to/your/experiment/clipper/      # CLIPper peak BED files
path2out=/path/to/your/experiment/peak_processing/  # All output goes here
path2genomes=/path/to/your/genomes/             # Genome files directory

mkdir -p "${path2out}"

# ---------------------------------------------------------------------------
# Condition 1 — replicate file prefixes and condition labels
# ---------------------------------------------------------------------------
PREFIX_EP_C1_R1=your_ip_replicate1_condition1_prefix   # eCLIP IP, condition 1, rep 1
PREFIX_INP_C1_R1=your_input_replicate1_condition1_prefix
PREFIX_EP_C1_R2=your_ip_replicate2_condition1_prefix
PREFIX_INP_C1_R2=your_input_replicate2_condition1_prefix

COND1_R1=condition1_R1    # Short label used for output file naming
COND1_R2=condition1_R2
COND1=condition1          # Condition-level label (both replicates combined)

# ---------------------------------------------------------------------------
# Condition 2 — replicate file prefixes and condition labels
# ---------------------------------------------------------------------------
PREFIX_EP_C2_R1=your_ip_replicate1_condition2_prefix
PREFIX_INP_C2_R1=your_input_replicate1_condition2_prefix
PREFIX_EP_C2_R2=your_ip_replicate2_condition2_prefix
PREFIX_INP_C2_R2=your_input_replicate2_condition2_prefix

COND2_R1=condition2_R1
COND2_R2=condition2_R2
COND2=condition2

# ---------------------------------------------------------------------------
# Genome / annotation files
# ---------------------------------------------------------------------------
CHROM_SIZES="${path2genomes}/your_genome/chrNameLength.txt"
GENE_BED="${path2genomes}/your_annotation.genes.bed"

# ---------------------------------------------------------------------------
# Read-strand flag used by bedtools coverage.
# Set to "-s" for strand-specific libraries (typical for eCLIP).
# ---------------------------------------------------------------------------
STRAND_FLAG="-s"

# Suffix for CLIPper output files (e.g. _r1, _r2 depending on read orientation)
CLIPPER_SUFFIX=_r1

# Suffix for deduplicated BAM files
BAM_SUFFIX=.dedup.bam


# =============================================================================
# STEP 1 — Sort and merge per-replicate peaks, then compute read coverage
#
# For each replicate:
#   a) Sort the CLIPper peak BED file
#   b) Merge overlapping peaks on the same strand (keeping name, score, strand)
#   c) Count input reads overlapping each merged peak
#   d) Count IP reads overlapping each merged peak (appended as extra columns)
# =============================================================================

for LABEL in "${COND1_R1}" "${COND1_R2}" "${COND2_R1}" "${COND2_R2}"; do

    # Resolve per-replicate prefixes from the label
    case "${LABEL}" in
        "${COND1_R1}") EP_PREFIX="${PREFIX_EP_C1_R1}";  INP_PREFIX="${PREFIX_INP_C1_R1}" ;;
        "${COND1_R2}") EP_PREFIX="${PREFIX_EP_C1_R2}";  INP_PREFIX="${PREFIX_INP_C1_R2}" ;;
        "${COND2_R1}") EP_PREFIX="${PREFIX_EP_C2_R1}";  INP_PREFIX="${PREFIX_INP_C2_R1}" ;;
        "${COND2_R2}") EP_PREFIX="${PREFIX_EP_C2_R2}";  INP_PREFIX="${PREFIX_INP_C2_R2}" ;;
    esac

    echo "[$(date '+%H:%M:%S')] Processing replicate: ${LABEL}"

    # Sort CLIPper peaks
    bedtools sort \
        -i "${path2bed}/${EP_PREFIX}${CLIPPER_SUFFIX}.peakClusters${CLIPPER_SUFFIX}.bed" \
        > "${path2out}/${LABEL}.sorted.bed"

    # Merge overlapping peaks; retain name (col4), min score (col5), strand (col6)
    bedtools merge ${STRAND_FLAG} \
        -c 4,5,6 -o collapse,min,distinct \
        -i "${path2out}/${LABEL}.sorted.bed" \
        > "${path2out}/${LABEL}.merged.bed"

    # Count input (SMInput) reads per merged peak
    bedtools coverage ${STRAND_FLAG} \
        -a "${path2out}/${LABEL}.merged.bed" \
        -b "${path2bam}/${INP_PREFIX}${BAM_SUFFIX}" \
        -split \
        > "${path2out}/${LABEL}.inputcov.bed"

    # Count IP reads per merged peak (input coverage columns already appended)
    bedtools coverage ${STRAND_FLAG} \
        -a "${path2out}/${LABEL}.inputcov.bed" \
        -b "${path2bam}/${EP_PREFIX}${BAM_SUFFIX}" \
        -split \
        > "${path2out}/${LABEL}.cov.bed"

done


# =============================================================================
# [R STEP 1] Normalise read counts and filter low-confidence peaks
#
# Script : peaks_r_treat.R  (section 1)
# Inputs : <replicate>.cov.bed  for each replicate
# Outputs: <replicate>.norm.bed — peaks passing filters with normalised scores
#
# Run this R script now before continuing.
# Example:
#   Rscript /path/to/your/scripts/peaks_r_treat.R \
#       --step norm \
#       --outdir "${path2out}" \
#       --cond1_r1 "${COND1_R1}" --cond1_r2 "${COND1_R2}" \
#       --cond2_r1 "${COND2_R1}" --cond2_r2 "${COND2_R2}"
# =============================================================================


# =============================================================================
# STEP 2 — Intersect replicates to identify reproducible peaks
#
# Peaks from replicate 1 are intersected with peaks from replicate 2.
# Only peaks overlapping on the same strand are kept (-s -wo).
# The resulting file contains paired peak rows for downstream R filtering.
# =============================================================================

for COND LABEL_R1 LABEL_R2 in \
    "${COND1}" "${COND1_R1}" "${COND1_R2}" \
    "${COND2}" "${COND2_R1}" "${COND2_R2}"; do

    echo "[$(date '+%H:%M:%S')] Intersecting replicates for condition: ${COND}"

    bedtools intersect ${STRAND_FLAG} -wo \
        -a "${path2out}/${LABEL_R1}.norm.bed" \
        -b "${path2out}/${LABEL_R2}.norm.bed" \
        > "${path2out}/${COND}.intersect.bed"

done


# =============================================================================
# [R STEP 2] Select reproducible peaks from intersected replicate pairs
#
# Script : peaks_r_treat.R  (section 2)
# Inputs : <cond>.intersect.bed
# Outputs: <cond>.reprod.bed — one row per reproducible peak locus
#
# Run this R script now before continuing.
# =============================================================================


# =============================================================================
# STEP 3 — Sort and merge reproducible peaks within each condition
#
# Columns 4–10 carry scores from both replicates; we keep the maximum value
# for each numeric column and the union of strand labels.
# =============================================================================

for COND in "${COND1}" "${COND2}"; do

    echo "[$(date '+%H:%M:%S')] Merging reproducible peaks for condition: ${COND}"

    bedtools sort ${STRAND_FLAG} \
        -i "${path2out}/${COND}.reprod.bed" \
        > "${path2out}/${COND}.reprod.sorted.bed"

    # Merge; retain max scores (cols 4–5, 7–10) and distinct strand (col 6)
    bedtools merge ${STRAND_FLAG} \
        -c 4,5,6,7,8,9,10 -o max,max,distinct,max,max,max,max \
        -i "${path2out}/${COND}.reprod.sorted.bed" \
        > "${path2out}/${COND}.reprod.merged.bed"

done


# =============================================================================
# [R STEP 3] Identify peaks common to both conditions
#
# Script : peaks_r_treat.R  (section 3)  or  peak_analysis.Rmd
# Inputs : <cond1>.reprod.merged.bed, <cond2>.reprod.merged.bed
# Outputs: common_peaks.unprocessed.bed
#
# Run this R script now before continuing.
# =============================================================================


# =============================================================================
# STEP 4 — Sort and merge common peaks
# =============================================================================

echo "[$(date '+%H:%M:%S')] Sorting and merging common peaks"

bedtools sort ${STRAND_FLAG} \
    -i "${path2out}/common_peaks.unprocessed.bed" \
    > "${path2out}/common_peaks.sorted.bed"

# Retain collapsed peak names, collapsed scores, and distinct strand
bedtools merge ${STRAND_FLAG} \
    -c 4,5,6 -o collapse,collapse,distinct \
    -i "${path2out}/common_peaks.sorted.bed" \
    > "${path2out}/common_peaks.merged.bed"


# =============================================================================
# STEP 5 — Generate random background regions via bedtools shuffle
#
# Shuffled peaks are confined to genic regions (-incl) to provide a
# matched background for motif enrichment analysis.
# N_SHUFFLES sets of randomised coordinates are produced, then concatenated
# into a single random_peaks.bed file.
# =============================================================================

N_SHUFFLES=15   # Total number of shuffled sets; adjust as needed
SHUFFLE_SOURCE="${path2out}/${COND1}.reprod.merged.bed"   # Peak set to shuffle

echo "[$(date '+%H:%M:%S')] Generating ${N_SHUFFLES} shuffled background sets"

# Remove any leftover file from a previous run
rm -f "${path2out}/random_peaks.bed"

for i in $(seq 1 "${N_SHUFFLES}"); do
    echo "  Shuffle iteration ${i}/${N_SHUFFLES}"
    bedtools shuffle \
        -i  "${SHUFFLE_SOURCE}" \
        -g  "${CHROM_SIZES}" \
        -incl "${GENE_BED}" \
        >> "${path2out}/random_peaks.bed"
done

echo "[$(date '+%H:%M:%S')] Pipeline complete."
