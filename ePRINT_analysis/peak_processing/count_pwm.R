# =============================================================================
# count_pwm.R
# Cluster-side script: score each peak sequence against all PWMs using
# position-weight matrix scanning, then compute per-position mean counts.
#
# Designed for execution on an HPC node with many cores available (SLURM submission).
# Inputs are FASTA files produced by peak_analysis.Rmd (Section 5).
# Output is an RData file loaded by peak_analysis.Rmd (Section 9).
#
# Usage: Rscript count_pwm.R
# =============================================================================

library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(universalmotif)
library(BiocParallel)


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Directory containing peak FASTA files (output of peak_analysis.Rmd)
peak_dir <- "/path/to/your/experiment/peak_processing/"

# Path to the combined PWM RData object (loads 'combined_pwms')
pwm_rdata <- "/path/to/your/resources/motifs/all_combined_pwms.RData"

# Condition labels — must match those used in peak_analysis.Rmd
cond1 <- "condition1"   # e.g. "noDOX"
cond2 <- "condition2"   # e.g. "DOX"

# Extension used when writing FASTA files in peak_analysis.Rmd
PEAK_EXT <- 50

# Minimum PWM match score threshold
PWM_MIN_SCORE <- "80%"

# Number of parallel workers (set to the number of cores allocated to the job)
N_WORKERS <- 20

# Fixed window length assumed when computing per-position counts.
# Must match the actual length of the (extended) peak sequences.
WINDOW_LEN <- 100


# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

param <- MulticoreParam(workers = N_WORKERS, progressbar = TRUE)
load(pwm_rdata)   # loads: combined_pwms

message("Starting PWM scoring: ", Sys.time())


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

#' Count PWM hits at each position across a set of sequences
#'
#' For each sequence, matchPWM is run at `PWM_MIN_SCORE`. The resulting match
#' positions are projected onto a 1:WINDOW_LEN position vector and summed,
#' giving a WINDOW_LEN × n_sequences matrix of per-position hit counts.
#'
#' @param pwm    a PWM object (from combined_pwms)
#' @param seqs   DNAStringSet
#' @return       matrix: sequences × positions
count_pwm <- function(pwm, seqs) {
  hits_per_seq <- lapply(seqs, function(x) {
    matches <- matchPWM(pwm, x, min.score = PWM_MIN_SCORE)
    countOverlaps(IRanges(seq_len(WINDOW_LEN), width = 1), IRanges(matches))
  })
  t(do.call(rbind, hits_per_seq))
}


# -----------------------------------------------------------------------------
# Load FASTA files
# -----------------------------------------------------------------------------

read_fasta <- function(label, ext = PEAK_EXT, dir = peak_dir) {
  readDNAStringSet(paste0(dir, label, "_seqs_ext", ext, ".fasta"))
}

seqs_c1       <- read_fasta(cond1)
seqs_c2       <- read_fasta(cond2)

# Randomised controls:
#   _shuff  = dinucleotide-shuffled versions of the actual peak sequences
#             (sequence composition matched, positional signal destroyed)
#   _rand   = sequences from genomically shuffled peak coordinates
#             (bedtools shuffle output, first shuffled replicate used)
seqs_c1_shuff <- shuffle_sequences(seqs_c1)
seqs_c2_shuff <- shuffle_sequences(seqs_c2)
seqs_c1_rand  <- read_fasta(paste0(cond1, "shuf1"))
seqs_c2_rand  <- read_fasta(paste0(cond2, "shuf1"))

message(cond1, " peaks: ", length(seqs_c1), " sequences")
message(cond2, " peaks: ", length(seqs_c2), " sequences")


# -----------------------------------------------------------------------------
# PWM scoring — parallelised over motifs
# -----------------------------------------------------------------------------

score_all <- function(seqs) {
  bplapply(combined_pwms, count_pwm, seqs, BPPARAM = param)
}

message("Scoring ", cond1, " peaks ...")
pwm_counts_c1       <- score_all(seqs_c1)
pwm_counts_c1_shuff <- score_all(seqs_c1_shuff)
pwm_counts_c1_rand  <- score_all(seqs_c1_rand)

message("Scoring ", cond2, " peaks ...")
pwm_counts_c2       <- score_all(seqs_c2)
pwm_counts_c2_shuff <- score_all(seqs_c2_shuff)
pwm_counts_c2_rand  <- score_all(seqs_c2_rand)


# -----------------------------------------------------------------------------
# Summarise: mean hit count per position across all peaks (used for heatmaps)
# -----------------------------------------------------------------------------

pwm_means_c1       <- sapply(pwm_counts_c1,       colMeans)
pwm_means_c1_shuff <- sapply(pwm_counts_c1_shuff, colMeans)
pwm_means_c1_rand  <- sapply(pwm_counts_c1_rand,  colMeans)
pwm_means_c2       <- sapply(pwm_counts_c2,       colMeans)
pwm_means_c2_shuff <- sapply(pwm_counts_c2_shuff, colMeans)
pwm_means_c2_rand  <- sapply(pwm_counts_c2_rand,  colMeans)

# Rename objects to use condition labels so peak_analysis.Rmd can retrieve them
# with get(paste0("pwm_means_", cond1, "ext50")) etc.
rename_obj <- function(obj, new_name) {
  assign(new_name, obj, envir = .GlobalEnv)
}

rename_obj(pwm_means_c1,       paste0("pwm_means_", cond1, "ext", PEAK_EXT))
rename_obj(pwm_means_c1_shuff, paste0("pwm_means_", cond1, "ext", PEAK_EXT, "_shuff"))
rename_obj(pwm_means_c1_rand,  paste0("pwm_means_", cond1, "ext", PEAK_EXT, "_rand"))
rename_obj(pwm_means_c2,       paste0("pwm_means_", cond2, "ext", PEAK_EXT))
rename_obj(pwm_means_c2_shuff, paste0("pwm_means_", cond2, "ext", PEAK_EXT, "_shuff"))
rename_obj(pwm_means_c2_rand,  paste0("pwm_means_", cond2, "ext", PEAK_EXT, "_rand"))


# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------

save(
  list = c(
    paste0("pwm_means_", cond1, "ext", PEAK_EXT),
    paste0("pwm_means_", cond1, "ext", PEAK_EXT, "_shuff"),
    paste0("pwm_means_", cond1, "ext", PEAK_EXT, "_rand"),
    paste0("pwm_means_", cond2, "ext", PEAK_EXT),
    paste0("pwm_means_", cond2, "ext", PEAK_EXT, "_shuff"),
    paste0("pwm_means_", cond2, "ext", PEAK_EXT, "_rand")
  ),
  file = paste0(peak_dir, "pwm_colmeans_ext", PEAK_EXT, ".RData")
)

message("Saved: pwm_colmeans_ext", PEAK_EXT, ".RData")
message("Done: ", Sys.time())



