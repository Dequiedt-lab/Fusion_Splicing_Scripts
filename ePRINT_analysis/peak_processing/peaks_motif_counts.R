# =============================================================================
# peaks_motif_counts.R
# Cluster-side script: count PWM hits (total, not positional) across all
# peak categories using countPWM for speed.
#
# Complements count_pwm.R: where count_pwm.R gives per-position profiles,
# this script gives total hit counts per peak per motif, enabling statistical
# comparisons between peak categories (exclusive, common, random).
#
# Usage: Rscript peaks_motif_counts.R
# =============================================================================

library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(GenomicRanges)
library(BiocParallel)


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Directory containing peak BED files and FASTA outputs
peak_dir <- "/path/to/your/experiment/peak_processing/"

# Path to the combined PWM RData object (loads 'combined_pwms')
pwm_rdata <- "/path/to/your/resources/motifs/all_combined_pwms.RData"

# Condition labels — must match those used in peak_processing.sh
cond1 <- "condition1"   # e.g. "noDOX"
cond2 <- "condition2"   # e.g. "DOX"

# Extension applied to peaks when extracting sequences for PWM scanning
PEAK_EXT <- 50

# Number of parallel workers
N_WORKERS <- 20

# Minimum PWM match score threshold
PWM_MIN_SCORE <- "80%"


# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

param <- MulticoreParam(workers = N_WORKERS, progressbar = TRUE)
load(pwm_rdata)   # loads: combined_pwms


# -----------------------------------------------------------------------------
# Helper: read a BED file, extend peaks, fetch sequences, count PWM hits
# -----------------------------------------------------------------------------

#' Count PWM hits across all peaks in a BED file
#'
#' @param bed_path  path to a BED file (columns 1-3 = chr/start/end, 6 = strand)
#' @param ext       number of bp to extend each peak on each side
#' @return          matrix: peaks × motifs, each cell = count of PWM hits
getcounts <- function(bed_path, ext = PEAK_EXT) {

  peaks <- read.table(bed_path)
  gr    <- GRanges(
    seqnames = peaks$V1,
    ranges   = IRanges(start = peaks$V2, end = peaks$V3),
    strand   = peaks$V6
  )

  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr + ext)
  names(seqs) <- paste0("peak_", seq_along(seqs))

  # Drop sequences containing ambiguous bases (N) which break countPWM
  has_n <- grep("N", seqs)
  if (length(has_n) > 0) seqs <- seqs[-has_n]

  # Count total PWM hits per sequence for every motif (parallelised over motifs)
  counts_list <- bplapply(combined_pwms, function(pwm) {
    vapply(seqs, function(x) countPWM(pwm, x, min.score = PWM_MIN_SCORE),
           numeric(1))
  }, BPPARAM = param)

  do.call(cbind, counts_list)
}


# -----------------------------------------------------------------------------
# Score all peak categories
# -----------------------------------------------------------------------------

message("Starting motif counting: ", Sys.time())

counts_c1_all       <- getcounts(paste0(peak_dir, cond1, ".reprod.merged.bed"))
message(cond1, " all peaks done: ", Sys.time())

counts_c1_exclusive <- getcounts(paste0(peak_dir, cond1, "_exclusive.bed"))
message(cond1, " exclusive peaks done: ", Sys.time())

counts_c2_all       <- getcounts(paste0(peak_dir, cond2, ".reprod.merged.bed"))
message(cond2, " all peaks done: ", Sys.time())

counts_c2_exclusive <- getcounts(paste0(peak_dir, cond2, "_exclusive.bed"))
message(cond2, " exclusive peaks done: ", Sys.time())

counts_common  <- getcounts(paste0(peak_dir, "common_peaks.merged.bed"))
counts_random  <- getcounts(paste0(peak_dir, "random_peaks.bed"))

message("All peak categories scored: ", Sys.time())


# -----------------------------------------------------------------------------
# Save output
# -----------------------------------------------------------------------------

save(
  counts_c1_all, counts_c1_exclusive,
  counts_c2_all, counts_c2_exclusive,
  counts_common, counts_random,
  file = paste0(peak_dir, "countsPWM_peaks.RData")
)

message("Saved: countsPWM_peaks.RData")
message("All done: ", Sys.time())



