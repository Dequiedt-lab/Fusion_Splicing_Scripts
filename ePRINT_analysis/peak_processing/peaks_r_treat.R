# =============================================================================
# peaks_r_treat.R
# ePRINT peak processing — R-side normalisation and reproducibility filtering.
#
# Called at two points interleaved with peak_processing.sh:
#
#   [R STEP 1]  Per-replicate coverage normalisation + enrichment test
#               Inputs:  <replicate>.cov.bed   (bedtools coverage output)
#               Outputs: <replicate>.norm.bed  → consumed by bedtools intersect
#
#   [R STEP 2]  Reproducibility filtering of intersected replicate pairs
#               Inputs:  <condition>.intersect.bed  (bedtools intersect output)
#               Outputs: <condition>.reprod.bed     → consumed by bedtools merge
#
# Usage: Rscript peaks_r_treat.R
# =============================================================================


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Directory containing all intermediate BED files
wd <- "/path/to/your/experiment/peak_processing/"

# Condition and replicate labels — must match labels used in peak_processing.sh
cond1    <- "condition1"     # e.g. "noDOX"
cond2    <- "condition2"     # e.g. "DOX"
cond1_r1 <- "condition1_R1"
cond1_r2 <- "condition1_R2"
cond2_r1 <- "condition2_R1"
cond2_r2 <- "condition2_R2"

# Total mapped read counts for each library.
# Obtained from samtools (see 
reads_ep_c1_r1  <- 0   # IP read count, condition 1, replicate 1
reads_ep_c1_r2  <- 0   # IP read count, condition 1, replicate 2
reads_ep_c2_r1  <- 0   # IP read count, condition 2, replicate 1
reads_ep_c2_r2  <- 0   # IP read count, condition 2, replicate 2

reads_inp_c1_r1 <- 0   # SMInput read count, condition 1, replicate 1
reads_inp_c1_r2 <- 0   # SMInput read count, condition 1, replicate 2
reads_inp_c2_r1 <- 0   # SMInput read count, condition 2, replicate 1
reads_inp_c2_r2 <- 0   # SMInput read count, condition 2, replicate 2

# Reproducibility thresholds (R STEP 2)
MIN_EPRINT_READS <- 10    # Minimum raw IP reads required in both replicates
MIN_FC           <- 1.5   # Minimum fold-change required in both replicates

# Replacement for Inf values in -log10(p) columns (when binom.test gives p = 0)
INF_SENTINEL <- 400


# =============================================================================
# [R STEP 1]  Per-replicate normalisation and enrichment testing
#
# For each merged peak in a replicate:
#   - Pseudocount of 1 added to both input and IP read counts.
#   - Counts converted to RPM using the respective library sizes.
#   - Fold-change = RPM(IP) / RPM(input).
#   - One-sided binomial test of IP enrichment over input.
#
# Output columns written to <replicate>.norm.bed:
#   chr | start | end | id | pval.peak | strand |
#   input.read | eprint.read | norm.input.read | norm.eprint.read |
#   FC | pval | log2FC | logPval
# =============================================================================

#' Normalise coverage and compute per-peak enrichment statistics
#'
#' @param dt       data.frame from a .cov.bed file (bedtools coverage output)
#' @param r_input  total mapped reads in the SMInput library
#' @param r_ep     total mapped reads in the IP library
#' @return         data.frame with normalised scores and enrichment p-values
test_fc <- function(dt, r_input, r_ep) {

  # Drop the fractional-coverage columns produced by bedtools coverage -split:
  # cols 8-10 = bases covered / feature length / fraction for input
  # cols 12-14 = the same trio for IP
  t1 <- dt[, -c(8, 9, 10, 12, 13, 14)]

  # Pseudocount
  t1$V7  <- t1$V7  + 1   # input
  t1$V11 <- t1$V11 + 1   # IP

  # RPM normalisation
  t1$V8  <- (t1$V7  / r_input) * 1e6
  t1$V12 <- (t1$V11 / r_ep)    * 1e6

  # Fold-change IP / input (in RPM space)
  t1$fc <- t1$V12 / t1$V8

  # One-sided binomial test: probability that IP reads are more frequent than
  # expected if IP and input were drawn from the same pool proportionally
  null_prob <- r_ep / (r_input + r_ep)
  p_vals <- vapply(seq_len(nrow(t1)), function(i) {
    binom.test(
      x           = t1$V11[i],
      n           = t1$V7[i] + t1$V11[i],
      p           = null_prob,
      alternative = "greater"
    )$p.value
  }, numeric(1))
  t1$pval <- p_vals

  colnames(t1) <- c(
    "chr", "start", "end", "id", "pval.peak", "strand",
    "input.read", "eprint.read",
    "norm.input.read", "norm.eprint.read",
    "FC", "pval"
  )

  t1$log2FC  <- log2(t1$FC)
  t1$logPval <- -log10(t1$pval)
  t1$logPval[is.infinite(t1$logPval)] <- INF_SENTINEL

  return(t1)
}

# Build a parameter list for all four replicates, then iterate
replicates <- list(
  list(label = cond1_r1, r_input = reads_inp_c1_r1, r_ep = reads_ep_c1_r1),
  list(label = cond1_r2, r_input = reads_inp_c1_r2, r_ep = reads_ep_c1_r2),
  list(label = cond2_r1, r_input = reads_inp_c2_r1, r_ep = reads_ep_c2_r1),
  list(label = cond2_r2, r_input = reads_inp_c2_r2, r_ep = reads_ep_c2_r2)
)

norm_results <- list()

for (rep in replicates) {
  message("Normalising replicate: ", rep$label)

  dt <- read.delim(paste0(wd, rep$label, ".cov.bed"), header = FALSE)

  result <- test_fc(dt, r_input = rep$r_input, r_ep = rep$r_ep)

  write.table(
    result,
    file      = paste0(wd, rep$label, ".norm.bed"),
    quote     = FALSE, sep = "\t",
    row.names = FALSE, col.names = FALSE
  )

  norm_results[[rep$label]] <- result
  message("  Written: ", rep$label, ".norm.bed  (", nrow(result), " peaks)")
}

# Retain column names for the intersect-parsing step below
norm_colnames <- colnames(norm_results[[1]])

message("\n[R STEP 1 complete] — run the bedtools intersect steps in the shell now.\n")


# =============================================================================
# [R STEP 2]  Reproducibility filtering
#
# The bedtools intersect -wo output concatenates R1 norm.bed columns and
# R2 norm.bed columns side by side, with a final overlap-width column.
#
# A peak pair is kept when BOTH replicates satisfy:
#   raw IP reads > MIN_EPRINT_READS  AND  fold-change > MIN_FC
#
# The reproducible locus spans the union of both replicate coordinates.
#
# Output columns written to <condition>.reprod.bed:
#   chr | start | end | peak.logPval.R1 | peak.logPval.R2 | strand |
#   log2FC.R1 | logPval.enrich.R1 | log2FC.R2 | logPval.enrich.R2
# =============================================================================

#' Filter an intersect BED for reproducible peaks
#'
#' @param intersect_file  path to <condition>.intersect.bed
#' @param col_names       column names from a norm.bed data.frame
#' @return                data.frame of reproducible peaks ready to write
filter_reproducible <- function(intersect_file, col_names) {

  int <- read.delim(intersect_file, header = FALSE)

  colnames(int) <- c(
    paste0(col_names, ".R1"),
    paste0(col_names, ".R2"),
    "overlap_width"
  )

  # Reproducibility filter
  keep <- with(int,
    eprint.read.R1 > MIN_EPRINT_READS &
    eprint.read.R2 > MIN_EPRINT_READS &
    FC.R1          > MIN_FC            &
    FC.R2          > MIN_FC
  )
  int2 <- int[keep, ]

  out <- data.frame(
    chr               = int2$chr.R1,
    start             = pmin(int2$start.R1, int2$start.R2),   # union start
    end               = pmax(int2$end.R1,   int2$end.R2),     # union end
    peak.id.R1        = int2$id.R1,
    peak.id.R2        = int2$id.R2,
    strand            = int2$strand.R1,
    peak.logPval.R1   = -log10(int2$pval.peak.R1),
    peak.logPval.R2   = -log10(int2$pval.peak.R2),
    log2FC.R1         = int2$log2FC.R1,
    logPval.enrich.R1 = int2$logPval.R1,
    log2FC.R2         = int2$log2FC.R2,
    logPval.enrich.R2 = int2$logPval.R2
  )

  out[out == Inf] <- INF_SENTINEL

  # Columns expected by downstream bedtools merge (cols 4-10 used as -c arguments)
  return(out[, c(1, 2, 3, 7, 8, 6, 9, 10, 11, 12)])
}

for (cond in c(cond1, cond2)) {
  message("Filtering reproducible peaks for condition: ", cond)

  reprod <- filter_reproducible(
    intersect_file = paste0(wd, cond, ".intersect.bed"),
    col_names      = norm_colnames
  )

  write.table(
    reprod,
    file      = paste0(wd, cond, ".reprod.bed"),
    quote     = FALSE, sep = "\t",
    row.names = FALSE, col.names = FALSE
  )

  message("  Written: ", cond, ".reprod.bed  (", nrow(reprod), " peaks)")
}

message("\n[R STEP 2 complete] — run the bedtools sort/merge steps in the shell now.\n")


