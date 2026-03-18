#!/usr/bin/env Rscript
# =============================================================================
# treating_SW_data.R
# Enrichment testing on sliding-window motif counts.
#
# Takes the per-window count matrices produced by sliding_window.R for a
# significant sequence set and a matched background set, and for each motif
# and each window position tests whether counts are higher in the significant
# set using a one-sided Wilcoxon rank-sum test. FDR correction is applied
# per motif across all window positions.
#
# Produces:
#   - An RData file with the full per-position statistics for every motif
#     (mean, SEM, FDR), consumed by plot_pwm_enrich.R
#   - Summary tables of the minimum FDR per motif and per RBP (PWM mode),
#     or per pentamer (5mer mode), for quick ranking of enriched elements
#
# Usage:
#   Rscript treating_SW_data.R \
#     --data       /path/to/sig_pwms.RData \
#     --background /path/to/back_pwms.RData \
#     --prefix     upstream_sig_vs_back_upstream \
#     --cpus       16 \
#     --outdir     /path/to/output
# =============================================================================

suppressWarnings(suppressMessages(library(dplyr)))
library(BiocParallel)
library(optparse)


# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------

option_list <- list(
  make_option(c("-d", "--data"),       type = "character", default = NULL,
    help = "Path to significant-set RData file (output of sliding_window.R)",
    metavar = "FILE"),
  make_option(c("-b", "--background"), type = "character", default = NULL,
    help = "Path to background-set RData file (output of sliding_window.R)",
    metavar = "FILE"),
  make_option(c("-p", "--prefix"),     type = "character", default = NULL,
    help = "Prefix for output file names", metavar = "STRING"),
  make_option(c("-c", "--cpus"),       type = "numeric",   default = NULL,
    help = "Number of CPUs to use", metavar = "INT"),
  make_option(c("-o", "--outdir"),     type = "character", default = NULL,
    help = "Output directory (must exist)", metavar = "DIR")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$data) || is.null(opt$background) || is.null(opt$prefix) ||
    is.null(opt$cpus) || is.null(opt$outdir)) {
  print_help(OptionParser(option_list = option_list))
  stop("All five arguments are required.", call. = FALSE)
}


# -----------------------------------------------------------------------------
# Load data
#
# Both RData files must contain either results_pwm or results_kmer; the type
# must be the same in both files (no mixing of PWM and k-mer results).
# -----------------------------------------------------------------------------

#' Load a sliding_window.R output file and return the result list
#'
#' @param path  path to an RData file containing results_pwm or results_kmer
#' @return      named list of per-motif count matrices
load_sw_result <- function(path) {
  e <- new.env()
  load(path, envir = e)
  if (exists("results_pwm",  envir = e)) return(list(data = e$results_pwm,  type = "pwm"))
  if (exists("results_kmer", envir = e)) return(list(data = e$results_kmer, type = "5mer"))
  stop("File '", path, "' must contain either results_pwm or results_kmer.", call. = FALSE)
}

message("Loading significant-set data from: ", opt$data)
sig_loaded  <- load_sw_result(opt$data)
res_sig     <- sig_loaded$data
motif_type  <- sig_loaded$type

message("Loading background data from: ", opt$background)
back_loaded <- load_sw_result(opt$background)
res_back    <- back_loaded$data

if (sig_loaded$type != back_loaded$type) {
  stop("Significant and background files contain different motif types (",
       sig_loaded$type, " vs ", back_loaded$type, ").", call. = FALSE)
}

n_motifs <- names(res_sig)
message("  ", length(n_motifs), " motifs to test across ",
        ncol(res_sig[[1]]), " window positions.")


# -----------------------------------------------------------------------------
# Per-window Wilcoxon test (parallelised over motifs)
#
# For each window position, a one-sided Wilcoxon rank-sum test asks whether
# motif counts in the significant set are greater than in the background.
# BH FDR correction is then applied across all window positions for that motif.
# -----------------------------------------------------------------------------

#' Run per-position Wilcoxon tests for one motif
#'
#' @param motif_name  character, name of the motif (key into res_sig / res_back)
#' @return            numeric vector of BH-adjusted p-values, one per window position
test_one_motif <- function(motif_name) {
  mat_sig  <- res_sig [[motif_name]]
  mat_back <- res_back[[motif_name]]
  raw_p <- vapply(seq_len(ncol(mat_sig)), function(i) {
    wilcox.test(mat_sig[, i], mat_back[, i], alternative = "greater")$p.value
  }, numeric(1))
  p.adjust(raw_p, method = "fdr")
}

message("Computing per-window means and SEMs...")
res_means_sig  <- lapply(res_sig,  colMeans)
res_means_back <- lapply(res_back, colMeans)
res_sem_sig    <- lapply(res_sig,  function(x) apply(x, 2, sd) / sqrt(nrow(x)))
res_sem_back   <- lapply(res_back, function(x) apply(x, 2, sd) / sqrt(nrow(x)))

message("Running Wilcoxon tests (parallelised over motifs)...")
param     <- MulticoreParam(workers = opt$cpus, progressbar = TRUE)
mot_fdr   <- bplapply(n_motifs, test_one_motif, BPPARAM = param)
names(mot_fdr) <- n_motifs


# -----------------------------------------------------------------------------
# Assemble the per-motif result list consumed by plot_pwm_enrich.R
# Each element is a data.frame: one row per window position
# -----------------------------------------------------------------------------

fin_list <- setNames(
  lapply(seq_along(n_motifs), function(i) {
    nm <- n_motifs[i]
    data.frame(
      Average_sig        = res_means_sig [[nm]],
      Average_background = res_means_back[[nm]],
      SEM_sig            = res_sem_sig   [[nm]],
      SEM_back           = res_sem_back  [[nm]],
      FDR                = mot_fdr[[nm]],
      logFDR             = -log10(mot_fdr[[nm]])
    )
  }),
  n_motifs
)

fdr_min <- vapply(fin_list, function(x) min(x$FDR), numeric(1))


# -----------------------------------------------------------------------------
# Write outputs
# -----------------------------------------------------------------------------

out_rdata <- file.path(opt$outdir, paste0(opt$prefix, "_enrichment_SW.RData"))
save(fin_list, file = out_rdata)
message("Saved: ", out_rdata)

if (motif_type == "pwm") {

  # Strip the trailing numeric index (e.g. "RBPNAME.1") to recover the RBP name
  rbp_names <- sub("[._][0-9]+$", "", names(fdr_min))

  summary_motifs <- data.frame(
    Motif_name     = names(fdr_min),
    RBP_name       = rbp_names,
    FDR_min_motif  = fdr_min
  )
  summary_rbps <- summary_motifs %>%
    group_by(RBP_name) %>%
    summarise(FDR_min_rbp = min(FDR_min_motif), .groups = "drop")

  write.table(summary_motifs,
    file.path(opt$outdir, paste0(opt$prefix, "_summary_motifs.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  write.table(summary_rbps,
    file.path(opt$outdir, paste0(opt$prefix, "_summary_rbps.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

} else {

  summary_5mers <- data.frame(
    Pentamer      = names(fdr_min),
    FDR_min_motif = fdr_min
  )
  write.table(summary_5mers,
    file.path(opt$outdir, paste0(opt$prefix, "_summary_5mers.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

message("All done.")
