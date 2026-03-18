#!/usr/bin/env Rscript
# =============================================================================
# sliding_window.R
# Sliding-window motif occurrence scoring across exon-intron junction sequences.
#
# For each sequence in a FASTA file, a fixed-size window is stepped across the
# full length. Within each window position, motif hits are counted using either:
#   - PWM matching (matchPWM, min score 80%) for RBP motif analysis
#   - Exact k-mer matching (matchPattern) for pentamer analysis
#
# The result is a matrix of per-window hit counts for every motif × sequence,
# which is then passed to treating_SW_data.R for enrichment testing.
#
# Designed to be run on the cluster; parallelised over motifs via BiocParallel.
#
# Usage:
#   Rscript sliding_window.R \
#     --fasta    sequences.fasta \
#     --prefix   upstream_sig \
#     --type     pwm \
#     --cpus     16 \
#     --window   50 \
#     --step     1 \
#     --outdir   /path/to/output
# =============================================================================

suppressWarnings(suppressMessages(library(Biostrings)))
library(BiocParallel)
library(optparse)


# -----------------------------------------------------------------------------
# Configuration — update once per environment
# -----------------------------------------------------------------------------

# Path to the combined PWM RData object (loads 'combined_pwms')
PWM_RDATA <- "/path/to/your/resources/motifs/all_combined_pwms.RData"

# Minimum PWM match score threshold
PWM_MIN_SCORE <- "80%"


# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------

option_list <- list(
  make_option(c("-f", "--fasta"),  type = "character", default = NULL,
    help = "Sequence file in FASTA format", metavar = "FILE"),
  make_option(c("-p", "--prefix"), type = "character", default = NULL,
    help = "Prefix for the output file name", metavar = "STRING"),
  make_option(c("-t", "--type"),   type = "character", default = NULL,
    help = "Motif type: 'pwm' (RBP PWMs) or '5mer' (pentamers)", metavar = "STRING"),
  make_option(c("-c", "--cpus"),   type = "numeric",   default = NULL,
    help = "Number of CPUs to use for parallelisation", metavar = "INT"),
  make_option(c("-w", "--window"), type = "numeric",   default = 50,
    help = "Sliding window size in bp [default: %default]", metavar = "INT"),
  make_option(c("-s", "--step"),   type = "numeric",   default = 1,
    help = "Step size between consecutive windows in bp [default: %default]",
    metavar = "INT"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
    help = "Output directory (must exist)", metavar = "DIR")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$fasta)  || is.null(opt$prefix) || is.null(opt$type) ||
    is.null(opt$cpus)   || is.null(opt$outdir)) {
  print_help(OptionParser(option_list = option_list))
  stop("All required arguments (--fasta, --prefix, --type, --cpus, --outdir) must be supplied.",
       call. = FALSE)
}

if (!opt$type %in% c("pwm", "5mer")) {
  stop("--type must be either 'pwm' or '5mer'. See --help.", call. = FALSE)
}


# -----------------------------------------------------------------------------
# Load sequences and motifs
# -----------------------------------------------------------------------------

message("Loading sequences from: ", opt$fasta)
myseqs <- unique(readDNAStringSet(opt$fasta))
message("  ", length(myseqs), " unique sequences loaded (length: ",
        width(myseqs)[1], " bp)")

load(PWM_RDATA)   # loads: combined_pwms

# Pre-generate all 4^5 = 1024 pentamers for k-mer mode
all_5mer <- apply(
  expand.grid(rep(list(c("A", "C", "T", "G")), 5)),
  1, paste, collapse = ""
)


# -----------------------------------------------------------------------------
# Window definitions
# These are shared by both counting functions via lexical scoping
# -----------------------------------------------------------------------------

win_size      <- opt$window
step_size     <- opt$step
n_windows     <- width(myseqs)[1] - win_size + 1
windows_ranges <- IRanges(
  start = seq(1, n_windows, by = step_size),
  width = win_size
)

message("Window size: ", win_size, " bp | Step: ", step_size,
        " bp | Windows per sequence: ", length(windows_ranges))


# -----------------------------------------------------------------------------
# Counting functions
#
# Both functions return a matrix: sequences (rows) × window positions (columns).
# A hit at window position w is counted if the motif match overlaps the window
# by at least the full motif length (minoverlap = match_length), ensuring only
# complete motif occurrences are scored within each window.
# -----------------------------------------------------------------------------

#' Count PWM hits per sliding window position across a set of sequences
#'
#' @param pwm         a PWM matrix (from combined_pwms)
#' @param seqs2count  DNAStringSet
#' @return            matrix: sequences × windows
count_pwm <- function(pwm, seqs2count) {
  match_length <- ncol(pwm)
  hits     <- lapply(seqs2count, function(x) matchPWM(pwm, x, min.score = PWM_MIN_SCORE))
  hit_irs  <- lapply(hits, IRanges)
  counts   <- lapply(hit_irs, function(x)
    countOverlaps(windows_ranges, x, minoverlap = match_length))
  t(do.call(rbind, counts))
}


#' Count exact k-mer hits per sliding window position across a set of sequences
#'
#' @param kmer        character string (k-mer sequence)
#' @param seqs2count  DNAStringSet
#' @return            matrix: sequences × windows
count_kmer <- function(kmer, seqs2count) {
  match_length <- nchar(kmer)
  hits     <- lapply(seqs2count, function(x) matchPattern(kmer, x))
  hit_irs  <- lapply(hits, IRanges)
  counts   <- lapply(hit_irs, function(x)
    countOverlaps(windows_ranges, x, minoverlap = match_length))
  t(do.call(rbind, counts))
}


# -----------------------------------------------------------------------------
# Run sliding window (parallelised over motifs / k-mers)
# -----------------------------------------------------------------------------

param <- MulticoreParam(workers = opt$cpus, progressbar = TRUE)
message("Starting sliding window analysis (type: ", opt$type, ")...")

if (opt$type == "pwm") {

  results_pwm <- bplapply(combined_pwms, count_pwm, myseqs, BPPARAM = param)
  names(results_pwm) <- names(combined_pwms)

  out_file <- file.path(opt$outdir, paste0(opt$prefix, "_pwms.RData"))
  save(results_pwm, file = out_file)

} else {

  results_kmer <- bplapply(all_5mer, count_kmer, myseqs, BPPARAM = param)
  names(results_kmer) <- all_5mer

  out_file <- file.path(opt$outdir, paste0(opt$prefix, "_5mer.RData"))
  save(results_kmer, file = out_file)

}

message("Done. Output written to: ", out_file)


