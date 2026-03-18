#!/usr/bin/env Rscript
# =============================================================================
# enrich_clip.R
# Computes RBP CLIP enrichment over a splicing landscape versus a matched
# background, using a permutation-based FDR approach.
#
#
# Analyses are run for:
#   - Skipped exon (SE) events only
#   - All AS event types combined (SE, RI, MXE, A5SS, A3SS)
# Each is tested against two CLIP databases:
#   - POSTAR3 (all_clip_ranges)
#   - ENCODE  (encode_ranges)
#
# Usage:
#   Rscript enrich_clip.R \
#     -l /path/to/landscape.RData \
#     -b /path/to/background.RData \
#     -o /path/to/output/directory
#
# Input RData objects:
#   landscape  — must contain: sig_se, sig_ri, sig_mxe,
#                              sig_a5ss, sig_a3ss
#   background — must contain: back_se, back_ri, back_mxe,
#                              back_a5ss, back_a3ss
#
# Each table must have columns: chr, strand, upstreamES, downstreamEE
# A5SS/A3SS tables require: longExonStart_0base, flankingEE
#
# Output files written to --outdir:
#   clip_enrich.RData                      — all four result data.frames
#   sig_vs_back_se_allclips.txt            — SE events vs POSTAR3
#   sig_vs_back_se_encode.txt              — SE events vs ENCODE
#   sig_vs_back_allevents_allclips.txt     — all AS events vs POSTAR3
#   sig_vs_back_allevents_encode.txt       — all AS events vs ENCODE
# =============================================================================

suppressWarnings(suppressMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(dplyr)
  library(optparse)
}))


# -----------------------------------------------------------------------------
# Configuration — paths to CLIP resource RData files
# These are stable resources shared across analyses; update once per environment
# -----------------------------------------------------------------------------

CLIP_ALL_RDATA    <- "/path/to/your/resources/clip_data/postar3/all_clip_ranges.RData"
CLIP_ENCODE_RDATA <- "/path/to/your/resources/clip_data/postar3/encode_ranges.RData"

# Number of permutation replicates for the background null distribution
N_PERMUTATIONS <- 10000


# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------

option_list <- list(
  make_option(c("-l", "--landscape"),  type = "character", default = NULL,
    help = "Path to landscape RData file (sig_se, sig_ri, etc.)",
    metavar = "FILE"),
  make_option(c("-b", "--background"), type = "character", default = NULL,
    help = "Path to background RData file (back_se, back_ri, etc.)",
    metavar = "FILE"),
  make_option(c("-o", "--outdir"),     type = "character", default = NULL,
    help = "Output directory (must exist)",
    metavar = "DIR")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$landscape) || is.null(opt$background) || is.null(opt$outdir)) {
  print_help(OptionParser(option_list = option_list))
  stop("All three arguments (--landscape, --background, --outdir) are required.",
       call. = FALSE)
}


# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------

message("Loading landscape and background data...")
load(opt$landscape)   # → sig_se, sig_ri, sig_mxe, sig_a5ss, sig_a3ss
load(opt$background)  # → back_se, back_ri, back_mxe, back_a5ss, back_a3ss

message("Loading CLIP reference ranges...")
load(CLIP_ALL_RDATA)    # → all_clip_ranges
load(CLIP_ENCODE_RDATA) # → encode_ranges


# -----------------------------------------------------------------------------
# Build GRanges objects for each event set
#
# For A5SS and A3SS the longer exon boundaries differ by strand:
#   +  strand: event spans longExonStart_0base → flankingEE
#   -  strand: event spans flankingEE          → longExonStart_0base
# -----------------------------------------------------------------------------

#' Build a GRanges covering the full window of each AS event
#'
#' @param se, ri, mxe, a5ss, a3ss  data.frames with AS event coordinates
#' @return GRanges with an `as_type` metadata column
build_as_ranges <- function(se, ri, mxe, a5ss, a3ss) {

  a5ss_start <- ifelse(a5ss$strand == "+", a5ss$longExonStart_0base, a5ss$flankingEE)
  a5ss_end   <- ifelse(a5ss$strand == "+", a5ss$flankingEE,          a5ss$longExonStart_0base)
  a3ss_start <- ifelse(a3ss$strand == "+", a3ss$flankingEE,          a3ss$longExonStart_0base)
  a3ss_end   <- ifelse(a3ss$strand == "+", a3ss$longExonStart_0base, a3ss$flankingEE)

  GRanges(
    seqnames = c(se$chr,  ri$chr,  mxe$chr,  a5ss$chr,  a3ss$chr),
    ranges   = IRanges(
      start = c(se$upstreamES, ri$upstreamES, mxe$upstreamES, a5ss_start, a3ss_start),
      end   = c(se$downstreamEE, ri$downstreamEE, mxe$downstreamEE, a5ss_end, a3ss_end)
    ),
    strand   = c(se$strand, ri$strand, mxe$strand, a5ss$strand, a3ss$strand),
    as_type  = c(
      rep("se",   nrow(se)),
      rep("ri",   nrow(ri)),
      rep("mxe",  nrow(mxe)),
      rep("a5ss", nrow(a5ss)),
      rep("a3ss", nrow(a3ss))
    )
  )
}

# Skipped-exon only (most common and biologically interpretable event type)
ranges_se_sig  <- GRanges(
  seqnames = sig_se$chr,
  ranges   = IRanges(start = sig_se$upstreamES, end = sig_se$downstreamEE),
  strand   = sig_se$strand
)
ranges_se_back <- GRanges(
  seqnames = back_se$chr,
  ranges   = IRanges(start = back_se$upstreamES, end = back_se$downstreamEE),
  strand   = back_se$strand
)

# All AS event types combined
ranges_all_sig  <- build_as_ranges(sig_se, sig_ri, sig_mxe,
                                    sig_a5ss, sig_a3ss)
ranges_all_back <- build_as_ranges(back_se, back_ri, back_mxe,
                                    back_a5ss, back_a3ss)


# -----------------------------------------------------------------------------
# Core enrichment analysis
#
# Strategy:
#   1. Find which RBPs overlap each significant event (deduplicated per event).
#   2. Compute observed overlap frequency per RBP = overlapping events / total events.
#   3. Build a null distribution by sampling N_PERMUTATIONS sets of background events
#      of the same size as the significant set, computing the same frequency each time.
#   4. Empirical p-value = proportion of permutations exceeding the observed frequency.
#   5. Apply BH FDR correction across all RBPs.
# -----------------------------------------------------------------------------

#' Summarise enrichment results as a data.frame
#'
#' @param hit_back   matrix: RBPs × permutations, each cell = background overlap frequency
#' @param hit_freq   named numeric vector: observed overlap frequency per RBP
#' @param rbp_names  character vector of RBP names (row order)
#' @return data.frame with columns: freq_sig, freq_random, dif, rel_dif, FDR
summarise_enrichment <- function(hit_back, hit_freq, rbp_names) {
  pval <- vapply(seq_along(rbp_names),
    function(i) mean(hit_back[i, ] > hit_freq[i]),
    numeric(1)
  )
  result <- cbind(
    freq_sig   = hit_freq,
    freq_random = rowMeans(hit_back),
    dif        = hit_freq - rowMeans(hit_back),
    rel_dif    = (hit_freq - rowMeans(hit_back)) / hit_freq,
    FDR        = p.adjust(pval, method = "fdr")
  )
  return(as.data.frame(result))
}


#' Compute RBP CLIP enrichment for a set of significant AS events
#'
#' @param ranges_sig   GRanges of significant AS events
#' @param ranges_back  GRanges of matched background AS events
#' @param clip_ranges  GRanges of RBP CLIP peaks with a `$rbp` metadata column
#' @return data.frame of enrichment statistics (one row per RBP)
get_clip_enrich <- function(ranges_sig, ranges_back, clip_ranges) {

  # Observed: deduplicated (event, RBP) pairs to avoid inflating counts
  # when a single event overlaps multiple peaks from the same RBP
  olap_sig       <- findOverlaps(ranges_sig, clip_ranges)
  olap_sig_df    <- as.data.frame(olap_sig)
  olap_sig_df$rbp <- clip_ranges$rbp[subjectHits(olap_sig)]
  olap_sig_dedup <- distinct(olap_sig_df[, c("queryHits", "rbp")])
  freq_rbp_sig   <- table(olap_sig_dedup$rbp) / length(ranges_sig)

  # Background: same deduplication, then permute
  olap_back       <- findOverlaps(ranges_back, clip_ranges)
  olap_back_df    <- as.data.frame(olap_back)
  olap_back_df$rbp <- clip_ranges$rbp[subjectHits(olap_back)]
  olap_back_dedup  <- distinct(olap_back_df[, c("queryHits", "rbp")])

  message("  Overlaps computed. Running ", N_PERMUTATIONS, " permutations...")

  perm_results <- replicate(N_PERMUTATIONS, {
    sampled_ids  <- sample(unique(olap_back_dedup$queryHits), length(ranges_sig))
    sampled_olap <- olap_back_dedup[olap_back_dedup$queryHits %in% sampled_ids, ]
    return(table(sampled_olap$rbp) / length(ranges_sig))
  })

  # Align all RBP names across permutations (some may be absent in some samples)
  all_names      <- sort(unique(unlist(lapply(perm_results, names))))
  perm_aligned   <- lapply(perm_results, function(v) v[all_names])
  back_mat       <- do.call(cbind, perm_aligned)
  back_mat[is.na(back_mat)] <- 0
  rownames(back_mat) <- all_names

  freq_rbp_sig        <- freq_rbp_sig[all_names]
  names(freq_rbp_sig) <- all_names

  message("  Permutations done. Computing FDR...")
  return(summarise_enrichment(back_mat, freq_rbp_sig, all_names))
}


# -----------------------------------------------------------------------------
# Run all four analyses
# -----------------------------------------------------------------------------

message("Analysis 1/4 — SE events vs all POSTAR3 CLIP...")
sig_vs_back_se_all    <- get_clip_enrich(ranges_se_sig,  ranges_se_back,  all_clip_ranges)

message("Analysis 2/4 — SE events vs ENCODE CLIP...")
sig_vs_back_se_encode <- get_clip_enrich(ranges_se_sig,  ranges_se_back,  encode_ranges)

message("Analysis 3/4 — All AS events vs all POSTAR3 CLIP...")
sig_vs_back_all       <- get_clip_enrich(ranges_all_sig, ranges_all_back, all_clip_ranges)

message("Analysis 4/4 — All AS events vs ENCODE CLIP...")
sig_vs_back_encode    <- get_clip_enrich(ranges_all_sig, ranges_all_back, encode_ranges)


# -----------------------------------------------------------------------------
# Write outputs
# -----------------------------------------------------------------------------

save(
  sig_vs_back_se_all, sig_vs_back_se_encode,
  sig_vs_back_all,    sig_vs_back_encode,
  file = file.path(opt$outdir, "clip_enrich.RData")
)

write.table(sig_vs_back_se_all,    file.path(opt$outdir, "sig_vs_back_se_allclips.txt"),
            quote = FALSE, sep = "\t")
write.table(sig_vs_back_se_encode, file.path(opt$outdir, "sig_vs_back_se_encode.txt"),
            quote = FALSE, sep = "\t")
write.table(sig_vs_back_all,       file.path(opt$outdir, "sig_vs_back_allevents_allclips.txt"),
            quote = FALSE, sep = "\t")
write.table(sig_vs_back_encode,    file.path(opt$outdir, "sig_vs_back_allevents_encode.txt"),
            quote = FALSE, sep = "\t")

message("All done. Results written to: ", opt$outdir)



