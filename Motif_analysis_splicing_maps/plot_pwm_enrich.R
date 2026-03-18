#!/usr/bin/env Rscript
# =============================================================================
# plot_pwm_enrich.R
# Renders splicing-map plots for one RBP across a specified junction type.
#
# For each PWM associated with the requested RBP, one page is added to a PDF
# showing:
#   - Mean per-window motif occurrence in the significant set (red) and
#     background (black), with 95% confidence intervals as shaded ribbons
#   - -log10(FDR) overlaid on a right-hand axis (dashed dark blue line) with
#     a significance threshold at -log10(0.05)
#
# Junction coordinate conventions:
#   downstream / s5  — splice donor:   x-axis spans position −250 to +50
#                      (exon body left, intron body right of the junction)
#   upstream   / s3  — splice acceptor: x-axis spans position −50 to +250
#                      (intron body left, exon body right of the junction)
#
# The 50-position gap around position 0 (the exon-intron boundary itself) is
# removed from the plot; the index vectors c(1:250,301:351) and c(1:50,101:351)
# implement this masking.
#
# Usage:
#   Rscript plot_pwm_enrich.R \
#     --rbp      FUS \
#     --enrich   /path/to/enrichment/folder \
#     --junction downstream \
#     --outdir   /path/to/plots
# =============================================================================

suppressWarnings(suppressMessages({
  library(universalmotif)
  library(optparse)
}))


# -----------------------------------------------------------------------------
# Configuration — update once per environment
# -----------------------------------------------------------------------------

# Path to the combined PWM RData object (loads 'combined_pwms')
PWM_RDATA <- "/path/to/your/resources/motifs/all_combined_pwms.RData"

# Significance threshold displayed as a horizontal dashed line on the FDR axis
FDR_THRESHOLD <- 0.05


# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------

option_list <- list(
  make_option(c("-r", "--rbp"),      type = "character", default = NULL,
    help = "RBP name (must match names in combined_pwms)", metavar = "STRING"),
  make_option(c("-e", "--enrich"),   type = "character", default = NULL,
    help = "Directory containing treating_SW_data.R output RData files",
    metavar = "DIR"),
  make_option(c("-j", "--junction"), type = "character", default = NULL,
    help = "Junction type: 'downstream' (donor/s5) or 'upstream' (acceptor/s3)",
    metavar = "STRING"),
  make_option(c("-o", "--outdir"),   type = "character", default = NULL,
    help = "Output directory for PDF plots", metavar = "DIR")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$rbp) || is.null(opt$enrich) ||
    is.null(opt$junction) || is.null(opt$outdir)) {
  print_help(OptionParser(option_list = option_list))
  stop("All four arguments are required.", call. = FALSE)
}

if (!opt$junction %in% c("downstream", "s5", "upstream", "s3")) {
  stop("--junction must be one of: downstream, s5, upstream, s3.", call. = FALSE)
}


# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------

load(PWM_RDATA)   # loads: combined_pwms

enrich_file <- file.path(opt$enrich,
  paste0("sig_", opt$junction, "_vs_back_", opt$junction, "_enrichment_SW.RData"))

if (!file.exists(enrich_file)) {
  stop("Enrichment file not found: ", enrich_file, call. = FALSE)
}
load(enrich_file)   # loads: fin_list


# -----------------------------------------------------------------------------
# Identify PWMs belonging to the requested RBP
# -----------------------------------------------------------------------------

pos_rbp <- grep(opt$rbp, names(combined_pwms))
if (length(pos_rbp) == 0) {
  stop("No PWMs found matching RBP name '", opt$rbp,
       "' in combined_pwms.", call. = FALSE)
}
message("Plotting ", length(pos_rbp), " PWM(s) for RBP: ", opt$rbp)


# -----------------------------------------------------------------------------
# Coordinate masking by junction type
#
# The enrichment vectors have 351 positions covering 300 bp of sequence plus
# a 50-position gap at the junction boundary (positions 251-300 for donor,
# 51-100 for acceptor).
# -----------------------------------------------------------------------------

if (opt$junction %in% c("downstream", "s5")) {
  x_coords  <- -250:50            # donor: exon (negative) → intron (positive)
  pos_mask  <- c(1:250, 301:351)  # skip the 50-position gap at position 0
} else {
  x_coords  <- -50:250            # acceptor: intron (negative) → exon (positive)
  pos_mask  <- c(1:50, 101:351)
}


# -----------------------------------------------------------------------------
# Render plots
# -----------------------------------------------------------------------------

out_pdf <- file.path(opt$outdir,
  paste0("plot_", opt$rbp, "_", opt$junction, ".pdf"))
pdf(file = out_pdf, width = 10, height = 6)

for (i in pos_rbp) {

  motif_name <- names(combined_pwms)[i]
  pwm_obj    <- combined_pwms[[i]]
  cons_seq   <- create_motif(pwm_obj)@consensus

  mot_data   <- fin_list[[motif_name]]
  if (is.null(mot_data)) {
    message("  Skipping ", motif_name, " — not found in enrichment data.")
    next
  }

  # Extract masked position vectors
  y_sig   <- mot_data$Average_sig       [pos_mask]
  sem_sig <- mot_data$SEM_sig           [pos_mask]
  y_back  <- mot_data$Average_background[pos_mask]
  sem_back <- mot_data$SEM_back         [pos_mask]
  fdr_vec <- mot_data$logFDR            [pos_mask]

  # 95% confidence intervals (±1.96 SEM)
  ci_lo_sig  <- y_sig  - 1.96 * sem_sig
  ci_hi_sig  <- y_sig  + 1.96 * sem_sig
  ci_lo_back <- y_back - 1.96 * sem_back
  ci_hi_back <- y_back + 1.96 * sem_back

  y_range  <- range(c(ci_lo_sig, ci_hi_sig, ci_lo_back, ci_hi_back))
  fdr_max  <- max(c(fdr_vec, 2))   # ensure axis spans at least up to 2

  # --- Left panel: motif occurrence -----------------------------------------
  op <- par(mar = c(5, 5, 4, 5) + 0.1)

  plot(x_coords, y_sig,
    type = "l", col = "red", lwd = 2,
    ylim = y_range,
    xlab = "Position relative to junction (bp)",
    ylab = "Motif occurrence (mean per window)",
    main = paste("Motif occurrence —", motif_name, "—", cons_seq),
    cex.lab = 1.3, las = 1
  )
  polygon(c(x_coords, rev(x_coords)),
          c(ci_lo_sig, rev(ci_hi_sig)),
          col = adjustcolor("red", alpha.f = 0.2), border = NA)

  lines(x_coords, y_back, col = "black", lwd = 2)
  polygon(c(x_coords, rev(x_coords)),
          c(ci_lo_back, rev(ci_hi_back)),
          col = adjustcolor("black", alpha.f = 0.2), border = NA)

  # --- Right axis: -log10(FDR) overlay --------------------------------------
  par(new = TRUE)
  plot(x_coords, fdr_vec,
    type = "l", col = "darkblue", lwd = 1, lty = 2,
    axes = FALSE, xlab = "", ylab = "",
    ylim = c(0, fdr_max)
  )
  axis(side = 4, col.axis = "darkblue", col = "darkblue")
  abline(h = -log10(FDR_THRESHOLD), lwd = 0.8, lty = 3, col = "darkblue")

  # Right-axis label, rotated to read bottom-to-top
  text(
    x      = par("usr")[2] + 30,
    y      = fdr_max / 2,
    labels = "-log10(FDR)",
    srt = 270, col = "darkblue", cex = 1.2, xpd = TRUE
  )

  par(op)
}

dev.off()
message("PDF written to: ", out_pdf)



