# Splicing Map — Sliding-Window Motif Enrichment

Identifies RNA-binding protein (RBP) binding motifs and pentamers that are
positionally enriched around exon-intron junctions in a set of alternatively
spliced events, relative to a matched background. Results are visualised as
splicing-map plots showing per-position motif occurrence alongside a
significance track.


---

## Method overview

Each junction FASTA file contains one sequence per splicing event, aligned
so that position 0 corresponds to the exon-intron boundary. A window of fixed
size (default 50 bp) is stepped one position at a time across the full
sequence length. At each position, the script counts how many sequences
contain a hit for each motif (PWM or exact pentamer) within that window.

This produces, for every motif, a matrix of hit counts with dimensions
*sequences × window positions*. Enrichment over background is then tested
per window position with a one-sided Wilcoxon rank-sum test, FDR-corrected
across positions within each motif. The resulting per-position means, SEMs,
and FDR values are plotted as splicing maps.

The regions to be plotted on the splicing maps directly depend on the width of junction sequences provided for the
`sliding window.R` script. Please adapt depending on this. Typically (see example in the file), the windows covering the 
50 bp immediately flanking the 0 position are excluded from plots to avoid plotting the uninformative boundary sequence biases.

---

## Naming conventions

In the examples provided here, the `upstream` / `s3` refer to the splicing acceptor sites (3′ SS) of the upstream and cassette exons, respectively. The `downstream` / `s5` refer to the donor sites (5′ SS) of the downstream and cassette exons. Moreover, these scripts are based on specific names for significant and background sequences in the input `.RData`. Please adapt the script depending on your own data frame names. 

---

## Dependencies

### R packages
```r
BiocManager::install(c("Biostrings", "BiocParallel", "universalmotif"))
install.packages(c("optparse", "dplyr"))
```

### Resource file

All three scripts share a single PWM resource file. Its path is set in the
**Configuration** block at the top of each script:
```r
PWM_RDATA <- "/path/to/your/resources/motifs/all_combined_pwms.RData"
```

This file must contain a named list called `combined_pwms`. PWM names are
expected to follow the convention `RBPNAME.N` (e.g. `FUS.1`, `FUS.2`), where
the numeric suffix distinguishes multiple PWMs for the same RBP.

---

## Workflow

### Step 1 — `sliding_window.R`

Run once per junction per sequence set (significant events and background).
Each call produces a single RData file containing either `results_pwm` or
`results_kmer`.
```bash
Rscript sliding_window.R \
  --fasta    /path/to/upstream_sig.fasta \
  --prefix   upstream_sig \
  --type     pwm \
  --cpus     16 \
  --window   50 \
  --step     1 \
  --outdir   /path/to/sw_counts
```

| Argument | Required | Description |
|---|---|---|
| `--fasta` | yes | FASTA file; all sequences must be the same length |
| `--prefix` | yes | Base name for the output file |
| `--type` | yes | `pwm` for RBP PWMs or `5mer` for pentamers |
| `--cpus` | yes | Number of parallel workers |
| `--window` | no | Window size in bp (default: 50) |
| `--step` | no | Step size in bp (default: 1) |
| `--outdir` | yes | Output directory |

**Output:** `<prefix>_pwms.RData` or `<prefix>_5mer.RData`

Run this script for all four junctions × two sequence sets (significant +
background), giving 8 RData files in total.

---

### Step 2 — `treating_SW_data.R`

Run once per junction, pairing the significant and background RData files
from Step 1.
```bash
Rscript treating_SW_data.R \
  --data       /path/to/sw_counts/upstream_sig_pwms.RData \
  --background /path/to/sw_counts/upstream_back_pwms.RData \
  --prefix     sig_upstream_vs_back_upstream \
  --cpus       16 \
  --outdir     /path/to/enrichment
```

| Argument | Required | Description |
|---|---|---|
| `--data` | yes | Significant-set RData from Step 1 |
| `--background` | yes | Background-set RData from Step 1 |
| `--prefix` | yes | Base name for output files |
| `--cpus` | yes | Number of parallel workers |
| `--outdir` | yes | Output directory |

**Output files:**

| File | Description |
|---|---|
| `<prefix>_enrichment_SW.RData` | Full per-position statistics for every motif; input to Step 3 |
| `<prefix>_summary_motifs.txt` | Minimum FDR per PWM across all positions (PWM mode only) |
| `<prefix>_summary_rbps.txt` | Minimum FDR per RBP (minimum across all its PWMs) (PWM mode only) |
| `<prefix>_summary_5mers.txt` | Minimum FDR per pentamer (5mer mode only) |

The output RData file must be named following the convention
`sig_<junction>_vs_back_<junction>_enrichment_SW.RData` for
`plot_pwm_enrich.R` to locate it automatically.

---

### Step 3 — `plot_pwm_enrich.R`

Run once per RBP per junction to render a PDF splicing map.
```bash
Rscript plot_pwm_enrich.R \
  --rbp      FUS \
  --enrich   /path/to/enrichment \
  --junction downstream \
  --outdir   /path/to/plots
```

| Argument | Required | Description |
|---|---|---|
| `--rbp` | yes | RBP name (matched against PWM names by `grep`) |
| `--enrich` | yes | Directory containing Step 2 RData output files |
| `--junction` | yes | `downstream`, `s5`, `upstream`, or `s3` |
| `--outdir` | yes | Directory where the PDF will be written |

**Output:** `plot_<RBP>_<junction>.pdf` — one page per PWM associated with
the RBP, each showing:
- Mean per-window motif occurrence in the significant set (red) and
  background (black), with 95% CI shaded ribbons
- −log10(FDR) overlaid on a right-hand axis (dashed dark blue), with a
  horizontal line at the FDR = 0.05 threshold

---

## Recommended SLURM submission

Steps 1 and 2 are the computationally intensive steps. Memory scales with the
number of sequences × window positions × motifs; 32–64 GB is typical for a
full PWM scan. Step 3 can run locally.
```bash
#!/bin/bash
#SBATCH --job-name=sw_motif
#SBATCH --partition=all_24hrs
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --output=/path/to/Logs/sw_%j.log

module load R/4.2.2-foss-2022b

JUNCTION=upstream
OUTDIR_SW=/path/to/sw_counts
OUTDIR_ENRICH=/path/to/enrichment

# Step 1 — significant set
Rscript sliding_window.R \
  --fasta   /path/to/${JUNCTION}_sig.fasta \
  --prefix  ${JUNCTION}_sig \
  --type    pwm --cpus 16 --outdir ${OUTDIR_SW}

# Step 1 — background set
Rscript sliding_window.R \
  --fasta   /path/to/${JUNCTION}_back.fasta \
  --prefix  ${JUNCTION}_back \
  --type    pwm --cpus 16 --outdir ${OUTDIR_SW}

# Step 2 — enrichment testing
Rscript treating_SW_data.R \
  --data       ${OUTDIR_SW}/${JUNCTION}_sig_pwms.RData \
  --background ${OUTDIR_SW}/${JUNCTION}_back_pwms.RData \
  --prefix     sig_${JUNCTION}_vs_back_${JUNCTION} \
  --cpus 16 --outdir ${OUTDIR_ENRICH}
```

---

## Contributors

**Loïc Ongena** - loic.ongena@uliege.be
