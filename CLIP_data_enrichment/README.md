# CLIP Enrichment Analysis

Tests whether the presence of RNA-binding proteins (RBPs) with published CLIP-seq data is
significantly enriched over alternatively spliced events in a specific set of significant events relative to a matched set of background events.

---

## Method

For each RBP represented in a CLIP database, the script computes the fraction
of significant AS events overlapping at least one CLIP peak for that RBP
(the observed frequency). A null distribution is built by drawing
10,000 random samples of background events — matched in number to the
significant set — and computing the same overlap frequency each time.
An empirical p-value is derived as the proportion of permutations that
equal or exceed the observed frequency, and BH FDR correction is applied
across all RBPs.

The analysis is run for two event scopes and two CLIP databases, producing
four result tables:

| Scope | CLIP database | Output file |
|---|---|---|
| Skipped exon (SE) | POSTAR3 | `sig_vs_back_se_allclips.txt` |
| Skipped exon (SE) | ENCODE | `sig_vs_back_se_encode.txt` |
| All AS event types | POSTAR3 | `sig_vs_back_allevents_allclips.txt` |
| All AS event types | ENCODE | `sig_vs_back_allevents_encode.txt` |

All four tables are also saved together as `clip_enrich.RData`.

---

## Dependencies

### R packages
```r
BiocManager::install(c("GenomicRanges", "Biostrings"))
install.packages(c("dplyr", "optparse"))
```

---

## Input data format

The script takes two RData files as arguments.

**Landscape file** (`--landscape`): must contain five filtered data.frames:
`sig_se`, `sig_ri`, `sig_mxe`, `sig_a5ss`, `sig_a3ss`

**Background file** (`--background`): must contain the five matched background data.frames objects:
`back_se`, `back_ri`, `back_mxe`, `back_a5ss`, `back_a3ss`

Each data.frame requires at minimum the full coordinates of the event (first 12 columns of the rMATS outputs).

---

## Usage
```bash
Rscript enrich_clip.R \
  --landscape  /path/to/landscape.RData \
  --background /path/to/background.RData \
  --outdir     /path/to/output/directory
```

The output directory must already exist.

### Recommended SLURM submission

The permutation step (10,000 replicates × 4 analyses) is the computational
bottleneck. 

```bash
#!/bin/bash
#SBATCH --job-name=clip_enrich
#SBATCH --partition=all_5hrs
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --ntasks=1
#SBATCH --output=/path/to/your/Logs/clip_enrich.%j.log

module load R/4.2.2-foss-2022b

Rscript /path/to/enrich_clip.R \
  --landscape  /path/to/landscape.RData \
  --background /path/to/background.RData \
  --outdir     /path/to/output/directory
```

---

## Output columns

Each of the four tab-separated output files contains one row per RBP with
the following columns:

| Column | Description |
|---|---|
| `freq_sig` | Observed overlap frequency in significant events |
| `freq_random` | Mean overlap frequency across permuted background samples |
| `dif` | Absolute difference (`freq_sig − freq_random`) |
| `rel_dif` | Relative difference (`dif / freq_sig`) |
| `FDR` | BH-adjusted p-value |

---

## Configuration

Two variables at the top of `enrich_clip.R` should be set once for your
environment:
```r
CLIP_ALL_RDATA    <- "/path/to/your/resources/clip_data/postar3/all_clip_ranges.RData"
CLIP_ENCODE_RDATA <- "/path/to/your/resources/clip_data/postar3/encode_ranges.RData"
```

`N_PERMUTATIONS` (default: 10,000) controls the resolution of the null
distribution and can be reduced for testing.

---

## Contributors

**Loïc Ongena** - loic.ongena@uliege.be

