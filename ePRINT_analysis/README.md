# ePRINT Analysis Pipeline

End-to-end analysis of [ePRINT](https://doi.org/10.1186/s13059-024-03271-1) data, from raw FASTQ files to reproducible peak sets, motif analysis, and differential binding. The pipeline depends on the latest library prep for ePRINT experiments in which R2 is the informative strand (check before applying to your data).

---

## Directory structure
```
eprint/
│
├── pipeline/                   # FASTQ → deduplicated BAMs → CLIPper peaks
│   ├── master_script.sh        # SLURM submission wrapper — start here
│   ├── umi2dedup.sh            # UMI extraction → trimming → alignment → deduplication
│   ├── qc_umi2dedup.sh         # FastQC / PreSeq / RSeQC quality controls
│   ├── clipper_eprint.sh      # CLIPper peak calling on deduplicated BAMs
│   └── make_bigwigs.sh         # Strand-specific normalised bigWig generation (interactive)
│
└── peak_processing/            # CLIPper peaks → reproducible peaks → motif analysis
    ├── peak_processing.sh      # bedtools workflow (interactive, inside Singularity shell)
    ├── peaks_r_treat.R         # R steps interleaved with bedtools (normalisation + reproducibility)
    ├── peak_analysis.Rmd       # R Markdown: annotations, CLIP validation, differential binding, motifs
    ├── count_pwm.R             # Cluster script: per-position PWM scoring (run before Rmd §9)
    └── peaks_motif_counts.R    # Cluster script: total PWM hit counts per peak category
```

---

## Dependencies

### Software containers (Singularity `.sif` files)

| Container | Used in |
|---|---|
| `star_2.7.11b.sif` | `umi2dedup.sh` |
| `preseq_latest.sif` | `qc_umi2dedup.sh` |
| `rseqc_5.0.1.sif` | `qc_umi2dedup.sh` |
| `clipper.sif` | `clipper_eprint2.sh` |
| `bedtools.sif` | `peak_processing.sh` |
| `makebigwigfiles_0.0.3.sif` | `make_bigwigs.sh` |

### System modules (loaded via `module load`)

- `cutadapt/4.4-GCCcore-12.2.0`
- `SAMtools/1.17-GCC-12.2.0`
- `fastqc/0.12.1`
- `R/4.2.2-foss-2022b`
- `singularity/3.7.1`

### Python environment

`umi_tools` must be installed in a Python virtual environment. The path is set via `VENV_PATH` in `umi2dedup.sh`.

### R packages
```r
BiocManager::install(c(
  "Biostrings", "GenomicRanges", "BSgenome.Hsapiens.UCSC.hg38",
  "rtracklayer", "universalmotif", "ChIPseeker",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "BiocParallel"
))
install.packages(c("pheatmap", "ape"))
```

### External binaries

- `fastq-sort` — path set via `FASTQSORT` in `umi2dedup.sh`

---

## Part 1 — Pipeline

### Execution model

All steps except bigWig generation are submitted to SLURM via `master_script.sh`. Steps run sequentially through job dependencies (`--dependency=afterok`). One invocation of `master_script.sh` processes a single sample; submit it once per sample.
```
master_script.sh
    └─[job 1]─ umi2dedup.sh
    └─[job 2, after job 1]─ qc_umi2dedup.sh
    └─[job 3, after job 2]─ clipper_eprint2.sh
```

`make_bigwigs.sh` is run interactively after all samples are processed.

### Configuration

Open `master_script.sh` and set the following variables before each submission:

| Variable | Description |
|---|---|
| `prefix` | Sample base name (shared by `_R1_001.fastq.gz` and `_R2_001.fastq.gz`) |
| `path2fastq` | Directory containing raw FASTQ files |
| `path2wd` | Working directory for all outputs |
| `path2genomes` | Directory containing STAR indices and annotation files |
| `path2sif` | Directory containing Singularity `.sif` files |
| `path2scripts` | Directory containing the pipeline scripts |

Also update tool-specific paths defined within individual scripts:

- `VENV_PATH` in `umi2dedup.sh` — Python virtual environment with `umi_tools`
- `FASTQSORT` in `umi2dedup.sh` — path to the `fastq-sort` binary
- `--mail-user` in `clipper_eprint.sh` — email address for SLURM notifications

### Genome and annotation files required

The following must exist under `path2genomes/`:
```
genomes/
├── repbase_STARindex/          # STAR index built from RepBase repeat sequences
├── gencode_star_index/         # STAR index built from GRCh38 + GENCODE annotation
│   └── chrNameLength.txt       # Chromosome sizes file
└── gencode.v48.bed             # Gene annotation in BED format (for RSeQC)
```

### Step-by-step overview

**`umi2dedup.sh`** — UMI extraction → adapter trimming → alignment → deduplication

1. **UMI extraction** (`umi_tools extract`): moves the 10 nt UMI from the 5′ end of each read into the read name for later deduplication.
2. **Adapter trimming** (`cutadapt`): removes Illumina TruSeq adapters and partial suffixes; discards reads shorter than 18 nt.
3. **FASTQ sorting** (`fastq-sort`): sorts R1 and R2 by read name to ensure paired order before alignment.
4. **RepBase alignment** (`STAR`): maps reads to repeat elements; unmapped reads are carried forward. Up to 30 multimappers allowed.
5. **Genome alignment** (`STAR`): maps RepBase-unmapped reads to GRCh38/GENCODE. Only uniquely mapping reads are retained (`--outFilterMultimapNmax 1`).
6. **BAM sorting and indexing** (`samtools`): name-sorts then coordinate-sorts each BAM; removes intermediate files.
7. **PCR deduplication** (`umi_tools dedup`): collapses PCR duplicates using UMI + mapping position. Run on both r1 and r2.

**`qc_umi2dedup.sh`** — Quality controls

- **FastQC**: run on raw and trimmed FASTQ files.
- **PreSeq** (`lc_extrap`): library complexity curves on aligned and deduplicated BAMs for r1 and r2.
- **RSeQC**: read distribution across genomic features (r1 and r2, pre- and post-dedup); read duplication plots (r2 only — the informative strand in ePRINT).

**`clipper_eprint2.sh`** — Peak calling

- Runs **CLIPper** on the deduplicated r2 BAM (`GRCh38_v40`, FDR 0.01).
- Counts mapped reads with `samtools view -cF 4`; count is saved to `${prefix}_readnum_r2.txt` for use as the normalisation denominator in peak processing.

**`make_bigwigs.sh`** — bigWig generation (interactive)

- Generates strand-specific, library-size-normalised bigWig files (`.norm.pos.bw` / `.norm.neg.bw`) using `makebigwigfiles`.

### Outputs per sample
```
${path2wd}/
├── umi_extract/    *_1.umi.fastq.gz, *_2.umi.fastq.gz, *.umi.log
├── cutadapt/       *_1.trim.sorted.fastq.gz, *_2.trim.sorted.fastq.gz, *.trim.log
├── star_align/     *.genome_mapped_r{1,2}.sorted.bam(.bai)
├── dedup/          *.dedup_r{1,2}.bam(.bai), *.umi_dedup_stats_r{1,2}.*
├── clipper/        *.peakClusters_r2.bed, *_readnum_r2.txt
├── bigwigs/        *.norm.pos.bw, *.norm.neg.bw
├── fastqc/         raw/ and trimmed/ FastQC reports
├── preseq/         *.ccurve.txt (pre- and post-dedup, r1 and r2)
└── rseqc/          read distribution and duplication reports
```

It is recommended to generate a MultiQC report (see: https://github.com/MultiQC/MultiQC/) on this directory and check all analysis results before going further. 

---

## Part 2 — Peak processing (`peak_processing/`)

### Execution model

This stage runs **interactively** and alternates between shell commands (inside a Singularity `bedtools` shell) and R scripts. The handoff points are marked with `[R STEP N]` in `peak_processing.sh` and with `message()` checkpoints in `peaks_r_treat.R`. The workflow assumes two conditions, each with two biological replicates.
```
peak_processing.sh ──► peaks_r_treat.R ──► peak_processing.sh ──► peaks_r_treat.R
    (steps 1–2)          (R steps 1–2)          (step 3)             (R step 2)
         │
         ▼
peak_processing.sh ──► peak_analysis.Rmd (§6) ──► peak_processing.sh
    (step 3 cont.)       (R step 3: common peaks)      (step 4)
         │
         ▼
peak_analysis.Rmd (§5) ──► count_pwm.R / peaks_motif_counts.R ──► peak_analysis.Rmd (§9–10)
  (sequence extraction)          (cluster, prerequisite)
```

### Configuration

All user-facing variables are defined at the top of each file.

**`peak_processing.sh`**

| Variable | Description |
|---|---|
| `path2sif`, `path2bam`, `path2bed`, `path2out`, `path2genomes` | Directory paths |
| `PREFIX_EP_C*`, `PREFIX_INP_C*` | Per-replicate file prefixes from the pipeline stage |
| `COND1`, `COND2`, `COND1_R1`, etc. | Condition and replicate labels |
| `CHROM_SIZES`, `GENE_BED` | Chromosome sizes file and gene annotation BED |
| `CLIPPER_SUFFIX`, `BAM_SUFFIX` | File suffixes matching pipeline outputs |
| `N_SHUFFLES` | Number of shuffled background sets to generate |

**`peaks_r_treat.R`**

| Variable | Description |
|---|---|
| `wd` | Peak processing output directory |
| `cond1`, `cond2`, `cond*_r*` | Labels (must match the shell script) |
| `reads_ep_*`, `reads_inp_*` | Total mapped reads per library (from `*_readnum_r2.txt`) |
| `MIN_EPRINT_READS`, `MIN_FC` | Reproducibility filter thresholds |

`peak_analysis.Rmd`, `count_pwm.R`, and `peaks_motif_counts.R` each contain an equivalent Configuration block. Condition labels and paths must be kept consistent across all four files.

### Step-by-step overview

**Shell step 1 — Sort, merge, and compute coverage**

For each replicate: sorts CLIPper peak BED files, merges overlapping same-strand peaks, then runs `bedtools coverage` twice to append SMInput and IP read counts to each peak.

**R step 1 — Normalisation and enrichment testing** (`peaks_r_treat.R`)

Adds a pseudocount, converts read counts to RPM, computes fold-change (IP/input), and runs a one-sided binomial test of IP enrichment per peak. Output: `<replicate>.norm.bed`.

**Shell step 2 — Replicate intersection**

`bedtools intersect -s -wo` identifies peaks overlapping between the two replicates of each condition.

**R step 2 — Reproducibility filtering** (`peaks_r_treat.R`)

Retains peak pairs where both replicates satisfy `raw IP reads > MIN_EPRINT_READS` AND `fold-change > MIN_FC`. The retained locus spans the union of both replicate coordinates. Output: `<condition>.reprod.bed`.

**Shell step 3 — Merge reproducible peaks**

Sorts and merges reproducible peaks within each condition, taking the maximum score across overlapping entries.

**R step 3 — Differential binding** (`peak_analysis.Rmd`, §6)

Classifies peaks as condition-exclusive or common by comparing `mean_lfc` across overlapping peaks between conditions (threshold: 2×). Writes `common_peaks.unprocessed.bed`.

**Shell step 4 — Merge common peaks**

Sorts and merges the common peak file.

**Shell step 5 — Random background generation**

`bedtools shuffle` confined to genic regions produces `N_SHUFFLES` randomised coordinate sets, concatenated into `random_peaks.bed`.

**`peak_analysis.Rmd`** — Full analysis notebook

Peak width distributions, genomic annotation (ChIPseeker), cross-validation against POSTAR3 CLIP data with permutation testing, sequence extraction (strict and ±50 bp extended), differential binding, annotation of exclusive and common peak sets, and motif z-score heatmaps.

**`count_pwm.R`** — Per-position PWM scoring (cluster)

Scans extended peak sequences against all PWMs at 80% score threshold using `matchPWM`; computes mean per-position hit counts. Submit to the cluster and wait for completion before knitting the motif section of the Rmd. Output: `pwm_colmeans_ext50.RData`.

**`peaks_motif_counts.R`** — Total PWM hit counts (cluster)

Uses `countPWM` for speed to count total motif occurrences across all peak categories (all peaks, condition-exclusive, common, random). Output: `countsPWM_peaks.RData`.

### Outputs
```
${peak_dir}/
├── <replicate>.cov.bed               bedtools coverage (input + IP reads)
├── <replicate>.norm.bed              normalised peaks with enrichment statistics
├── <condition>.intersect.bed         replicate intersection
├── <condition>.reprod.bed            reproducible peaks (R-filtered)
├── <condition>.reprod.merged.bed     merged reproducible peaks
├── common_peaks.unprocessed.bed      common peaks before merge
├── common_peaks.merged.bed           final common peak set
├── <condition>_exclusive.bed         condition-exclusive peaks
├── random_peaks.bed                  pooled shuffled background regions
├── *_seqs.fasta                      peak sequences (strict)
├── *_seqs_ext50.fasta                peak sequences (±50 bp extended)
├── pwm_colmeans_ext50.RData          per-position PWM mean counts
└── countsPWM_peaks.RData             total PWM counts per peak category
```

NOTE: these scripts should be adapted to what you want to say depending your experimental conditions and your hypotheses.

---

## Contributors

**Loïc Ongena** - loic.ongena@uliege.be
