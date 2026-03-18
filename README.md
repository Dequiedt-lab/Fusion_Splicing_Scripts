# Homemade scripts - FET oncogenic fusions in the regulation of alternative splicing
This repository contains homemade original scripts developped for the investigation of the role of FET oncogenic fusions in the regulation of pre-mRNA splicing.

## Repository structure

```

├── ePRINT_analysis/              # eCLIP pipeline — from FASTQ to peaks and motif analysis
│   ├── pipeline/                 #   FASTQ → deduplicated BAMs → CLIPper peak calling
│   └── peak_processing/          #   Peak reproducibility, differential binding, motifs
├── CLIP_data_enrichment/         # RBP CLIP enrichment over splicing landscapes
└── Motif_analysis_splicing_maps/ # Sliding-window motif enrichment and splicing-map plots

```

---

## General conventions

**Execution environments.** SLURM-submitted jobs use `sbatch`.
Interactive steps (e.g. bedtools peak processing, bigWig generation) are run
directly in the terminal, sometimes inside a Singularity shell.

**Singularity containers.** All third-party tools are encapsulated in `.sif`
files stored in a shared directory referenced by the `path2sif` variable,
which is defined at the top of each script.

**Path configuration.** Every script defines its input/output paths in a
clearly marked configuration block at the top of the file. Hardcoded
institution-specific paths have been replaced with `/path/to/your/...`
placeholders — update these before running.

**Parallelisation.** Cluster scripts that use BiocParallel expose a `--cpus`
argument. Set this to match the `--cpus-per-task` value in the corresponding
SLURM header.

---

## Contributors

Loïc Ongena - loic.ongena@uliege.be
Main contact: Franck Dequiedt - fdequiedt@uliege.be
