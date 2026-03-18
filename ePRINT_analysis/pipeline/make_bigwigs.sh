#!/bin/bash
# Interactive script — run directly, no SLURM submission needed
# Generates strand-specific normalised bigWig files from deduplicated BAMs
# using the makebigwigfiles tool (Singularity container)
#
# Run from any directory; paths are defined below

# ============================================================
# Update these paths for your environment
# ============================================================

path2bams=${path2wd}/dedup/
path2bigwigs=${path2wd}/bigwigs/
path2sif=/path/to/your/SIF_files
path2chr=/path/to/your/genomes/gencode_star_index/   # directory containing chrNameLength.txt

mkdir -p ${path2bigwigs}

# ============================================================
# Discover all deduplicated r2 BAM files (informative strand)
# ============================================================

cd ${path2bams}
prefix=($(ls *dedup_r2.bam | cut -d'.' -f1))

printf "Found %d BAM files:\n" "${#prefix[@]}"
for t in ${prefix[@]}; do printf "  %s\n" "${t}"; done
printf "\n"

# ============================================================
# Generate bigWigs
# --bw_pos / --bw_neg : output files for + and - strand
# Reads are normalised to library size internally by makebigwigfiles
# ============================================================

for t in ${prefix[@]}; do
    printf "Processing %s...\n" "${t}"
    singularity exec --no-home \
        --bind ${path2bams},${path2bigwigs},${path2chr} \
        ${path2sif}/makebigwigfiles_0.0.3.sif \
        makebigwigfiles \
        --bw_pos ${path2bigwigs}/${t}.norm.pos.bw \
        --bw_neg ${path2bigwigs}/${t}.norm.neg.bw \
        --bam    ${path2bams}/${t}.dedup_r2.bam \
        --genome ${path2chr}/chrNameLength.txt
done

printf "\nAll bigWigs written to %s\n" "${path2bigwigs}"



