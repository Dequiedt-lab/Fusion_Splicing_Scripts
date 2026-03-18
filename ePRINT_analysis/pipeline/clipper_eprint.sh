#!/bin/bash
#
#SBATCH --job-name=clipper
#SBATCH --partition=all_24hrs
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=1
#SBATCH --mail-user=your.email@institution.be
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/path/to/your/eprint/Logs/clipper/clipper.%j.log

printf " ############################################# \n # Starting CLIPper for $prefix # \n ############################################# \n\n"

source ~/.bashrc
module load SAMtools/1.17-GCC-12.2.0

path2bam=${path2wd}/dedup
path2clipper=${path2wd}/clipper
mkdir -p ${path2clipper}

# ============================================================
# CLIPper — peak calling
# --species GRCh38_v40 : genome/annotation version used internally by CLIPper
# --FDR 0.01           : false discovery rate threshold for peak significance
# Run on r2 only for duplication, as r2 is the informative strand in ePRINT
# ============================================================

printf "\nStarting CLIPper...\n\n"

singularity exec --bind ${path2bam},${path2clipper} \
    ${path2sif}/clipper.sif \
    clipper --species GRCh38_v40 --FDR 0.01 \
    --bam     ${path2bam}/${prefix}.dedup_r2.bam \
    --outfile ${path2clipper}/${prefix}.peakClusters_r2.bed

# ============================================================
# Count mapped reads in deduplicated BAM
# These numbers are used for normalisation in downstream analysis
# ============================================================

printf "\nCounting mapped reads...\n\n"
samtools view -cF 4 ${path2bam}/${prefix}.dedup_r2.bam > ${path2clipper}/${prefix}_readnum_r2.txt

printf "\nAll done!\n\n"

sacct --format="JobId,JobName,NodeList,State,Elapsed,Timelimit,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveVMSize,ReqMem,Submit,Eligible" -j ${SLURM_JOB_ID}


