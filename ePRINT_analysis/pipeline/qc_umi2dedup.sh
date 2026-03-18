#!/bin/bash
#
#SBATCH --job-name=qc_umi2dedup
#SBATCH --partition=all_5hrs
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --ntasks=1
#SBATCH --output=/path/to/your/eprint/Logs/qc_umi2dedup/qc_umi2dedup.%j.log

printf " ##################################################### \n # Starting QC script for $prefix # \n ##################################################### \n\n"

module load fastqc/0.12.1
module load singularity/3.7.1
module load R/4.2.2-foss-2022b

# ============================================================
# FastQC — raw and trimmed reads
# Note: raw FASTQ files are expected at:
#   ${path2fastq}/${prefix}_R1_001.fastq.gz (and _R2_)
# If your raw files use a different naming convention, update
# the filenames below accordingly
# ============================================================

printf "FastQC version: "; fastqc --version
mkdir -p ${path2wd}/fastqc/raw
mkdir -p ${path2wd}/fastqc/trimmed

fastqc ${path2fastq}/${prefix}_R1_001.fastq.gz --threads 4 -o ${path2wd}/fastqc/raw
fastqc ${path2fastq}/${prefix}_R2_001.fastq.gz --threads 4 -o ${path2wd}/fastqc/raw
fastqc ${path2wd}/cutadapt/${prefix}_1.trim.sorted.fastq.gz --threads 4 -o ${path2wd}/fastqc/trimmed
fastqc ${path2wd}/cutadapt/${prefix}_2.trim.sorted.fastq.gz --threads 4 -o ${path2wd}/fastqc/trimmed

printf "Done for FastQC!\n\n"

# ============================================================
# PreSeq — library complexity estimation
# lc_extrap: extrapolates the library yield at higher sequencing depths
# Run on both the aligned BAM (before dedup) and the deduplicated BAM
# to assess the impact of PCR duplication
# ============================================================

printf "Starting PreSeq...\n\n"
mkdir -p ${path2wd}/preseq

singularity exec --bind ${HOME},${path2genomes},${path2wd} \
    ${path2sif}/preseq_latest.sif preseq --version

for READ in 1 2; do
    singularity exec --bind ${HOME},${path2genomes},${path2wd} \
        ${path2sif}/preseq_latest.sif preseq lc_extrap \
        -output  ${path2wd}/preseq/${prefix}_r${READ}.ccurve.txt \
        -verbose -bam \
        ${path2wd}/star_align/${prefix}.genome_mapped_r${READ}.sorted.bam

    singularity exec --bind ${HOME},${path2genomes},${path2wd} \
        ${path2sif}/preseq_latest.sif preseq lc_extrap \
        -output  ${path2wd}/preseq/${prefix}_dedup_r${READ}.ccurve.txt \
        -verbose -bam \
        ${path2wd}/dedup/${prefix}.dedup_r${READ}.bam
done

printf "Done for PreSeq!\n\n"

# ============================================================
# RSeQC — read distribution and duplication
# read_distribution.py : fraction of reads in each genomic feature
#                        (CDS exon, UTR, intron, intergenic…)
# read_duplication.py  : duplication rate by read position and sequence
#                        (produces R plots automatically)
# Run on r2 only for duplication (informative strand)
# ============================================================

printf "Starting RSeQC...\n\n"
mkdir -p ${path2wd}/rseqc

for READ in 1 2; do
    singularity exec --bind ${HOME},${path2genomes},${path2wd} \
        ${path2sif}/rseqc_5.0.1.sif read_distribution.py \
        -i ${path2wd}/star_align/${prefix}.genome_mapped_r${READ}.sorted.bam \
        -r ${path2genomes}/gencode.v48.bed \
        > ${path2wd}/rseqc/${prefix}.read_distribution_r${READ}.txt

    singularity exec --bind ${HOME},${path2genomes},${path2wd} \
        ${path2sif}/rseqc_5.0.1.sif read_distribution.py \
        -i ${path2wd}/dedup/${prefix}.dedup_r${READ}.bam \
        -r ${path2genomes}/gencode.v48.bed \
        > ${path2wd}/rseqc/${prefix}.read_distribution_dedup_r${READ}.txt
done

# Read duplication plots on r2 (informative strand)
singularity exec --bind ${HOME},${path2genomes},${path2wd} \
    ${path2sif}/rseqc_5.0.1.sif read_duplication.py \
    -i ${path2wd}/star_align/${prefix}.genome_mapped_r2.sorted.bam \
    -o ${path2wd}/rseqc/${prefix}.read_duplication_r2
Rscript --vanilla ${path2wd}/rseqc/${prefix}.read_duplication_r2.DupRate_plot.r

singularity exec --bind ${HOME},${path2genomes},${path2wd} \
    ${path2sif}/rseqc_5.0.1.sif read_duplication.py \
    -i ${path2wd}/dedup/${prefix}.dedup_r2.bam \
    -o ${path2wd}/rseqc/${prefix}.read_duplication_dedup_r2
Rscript --vanilla ${path2wd}/rseqc/${prefix}.read_duplication_dedup_r2.DupRate_plot.r

printf "All done!\n\n"

sacct --format="JobId,JobName,NodeList,State,Elapsed,Timelimit,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveVMSize,ReqMem,Submit,Eligible" -j ${SLURM_JOB_ID}



