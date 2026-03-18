#!/bin/bash
#
#SBATCH --job-name=umi2dedup
#SBATCH --partition=all_5hrs
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=32G
#SBATCH --ntasks=1
#SBATCH --output=/path/to/your/eprint/Logs/umi2dedup/umi2dedup.%j.log

printf " ################################################################## \n # Starting UMI2DEDUP script for $prefix # \n ################################################################## \n\n"

source ~/.bashrc
module load cutadapt/4.4-GCCcore-12.2.0
module load SAMtools/1.17-GCC-12.2.0

# ============================================================
# UMI extraction
# Moves the 10nt UMI from the 5' end of each read into the
# read name, where umi_tools dedup can later use it
# ============================================================

# Python virtual environment containing umi_tools
# Update this path to your environment location
VENV_PATH=/path/to/your/eprint/MyPythonEnv

python3 -m venv ${VENV_PATH}
source ${VENV_PATH}/bin/activate

mkdir -p ${path2wd}/umi_extract
cd ${path2wd}

printf "Launching UMI-EXTRACT...\n"
umi_tools --version

umi_tools extract --random-seed 1 --bc-pattern NNNNNNNNNN \
    --stdin  ${path2fastq}/${prefix}_R1_001.fastq.gz \
    --stdout ${path2wd}/umi_extract/${prefix}_1.umi.fastq.gz \
    --log    ${path2wd}/umi_extract/${prefix}_1.umi.log

umi_tools extract --random-seed 1 --bc-pattern NNNNNNNNNN \
    --stdin  ${path2fastq}/${prefix}_R2_001.fastq.gz \
    --stdout ${path2wd}/umi_extract/${prefix}_2.umi.fastq.gz \
    --log    ${path2wd}/umi_extract/${prefix}_2.umi.log

deactivate
printf "\nUMI-EXTRACT done!\n\n"

# ============================================================
# Adapter trimming with cutadapt
# The adapter list covers the Illumina TruSeq sequence and its
# partial suffixes, ensuring truncated adapters are also removed
# --quality-cutoff 6 : trim low-quality 3' bases before adapter removal
# -m 18              : discard reads shorter than 18nt after trimming
# ============================================================

mkdir -p ${path2wd}/cutadapt
printf "cutadapt version: "; cutadapt --version
printf "Starting trimming...\n\n"

ADAPTERS="\
    -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -a CTTCCGATCTACAAGTT \
    -a CTTCCGATCTTGGTCCT \
    -a AACTTGTAGATCGGA \
    -a AGGACCAAGATCGGA \
    -a ACTTGTAGATCGGAA \
    -a GGACCAAGATCGGAA \
    -a TTGTAGATCGGAAGA \
    -a ACCAAGATCGGAAGA \
    -a TGTAGATCGGAAGAG \
    -a CCAAGATCGGAAGAG \
    -a GTAGATCGGAAGAGC \
    -a CAAGATCGGAAGAGC \
    -a TAGATCGGAAGAGCG \
    -a AAGATCGGAAGAGCG \
    -a AGATCGGAAGAGCGT \
    -a GATCGGAAGAGCGTC \
    -a ATCGGAAGAGCGTCG \
    -a TCGGAAGAGCGTCGT \
    -a CGGAAGAGCGTCGTG \
    -a GGAAGAGCGTCGTGT"

cutadapt --match-read-wildcards --times 1 --cores 0 \
    -e 0.1 -O 1 --quality-cutoff 6 -m 18 \
    ${ADAPTERS} \
    -o ${path2wd}/cutadapt/${prefix}_1.trim.fastq.gz \
    ${path2wd}/umi_extract/${prefix}_1.umi.fastq.gz \
    > ${path2wd}/cutadapt/${prefix}_1.trim.log

cutadapt --match-read-wildcards --times 1 --cores 0 \
    -e 0.1 -O 1 --quality-cutoff 6 -m 18 \
    ${ADAPTERS} \
    -o ${path2wd}/cutadapt/${prefix}_2.trim.fastq.gz \
    ${path2wd}/umi_extract/${prefix}_2.umi.fastq.gz \
    > ${path2wd}/cutadapt/${prefix}_2.trim.log

printf "Done trimming!\n\n"

# ============================================================
# FASTQ sorting by read name
# Required so that R1 and R2 reads are paired in the same order
# before alignment. Uses fastq-sort (update path below).
# ============================================================

# Update this path to your fastq-sort binary
FASTQSORT=/path/to/your/bin/fastq-sort

printf "Sorting FASTQ files...\n\n"

for READ in 1 2; do
    printf "Processing read %s...\n" "${READ}"
    gunzip -c ${path2wd}/cutadapt/${prefix}_${READ}.trim.fastq.gz \
        > ${path2wd}/cutadapt/${prefix}_${READ}.trim.fastq
    ${FASTQSORT} --id ${path2wd}/cutadapt/${prefix}_${READ}.trim.fastq \
        > ${path2wd}/cutadapt/${prefix}_${READ}.trim.sorted.fastq
    gzip -f ${path2wd}/cutadapt/${prefix}_${READ}.trim.sorted.fastq
    rm   ${path2wd}/cutadapt/${prefix}_${READ}.trim.fastq
    printf "Done for read %s\n\n" "${READ}"
done

# ============================================================
# Alignment with STAR — step 1: RepBase (repeat elements)
# Reads mapping to repeat elements are discarded; unmapped reads
# are carried forward to the genome alignment step
# --outFilterMultimapNmax 30 : allow up to 30 multimappers for RepBase
# ============================================================

path2out=${path2wd}/star_align
mkdir -p ${path2out}

printf "STAR version: "
singularity exec --bind ${path2genomes},${path2wd},${path2fastq},${path2out} \
    ${path2sif}/star_2.7.11b.sif STAR --version

STAR_COMMON_ARGS="\
    --alignEndsType EndToEnd \
    --genomeLoad NoSharedMemory \
    --outBAMcompression 10 \
    --outFilterMultimapScoreRange 1 \
    --outFilterScoreMin 10 \
    --outFilterType BySJout \
    --outReadsUnmapped Fastx \
    --outSAMattrRGline ID:foo \
    --outSAMattributes All \
    --outSAMmode Full \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --outStd Log \
    --runMode alignReads \
    --runThreadN 8"

for READ in 1 2; do
    printf "\nAligning read %s to RepBase...\n\n" "${READ}"
    singularity exec --bind ${path2genomes},${path2wd},${path2fastq},${path2out} \
        ${path2sif}/star_2.7.11b.sif STAR \
        ${STAR_COMMON_ARGS} \
        --genomeDir             ${path2genomes}/repbase_STARindex/ \
        --outFilterMultimapNmax 30 \
        --outFileNamePrefix     ${path2out}/${prefix}.mapped_repbase_r${READ}. \
        --readFilesCommand      zcat \
        --readFilesIn           ${path2wd}/cutadapt/${prefix}_${READ}.trim.sorted.fastq.gz
    printf "Mapped to RepBase (r%s): " "${READ}"
    samtools view -c -F 260 ${path2out}/${prefix}.mapped_repbase_r${READ}.Aligned.out.bam
done

# ============================================================
# Alignment with STAR — step 2: genome (GRCh38 / GENCODE)
# Only reads that did NOT map to RepBase are used as input
# --outFilterMultimapNmax 1 : keep only uniquely mapping reads
# ============================================================

for READ in 1 2; do
    printf "\nAligning read %s to genome...\n\n" "${READ}"
    singularity exec --bind ${path2genomes},${path2wd},${path2fastq},${path2out} \
        ${path2sif}/star_2.7.11b.sif STAR \
        ${STAR_COMMON_ARGS} \
        --genomeDir             ${path2genomes}/gencode_star_index \
        --outFilterMultimapNmax 1 \
        --outFileNamePrefix     ${path2out}/${prefix}.genome_mapped_r${READ}. \
        --readFilesIn           ${path2out}/${prefix}.mapped_repbase_r${READ}.Unmapped.out.mate1
    printf "Mapped to genome (r%s): " "${READ}"
    samtools view -c -F 260 ${path2out}/${prefix}.genome_mapped_r${READ}.Aligned.out.bam
done

# ============================================================
# BAM sorting and indexing
# Name-sort first (required by umi_tools dedup),
# then coordinate-sort and index for downstream tools
# ============================================================

printf "\nSorting and indexing BAM files...\n\n"
samtools --version | head -n 1

for READ in 1 2; do
    samtools sort -n \
        -o ${path2out}/${prefix}.genome_mapped_r${READ}.namesorted.bam \
           ${path2out}/${prefix}.genome_mapped_r${READ}.Aligned.out.bam
    samtools sort \
        -o ${path2out}/${prefix}.genome_mapped_r${READ}.sorted.bam \
           ${path2out}/${prefix}.genome_mapped_r${READ}.namesorted.bam
    samtools index ${path2out}/${prefix}.genome_mapped_r${READ}.sorted.bam
    rm ${path2out}/${prefix}.genome_mapped_r${READ}.namesorted.bam
done

# ============================================================
# PCR duplicate removal with umi_tools dedup
# --method unique      : only deduplicate reads with identical UMI + position
# --spliced-is-unique  : treat spliced reads as distinct from unspliced
# ============================================================

printf "\nStarting deduplication...\n\n"

source ${VENV_PATH}/bin/activate
mkdir -p ${path2wd}/dedup

for READ in 1 2; do
    umi_tools dedup \
        --method unique \
        --spliced-is-unique \
        --log    ${path2wd}/dedup/${prefix}.process_r${READ}.log \
        -I       ${path2out}/${prefix}.genome_mapped_r${READ}.sorted.bam \
        --output-stats ${path2wd}/dedup/${prefix}.umi_dedup_stats_r${READ} \
        -S       ${path2wd}/dedup/${prefix}.dedup_r${READ}.bam
    printf "Deduplicated reads (r%s): " "${READ}"
    samtools view -c -F 260 ${path2wd}/dedup/${prefix}.dedup_r${READ}.bam
    samtools index ${path2wd}/dedup/${prefix}.dedup_r${READ}.bam
done

deactivate
printf "\nAll done!\n\n"

sacct --format="JobId,JobName,NodeList,State,Elapsed,Timelimit,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveVMSize,ReqMem,Submit,Eligible" -j ${SLURM_JOB_ID}


