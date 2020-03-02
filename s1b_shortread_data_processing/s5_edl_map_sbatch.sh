#! /bin/bash
# this file is s5_edl_sbatch.sh
set -o pipefail
set -e

function fail {
    echo "$@" >&2
    exit 1
}

module load samtools/1.6 || fail "could not load samtools module"
module load STAR/2.5.2a        || fail "could not load STAR module"

GENOME=/isoseq_project/references
FASTQ=/isoseq_project/shortread/results/s2_trimmed
OUTPUT=/isoseq_project/shortread/results/s5_mapped

cd $OUTPUT || fail "no such directory"

STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $GENOME \
    --sjdbOverhang 100 \
    --readFilesIn $FASTQ/EDL_CTR_SRR7415846_1_paired.fastq $FASTQ/EDL_CTR_SRR7415846_2_paired.fastq \
    --outSAMtype BAM SortedByCoordinate \
    --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp \
    --outFileNamePrefix $OUTPUT/EDL_CTR_SRR7415846_sorted