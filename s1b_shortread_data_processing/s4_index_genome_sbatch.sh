#! /bin/bash
# this file is s4_index_genome_sbatch.sh
set -o pipefail
set -e

function fail {
    echo "$@" >&2
    exit 1
}

module load samtools/1.6 || fail "could not load samtools module"
module load STAR/2.5.2a        || fail "could not load STAR module"

GENOME=/isoseq_project/references

cd $GENOME || fail "no such directory"

STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --runMode genomeGenerate \
    --genomeDir $GENOME \
    --genomeFastaFiles $GENOME/GRCm38.primary_assembly.genome.fa \
    --sjdbGTFfile $GENOME/gencode.vM10.primary_assembly.annotation.gtf \
    --sjdbOverhang 100 \
    --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp