#! /bin/bash
# this file is s6_index_bams_sbatch.sh
set -o pipefail
set -e

function fail {
    echo "$@" >&2
    exit 1
}

module load samtools/1.6 || fail "could not load samtools module"

BAMS=/isoseq_project/shortread/results/s5_mapped

cd $BAMS || fail "no such directory"

samtools index EDL_CTR_SRR7415846_sorted.bam
samtools index SOLEUS_CTR_SRR7415850_sorted.bam