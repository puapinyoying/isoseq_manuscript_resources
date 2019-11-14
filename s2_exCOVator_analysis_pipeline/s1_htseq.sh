#! /bin/bash

# Activate enviornment before running script
# conda activate excovator-env

htseq-count -f bam -i gene_name ./bams/cardiac.bam ./refs/mm10.gtf > ./output/cardiac_htseq-count.txt
htseq-count -f bam -i gene_name ./bams/edl.bam ./refs/mm10.gtf > ./output/edl_htseq-count.txt
htseq-count -f bam -i gene_name ./bams/soleus.bam ./refs/mm10.gtf > ./output/soleus_htseq-count.txt