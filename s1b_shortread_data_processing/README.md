## Downloading and extracting shortread data from SRA

1. Download mouse EDL and soleus shortread data from SRA from this study
    - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6070147/
    
2. The study has illumina poly-A selected RNA-Seq data from soleus and EDL muscles
    - 2x101bp EDL wt (TC43508, EDL_Control replicate 1)
        - https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR7415846
        - https://www.ncbi.nlm.nih.gov/sra/SRX4286810[accn]
    - 2x101bp Soleus wt
        - https://www.ncbi.nlm.nih.gov/sra/SRX4286806[accn] 
        - https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR7415850 

3. On you HPC (SLURM), start an interactive session to download the data, and load the sratoolkit
```bash
sinteractive --mem=5g --cpus-per-task=8
module load sratoolkit
```

4. Check the size of the file you plan to download
```bash
## File 1 (TC43508, EDL_Control replicate 1)
vdb-dump --info SRR7415846

# Output
acc    : SRR7415846
path   : https://sra-download.ncbi.nlm.nih.gov/traces/sra65/SRR/007242/SRR7415846
size   : 17,530,199,642
type   : Table
platf  : SRA_PLATFORM_ILLUMINA
SEQ    : 168,602,712
SCHEMA : NCBI:SRA:Illumina:tbl:phred:v2#1.0.4
TIME   : 0x000000005b2e8253 (06/23/2018 13:24)
FMT    : Fastq
FMTVER : 2.8.2
LDR    : fastq-load.2.8.2
LDRVER : 2.8.2
LDRDATE: Mar  2 2017 (3/2/2017 0:0)

## File 2
vdb-dump --info SRR7415850

#Output
acc    : SRR7415850
path   : https://sra-download.ncbi.nlm.nih.gov/traces/sra65/SRR/007242/SRR741585
0
size   : 17,956,515,226
type   : Table
platf  : SRA_PLATFORM_ILLUMINA
SEQ    : 171,767,706
SCHEMA : NCBI:SRA:Illumina:tbl:phred:v2#1.0.4
TIME   : 0x000000005b2e7c59 (06/23/2018 12:59)
FMT    : Fastq
FMTVER : 2.8.2
LDR    : fastq-load.2.8.2
LDRVER : 2.8.2
LDRDATE: Mar  2 2017 (3/2/2017 0:0)

/scratch/$USER
```

5. Download the files using fasterq-dump and automatically convert to fastq files.
    - the `-t` is a temp folder for the files, but the current working dir is where the files will ultimately end up
```bash
fasterq-dump -p -t /scratch/$USER SRR7415846
fasterq-dump -p -t /scratch/$USER SRR7415850
```

6. Run the rest of the scripts in the directory to QC, trim, align and index the short-read data. Note, this only works if your HPC has the Slurm workload manager and the software modules installed.
	- Most scripts are run with `sbatch <script>_sbatch.sh` or `bash <script.swarm>`
```bash
# Step1: Initial FastQC
swarm -f s1_fastqc.swarm --module fastqc -g 4 -t 8

# Step 2: Trim adapters and low quality bases
swarm -f s2_trim.swarm --time 08:00:00 --module trimmomatic/0.36 -g 4 -t 8

# Step 3: Do another FastQC check to see if quality has improved
swarm -f s3_fastqc_again.swarm --module fastqc -g 4 -t 8

# Step 4: Generate the genome index files for STAR
sbatch s4_index_genome_sbatch.sh

# Step 5: Use STAR to align short-reads to genome
sbatch s5_edl_map_sbatch.sh
sbatch s5_soleus_map_sbatch.sh

# Step 6: Index the bam files (already pre-sorted using STAR flags)
sbatch s6_index_bams_sbatch.sh
```