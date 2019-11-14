## exCOVator - exon-based analysis

#### Prerequisites:
1. This is a Python 2.7 pipeline I wrote a few years ago. Its recommended that you use Anaconda to setup the enviornment.
	- Navigate to `https://www.anaconda.com/distribution/#download-section` and download the commmand line installer (NOT the GUI version) appropriate for your operating system.  I would get the Python 3 version. We can still setup a python2 enviornment as well as we will do below.

2. Install Anaconda following the command line installation for your operating system
	- https://docs.anaconda.com/anaconda/install/
	- e.g. (https://docs.anaconda.com/anaconda/install/mac-os/#using-the-command-line-install)


3. Once complete. Run the `enviornment.yaml` file to create a virtual environment in python 2 will all the required python packages to run the exCOVator pipeline.
	- On the command line run:
```bash
conda env create --file environment.yml
```

	- This should install the following packages:
		- python=2.7
		- biopython (conda)
		- pandas (conda)
		- tqdm (conda)
		- matplotlib (conda)
		- argparse (pip)
	  	- htseq (pip)
		- pyensembl (pip)
		- mygene (pip)

4. Now activate the enviornment when you are ready to run the pipeline
	- may take time to create the enviornment, it often seems like its stuck for 5 mins at `Executing transaction: done`, but keep waiting; its working under the hood.
```bash
conda activate excovator-env
```

5. I have included [dexseq_prepare_annotation.py](https://github.com/Bioconductor-mirror/DEXSeq/blob/release-3.4/inst/python_scripts/dexseq_prepare_annotation.py) script from the DEXSeq R package.  You can read more about it here: [DEXSeq R package from bioconductor](http://bioconductor.org/packages/release/bioc/html/DEXSeq.html).

6. Download the primary Gencode annotation for the species you are interested in:
	- e.g. I used: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.primary_assembly.annotation.gtf.gz
	- `gencode.vM10.primary_assembly.annotation.sorted.gtf --> mm10.gtf`

7. Download the bam files from GEO/SRA
	- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138362
		- cardiac_muscle_5-10kb_mm10_sorted.bam, cardiac_muscle_5-10kb_mm10_sorted.bam.bai
		- edl_muscle_5-10kb_mm10_sorted.bam, edl_muscle_5-10kb_mm10_sorted.bam.bai
		- soleus_muscle_5-10kb_mm10_sorted.bam, soleus_muscle_5-10kb_mm10_sorted.bam.bai

8. Rename or create softlinks to simplify
```bash
ln -s cardiac_muscle_5-10kb_mm10_sorted.bam cardiac.bam
ln -s cardiac_muscle_5-10kb_mm10_sorted.bam.bai cardiac.bam.bai
ln -s edl_muscle_5-10kb_mm10_sorted.bam edl.bam
ln -s edl_muscle_5-10kb_mm10_sorted.bam.bai edl.bam.bai
ln -s soleus_muscle_5-10kb_mm10_sorted.bam soleus.bam
ln -s soleus_muscle_5-10kb_mm10_sorted.bam.bai soleus.bam.bai
```

#### Step 1: Get isoform counts for each gene, based on gene_name

1. Use the uncollapsed bam files from the IsoSeq method

```bash
htseq-count -f bam -i gene_name ./bams/cardiac.bam ./refs/mm10.gtf > ./output/s1_cardiac_htseq-count.txt
htseq-count -f bam -i gene_name ./bams/edl.bam ./refs/mm10.gtf > ./output/s1_edl_htseq-count.txt
htseq-count -f bam -i gene_name ./bams/soleus.bam ./refs/mm10.gtf > ./output/s1_soleus_htseq-count.txt
```

#### Step 2: Get a list of genes with the most abundant isoforms
Without this selection step, it would take too long analyzing uncovered genes (probably days).

1. In copy each of the isoform counts for each sample and place into an Excel sheet according to output order. Move the last 5 rows of data to another sheet. (`__no_feature` etc)
	- `s2_combined_filtered_htseq-count.xlsx`
	- Be careful of the `March` genes which will be converted by Excel into dates. Highlight the gene column before pasting into Excel and change the column from "General" data to "Text" 
2. Sort each sample's isoforms in decending order 
3. Copy the names of the top expressed genes (e.g. >= 10 isoforms per gene etc) onto a new sheet and place in the same column
4. Sort the column by name in acending order
5. Filter for unique entries only
6. Remove all the Rik genes
7. Export the Excel sheet with the single column of gene names to a csv file (omit headers)
	- `s3_686_gene_list.csv`

#### Step 3: Select the genes from the original annotation gtf file
We grab the 686 genes and all their transcript annotations from the original reference gtf
- The reference gtf can be downloaded using the script in the `refs` folder.
```bash
bash ./refs/download_gencode_gtf.sh
```
- Or download it manually here: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.primary_assembly.annotation.gtf.gz
- You could also use the latest primary assembly if you would like
- `mm10.gtf -> gencode.vM10.primary_assembly.annotation.gtf`
```bash
python s3_selectGencodeAnnot.py -i ./refs/mm10.gtf -g s3_686_gene_list.csv -o ./refs/mm10_686.gtf
```

#### Step 4: Collapse gtf using the DEXSeq script
Collapse the transcripts down into a single meta transcript per gene.
```bash
python s4_dexseq_prepare_annotation.py ./refs/mm10_686.gtf ./refs/mm10_686.gff
```

#### Step 5: Get coverage data & plots
Determines the gene counts and percent spliced/fratios using the full-length read counts of all exons found with the collapsed annotation and unannotated exons found in the read files. Outputs a CSV file of the processed data and a PDF of gene graphs
```bash
python s5_exCOVator.py --bam ./bams/cardiac.bam ./bams/edl.bam ./bams/soleus.bam --gff ./refs/mm10_686.gff --species mouse --release 85 -o ./output/686_genes_exon_counts
```
- This step will take around 10-15 mins to run on a modern computer.
- Its not very efficient (I'm an amateur a developer)

#### Step 6: Filter for differentially used exons
Script will output a csv with the filtered list of exons and a simple list of their gene symbols for step 7
- Filtered by 10% difference: took 686 genes --> 443 genes differentially used, 2630 exons
```bash
python s6_filterDiffUsedExons.py -i ./output/686_genes_exon_ratios.csv -o ./output/filtered_10p_diff_used --ratio_diff 0.10
```

#### Step 7: Gene summaries
Grabs full gene names and their summaries from NCBI to simplify analysis
- provide gene list output from step 6
- need to specify animal for grabbing the summaries from. Human will have the most gene summaries compared to mouse, but will not have mouse specific genes. Can run again with diff animal 
- NCBI requires an email to contact you if you overload their servers with too many requests over a short period of time. Should not be an issue here.
```bash
python s7_getGeneSummaries.py -g ./output/filtered_10p_diff_used_gene_list.csv -a human -o ./output/443_diff_used -e your.email@institution.edu
```


