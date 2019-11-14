exPhaser protocol
Brian Uapinyoying

#### Prerequisites:
1. This is a Python >= 3.6 pipeline. Its recommended that you use Anaconda to setup the enviornment.
	- Navigate to `https://www.anaconda.com/distribution/#download-section` and download the commmand line installer (NOT the GUI version) appropriate for your operating system. Get the Python 3 version. 

2. Install Anaconda following the command line installation for your operating system
	- https://docs.anaconda.com/anaconda/install/
	- e.g. (https://docs.anaconda.com/anaconda/install/mac-os/#using-the-command-line-install)

3. Once complete. Run the `enviornment.yaml` file to create a virtual environment in python 2 will all the required python packages to run the exCOVator pipeline.
	- On the command line run:
```bash
conda env create --file environment.yml
```
	- This should install the following packages:
		- python=3.6
		- pandas (conda)
		- tqdm (conda)
		- matplotlib (conda)
		- argparse (pip)
	  	- htseq (pip)

4. Make some output directories
```bash
# from root folder
mkdir output output/neb output/nrap output/ttn
```

5. Download the gencode gtf annotation file to `refs` directory
The script will download and create a softlink to the file called mm10.gtf
- original reference used in publication: `gencode.vM10.primary_assembly.annotation.gtf`
```bash
cd ./refs
bash download_gencode_gtf.sh
```

6. Download the bam files from GEO/SRA and place in `bams` directory
	- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138362
		- cardiac.bam, cardiac.bam.bai
		- edl.bam, edl.bam.bai
		- soleus.bam, soleus.bam.bai

7. Rename or create softlinks to simplify
```bash
# from bam directory
ln -s cardiac.bam cardiac.bam
ln -s cardiac.bam.bai cardiac.bam.bai
ln -s edl.bam edl.bam
ln -s edl.bam.bai edl.bam.bai
ln -s soleus.bam soleus.bam
ln -s soleus.bam.bai soleus.bam.bai
```

#### Step 1. Subset the GTF files for genes of interest to speed up processing 
Do this by using grep and adding the ensembl gene id to the query.
- Nebulin	= ENSMUSG00000026950
- Nrap 		= ENSMUSG00000049134
- Titin 	= ENSMUSG00000051747

```bash
cd ./refs
# in refs directory

grep "ENSMUSG00000026950" mm10.gtf > neb_mm10.gtf
grep "ENSMUSG00000049134" mm10.gtf > nrap_mm10.gtf
grep "ENSMUSG00000051747" mm10.gtf > ttn_mm10.gtf
```

#### Step 2. Phase exons
Use included files in `bed_files` directory to recreate analysis in publication

1. Nrap exons 2-40
```bash
cd .. # Back to root folder
python exPhaser.py \
--bams bams/cardiac.bam bams/edl.bam bams/soleus.bam \
--gtf ./refs/nrap_mm10.gtf \
--bed ./bed_files/nrap_2-40.bed \
--outputPrefix ./output/nrap/nrap_2-40
```

2. Nebulin exons 122-152 (these are a lot of exons, it will take very long about 30 mins and seem like it hangs when checking transcript patterns, but I just forgot to add another progress bar :/)
```bash
python exPhaser.py \
--bams bams/cardiac.bam bams/edl.bam bams/soleus.bam \
--gtf ./refs/neb_mm10.gtf \
--bed ./bed_files/neb_122-152.bed \
--outputPrefix ./output/neb/neb_122-152
```

3. Titin exons 45-50 and 45-169 (skeletal vs cardiac Titin isoforms)
```bash
python exPhaser.py \
--bams bams/cardiac.bam bams/edl.bam bams/soleus.bam \
--gtf ./refs/ttn_mm10.gtf \
--bed ./bed_files/ttn_45-50.bed \
--outputPrefix ./output/ttn/ttn_45-50.bed

python exPhaser.py \
--bams bams/cardiac.bam bams/edl.bam bams/soleus.bam \
--gtf ./refs/ttn_mm10.gtf \
--bed ./bed_files/ttn_45-169.bed \
--outputPrefix ./output/ttn/ttn_45-169.bed
```

4. Titin exons 10-14 (targeting 11, 12 & 13)
```bash
python exPhaser.py \
--bams bams/cardiac.bam bams/edl.bam bams/soleus.bam \
--gtf ./refs/ttn_mm10.gtf \
--bed ./bed_files/ttn_10-14.bed \
--outputPrefix ./output/ttn/ttn_10-14.bed
```

5. Titin exons 190-192 (targeting 191)
```bash
python exPhaser.py \
--bams bams/cardiac.bam bams/edl.bam bams/soleus.bam \
--gtf ./refs/ttn_mm10.gtf \
--bed ./bed_files/ttn_190-192.bed \
--outputPrefix ./output/ttn/ttn_190-192.bed
```

6. Titin exon 311-313 (targeting 312)
```bash
python exPhaser.py \
--bams bams/cardiac.bam bams/edl.bam bams/soleus.bam \
--gtf ./refs/ttn_mm10.gtf \
--bed ./bed_files/ttn_311-313.bed \
--outputPrefix ./output/ttn/ttn_311-313.bed
```

7. Titin exon 44 and alternate 3' exon 45 (targeting exon alt3' 45)
```bash
python exPhaser.py \
--bams bams/cardiac.bam bams/edl.bam bams/soleus.bam \
--gtf ./refs/ttn_mm10.gtf \
--bed ./bed_files/ttn_alt3prime45-44.bed \
--outputPrefix ./output/ttn/ttn_alt3prime45-44.bed
```