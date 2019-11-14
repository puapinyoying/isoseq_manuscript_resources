# IsoSeq bioinformatics workflow:
-------------------------------

Brian Uapinyoying
Updated: 10/02/2019

---

### Index

*   [Part I: Introduction and setup](#part-i-introduction-and-setup)
    *   [Introduction](#introduction)
    *   [Installation Requirements](#installation-requirements)
    *   [Installing pbtranscript-TOFU](#installing-pbtranscript-tofu)
*   [Part II: Obtaining reads of insert](#part-ii-obtaining-reads-of-insert)
*   [Part III: Classifying the reads of insert](#part-iii-classifying-the-reads-of-insert)
*   [Part IIIb: Demultiplexing the samples (OPTIONAL)](#part-iiib-demultiplexing-the-samples)
*   [Part IV: Cluster (ICE & QUIVER)](#part-iv-cluster)
*   [Part V: Setup and alignment with STARlong](#part-v-setup-and-alignment-with-starlong)
*   [Part VI: File conversions with SamTools](#part-vi-file-conversions-with-samtools)
*   Optional steps:
    * ***NOTE***: IsoSeq analysis using exCOVator.py and exPhaser made use of uncollapsed data. We found that collapsing the isoforms results in loss of quantitative ability. If following protocol in our publication, skip to Part X.
    *   [Part VII: Collapse Redundant Isoforms](part-vii-collapse-redundant-isoforms)
    *   [Part VIII: STARlong and Samtools revisited](#part-viii-starlong-and-samtools-revisited)
    *   [Part IX: MatchAnnot](#part-ix-matchannot)
*   [Part X: HTSeq-count](#part-x-htseq-count)
*   [Part XI: Visualizing Data](#part-xi-visualizing-data)

### Part I: Introduction and setup
---

#### Introduction
This protocol is my personal/streamlined version of [Pacbio's IsoSeq cDNA_primer tutorials.](https://github.com/PacificBiosciences/cDNA_primer/wiki#rs_isoseq-official-pipeline-tutorial) I use all the outlined tools provided by PacBio except for the sample demultiplexing step, which they have not provided an official solution to yet (at time of writing). So I created a custom python script to do the job. The protocol also assumes that you would like to use STAR to do genome alignments, instead of PacBio's native GMAP, in order to get the best possible results. PacBio has done a comparison between the two aligners and found that STAR performs better (fewer artifacts/false positives) than GMAP. The study and run parameters are detailed [here.](https://github.com/PacificBiosciences/cDNA_primer/wiki/Bioinfx-study%3A-Optimizing-STAR-aligner-for-Iso-Seq-data) However, if you plan to use the native GMAP, there are benefits. The developer's version of the isoseq protocol ([pbtranscript-tofu](https://github.com/PacificBiosciences/cDNA_primer/wiki/What-is-pbtranscript-tofu%3F-Do-I-need-it%3F)) has a wrapper script (tofu_wrap.py) that can streamline this entire protocol in a few short commands. I will not be covering these in this protocol so please read their [official notes.](https://github.com/PacificBiosciences/cDNA_primer/wiki/tofu-Tutorial-%232.-Isoform-level-clustering-%28ICE-and-Quiver%29)

I am also using a Debian based Linux enviornment, so all instructions will reflect that.  It's not too much different to use RedHat Linux or MacOS X (Unix).


#### Installation Requirements
---

##### 1.  **PacBio SmrtAnalysis Version 2.3.0** (as of 02/07/2016)
*   [Install Documentation](http://www.pacb.com/wp-content/uploads/SMRT-Analysis-Software-Installation-v2-3-0.pdf)
*   [Website](http://www.pacb.com/support/software-downloads/)
*   [Download link](http://files.pacb.com/software/smrtanalysis/2.3.0/smrtanalysis_2.3.0.140936.run)
    
##### 2.  **SmrtAnalysis Update Patch 5** (latest patch as of 02/07/2016)
*   [Download link](https://s3.amazonaws.com/files.pacb.com/software/smrtanalysis/2.3.0/smrtanalysis-patch_2.3.0.140936.p5.run)  
*   You can install smrtAnalysis and its latest patch in one [command](http://www.pacb.com/wp-content/uploads/SMRT-Analysis-Software-Installation-v2-3-0.pdf)
*   `bash smrtanalysis_2.3.0.140936.run -p smrtanalysis-patch_2.3.0.140936.p5.run --rootdir`

##### 3.  **Git/Github** access (many tools are hosted on github, it's helpful to )
*   `sudo apt-get install git`
*   OPTIONAL:       Setup your own [github account](https://help.github.com/articles/set-up-git/#platform-linux)

##### 4.  [PBtranscript-TOFU](https://github.com/PacificBiosciences/cDNA_primer/wiki/Setting-up-virtualenv-and-installing-pbtranscript-tofu) - Developmental version of smrtAnalysis
*   Installing this software can be tricky. It is a little buggy and requires very specific versions of dependencies, so I am going to go over this in more detail. See the [next section.](#installing-pbtranscript-tofu)
*   This is software is optional if following the publication protocol.
                        
##### 5.  **Python 2.7**
*   [Anaconda](http://docs.continuum.io/anaconda/install) - Simplifies downloading and installing scientific python packages.  Recommended. Install locally (not to smrtanalysis)         

*   [virtualenv](https://virtualenv.readthedocs.org/en/latest/) - Required for pbtranscript-TOFU installation and use. Allows you to run a virtual instance of python with custom settings without causing problems with original installation.   
    *   `pip install virtualenv`

*   [BioPython](https://github.com/biopython/biopython) - For reading fasta files and much more. Required for running custom demultiplexing script
    *   `pip install Bio` (or `conda install Bio`)

*   [tqdm](https://github.com/tqdm/tqdm) - Easy progress bar library, used in custom demultiplexingscript.
    *   `pip install tqdm` 

*   [HTSeq](http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html) - Counting gene coverage. Useful for quantifying isoforms per gene
    *   [Installation instructions](http://www-huber.embl.de/users/anders/HTSeq/doc/install.html)

##### 6.  [Samtools](http://www.htslib.org/) - 
*   For manipulating and converting SAM and BAM alignment files  

##### 7.  [STAR aligner](https://github.com/alexdobin/STAR) (get version => 2.5)
*   Download:    `wget https://github.com/alexdobin/STAR/archive/STAR_2.5.0a.tar.gz`
*   [Documentation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

##### 8.  **Primary assembly genomes and annotation (GTF)**
*   Genomes should have all chromosomes and scaffolds, but not haplotypes.
*   Generally called the primary assemblies on ensembl, ucsc and gencode.
*   The toplevel assemblies include haplotype information.
*   Grab the latest version for your organism (if available)
*   Look to [ensembl](http://useast.ensembl.org/info/data/ftp/index.html) or [UCSC](https://genome.ucsc.edu/) for references for other organisms  
*   [Gencode project](http://www.gencodegenes.org/about.html) has comprehensive annotations for human and mouse genomes (recommended). For latest mouse as of (02/10/2016). Be sure to grab the matching primary comprehensive gene annotation files (GTF)

```
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M8/GRCm38.primary_assembly.genome.fa.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M8/gencode.vM8.primary_assembly.annotation.gtf.gz
```

##### 9.  [MatchAnnot software](https://github.com/TomSkelly/MatchAnnot)
*   For comparing splice isoforms in our results to known gene annotations. Comes with ClusterView. Both are python scripts and easy to install.  

##### 10.  [Integrated Genome Visualizer (IGV)](https://www.broadinstitute.org/igv/) 
*   An easy to use, java-based genome browser


####  Installing pbtranscript-TOFU
---
[Github page](https://github.com/PacificBiosciences/cDNA_primer)

It's best to use the install script outlined in the [instructions](https://github.com/PacificBiosciences/cDNA_primer/wiki/Setting-up-virtualenv-and-installing-pbtranscript-tofu). If you try to use newer versions of the tools outlined in the instructions such as virtualenv, setuptools or bxpython it will likely not work and you won't know why.  I learned the hard way.

The provided script on the instruction page works best. However I had to tweak the script a bit to get it to work on my Ubuntu x64 14.04 LTS setup. 

##### 1.  **Create new folder in home directory**. Place the shell script below in the directory. 

Load the virtual enviornment and smrtShell before running the script below
    *   `/opt/smrtanalysis/current/etc/setup.sh`
    *   `/opt/smrtanalysis/smrtcmds/bin/smrtshell`

**Note**: the header should be `#!</path/to/smrtanalysis>/smrtcmds/bin/smrtshell`
```bash
!#/opt/smrtanalysis/smrtcmds/bin/smrtshell
set -vex
#module load ccache
#module load zlib
export VENV_TOFU=${PWD}/VENV_TOFU
wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz
tar zxf virtualenv-1.11.6.tar.gz -C /tmp/
python /tmp/virtualenv-1.11.6/virtualenv.py --system-site-packages $VENV_TOFU
source ${VENV_TOFU}/bin/activate

pip install urllib3 # added to script
pip install hmac # added to script

pip install --upgrade pip
pip install --upgrade setuptools
pip install Cython
pip install numpy
pip install bx-python
pip install pysam
pip install pbcore
```

##### 2.  **Run the script**
*   `bash <script_name.sh>`

##### 3.  **Exit smrtshell and then run git clone to download install files**
*   `git clone https://github.com/PacificBiosciences/cDNA_primer.git`

##### 4.  **Go back into the VENV_TOFU and smrtshell**
*   `source ${VENV_TOFU}/bin/activate`
*   `/opt/smrtanalysis/smrtcmds/bin/smrtshell`


##### 5. **Create a second shell script with below in it and run**

```bash
set -vex
export VENV_TOFU=${PWD}/VENV_TOFU
source ${VENV_TOFU}/bin/activate  

cd cDNA_primer/pbtranscript-tofu/pbtranscript/pbtools/pbtranscript/branch/C/
python setup2.py build_ext --inplace
cd ../../../../
make

cd ../external_daligner/
unzip  DALIGNER-d4aa4871122b35ac92e2cc13d9b1b1e9c5b5dc5c-ICEmod.zip
cd  DALIGNER-d4aa4871122b35ac92e2cc13d9b1b1e9c5b5dc5c-ICEmod/
make
cp HPCdaligner HPCmapper LA4Ice LAcat LAcheck LAmerge LAshow LAsort LAsplit daligner $VENV_TOFU/bin

cd ../
unzip DAZZ_DB-40bb7e4b2041cdfd3c50b22201301b8df06342fa.zip
cd  DAZZ_DB-40bb7e4b2041cdfd3c50b22201301b8df06342fa/
make
cp Catrack DAM2fasta DB2fasta DB2quiva DBdust DBrm DBshow DBsplit DBstats fasta2DAM fasta2DB quiva2DB simulator ${VENV_TOFU}/bin
```


##### 6.  Test install and troubleshooting. There is a [known issue](https://github.com/PacificBiosciences/cDNA_primer/issues/3) post-installation

1.  After TOFU installation, load the virtual enviornment and smrtShell
    *   `source ${VENV_ENV}/bin/activate`
    *   `/opt/smrtanalysis/smrtcmds/bin/smrtshell`

2.  Test if it works by trying:
    *   `collapse_isoforms_by_sam.py --help`

3.  If you get this error: 'ImportError: No module named modified_bx_intervals.intersection_unique'

4.  Go into the modified_bx_intervals directory
    *   `cd ${VENV_ENV}/lib/python2.7/site-packages/pbtools.pbtranscript-0.4.0-py2.7-linux-x86_64.egg/pbtools/pbtranscript/modified_bx_intervals/`

5.  See if you have these two files:
    * 'intersection_unique.so'
    * '__init__.py'

6.  If you are only missing the '__init__.py', just create it in the folder using the touch command
    *   `touch __init__.py`

7.  Rerun the command with help to see if the error is gone

8.  Worst case senario, you are missing 'intersection_unique.so'. Youwill need to reinstall pbtranscript-TOFU

Other installation issues see [github issues for cDNA_primer](https://github.com/PacificBiosciences/cDNA_primer/issues/17)


### Part II:  Obtaining reads of insert
---

[Official Documentation](https://github.com/PacificBiosciences/SMRT-Analysis/wiki/ConsensusTools-v2.3.0-Documentation#CCS_CL)

##### 1.  Create a file of file names (fofn) for the bax.h5 raw sequence files planned for processing. 

*   List the full path to each file. Example data from a single smrtCell:
```
/data/2015-11-02_M2_T1/E01_1/Analysis_Results/m151103_154521_42136_c100889062550000001823203204021684_s1_p0.1.bax.h5
/data/2015-11-02_M2_T1/E01_1/Analysis_Results/m151103_154521_42136_c100889062550000001823203204021684_s1_p0.2.bax.h5
/data/2015-11-02_M2_T1/E01_1/Analysis_Results/m151103_154521_42136_c100889062550000001823203204021684_s1_p0.3.bax.h5
```

##### 2.	**The command parameters**

```bash
ConsensusTools.sh CircularConsensus \
--minFullPasses 0  \
--minPredictedAccuracy 75 \
--parameters <parameters_directory>  # see notes below
--numThreads 24 --fofn /MYHOME/test_dir/input.fofn \
-o /MYHOME/test_dir/data
```

**Notes:**

*   For IsoSeq set minimum full passes to 0. We want to use every subread not just CCS reads. Long transcripts may only be able to get 1 pass or no full passes at all. 
*   Algorythm parameters are located in the directory below.  
*   `/<smrtanalysis_directory>/current/analysis/etc/algorithm_parameters/<YYYY-DD>/`
    -   `<YYYY-DD>` is the latest protocol version (as of writing it is '2015-11')

**Output:**

*   Each 'bax.h5' file will output: one 'fast**a**', one 'fast**q**'' and one 'ccs.h5' file 

##### 3.  Start smrt enviornment and shell
*    `source /opt/smrtanalysis/current/etc/setup.sh`
*    `/opt/smrtanalysis/smrtcmds/bin/smrtshell`

##### 4.  Run the command

**Example**:
```bash
ConsensusTools.sh CircularConsensus  --minFullPasses 0  \
--minPredictedAccuracy 75 \
--parameters /opt/smrtanalysis/current/analysis/etc/algorithm_parameters/2015-11/ \
--numThreads 15 \
--fofn LT1_input.fofn -o /data/brian/test/consensus/
```

**Test run start 01/29/2015 @ 1:30pm**

*   First bax.h5 file processed in 1:06 hours with 15 cores and 64gb ram using Ubuntu linux on a VM and data on an NFS drive
*   Each cell has 3 bax.h5 files so ~ 3 hours per cell.
*   for 24 smrtCells * 3 = 72 hours on a single node

##### 5.  **reads_of_insert**
If you had multiple bax.h5 files in your run, concatinate the output into reads_of_insert files for the next part.

```bash
cat /data/brian/test/consensus/*.ccs.fasta > reads_of_insert.fasta
cat /data/brian/test/consensus/*.ccs.fastq > reads_of_insert.fastq
find $PWD -name *.ccs.h5  > reads_of_insert.fofn
```

***WARNING*** The reads_of_insert.fofn is actually the fullpath *names* of the *.ccs.h5* files you are grouping to form the reads_of_insert. In the official tutorial, it makes it look like you are concatinating all of the data in the ccs.h5 files into a single file. That is incorrect.

### Part III:  Classifying the reads of insert
---

[cDNA_primer Github Wiki](https://github.com/PacificBiosciences/cDNA_primer/wiki/RS_IsoSeq-%28v2.3%29-Tutorial-%231.-Getting-full-length-reads) for official instructions.

##### 0.  **Demultiplexing**
If you multiplexed/barcoded your samples during library prep you must create a custom primers fasta file containing your barcode sequences. This will be key for demultiplexing your samples BEFORE the going into the cluster step. Otherwise, there are no subsequent ways to determine which reads came from which sample.

**Example**: Sequence header from an isoseq_flnc.fasta (classify output) file

```
>m140121_100730_42141_c100626750070000001823119808061462_s1_p0/119/30_1067_CCS strand=+;fiveseen=1;polyAseen=1;threeseen=1;fiveend=30;polyAend=1067;threeend=1096;primer=1;chimera=0
```
    
In the above header, the `primer=<int>;` field will only exist if you provide `pbtranscript.py classify` the `-p <custom_primers.fasta>` flag.

Here are all 6 possible (as of 02-02-2016) IsoSeq barcode sequences in order as they should appear in the file. Samples are numbered starting with zero first (e.g. sample 1 = F0, R0; sample 2 = F1, R1 ...)

***WARNING***: 

*   Do not rename the sequence headers (>F0, >R0, ...) or it will break! 
*   Delete any barcodes/primer pairs not used in your experiment or you may get reads misidentified to non-existant samples due to read errors.

**custom_primers.fasta**

```
>F0
AAGCAGTGGTATCAACGCAGAGTAC
>R0
AAATGACGCATCGTCTGAGTACTCTGCGTTGATACCACTGCTT
>F1
AAGCAGTGGTATCAACGCAGAGTAC
>R1
AAGCAGAGTCATGTATAGGTACTCTGCGTTGATACCACTGCTT
>F2
AAGCAGTGGTATCAACGCAGAGTAC
>R2
AAGAGTGCTACTCTAGTAGTACTCTGCGTTGATACCACTGCTT
>F3
AAGCAGTGGTATCAACGCAGAGTAC
>R3
AACATGTACTGATACACAGTACTCTGCGTTGATACCACTGCTT
>F4
AAGCAGTGGTATCAACGCAGAGTAC
>R4
AAGCATATAGTAGAGATCGTACTCTGCGTTGATACCACTGCTT
>F5
AAGCAGTGGTATCAACGCAGAGTAC
>R5
AACAGCAGTATAGACTGTGTACTCTGCGTTGATACCACTGCTT
```

##### 1.  Copy or create a symbolic link of reads_of_insert.fasta from Part I to a new directory.

##### 2.  Run parameters for `pbtranscript.py classify`

```bash
pbtranscript.py classify 
    reads_of_insert.fasta       # input file
    isoseq_draft.fasta          # intermediate file (won't use downstream)
    --min_seq_len 300           # minimum length threshold
    --cpus 15                   # number of threads to use (max-1)
    --flnc isoseq_flnc.fasta    # full-length output fasta filename
    --nfl isoseq_nfl.fasta      # non-full-length output fasta filename
    -p custom_primers.fasta     # Only if multiplexed samples, see step 0.
    --report PRIMERREPORTFN     # CSV file to output primer info (default: *.primer_info.csv
    --summary SUMMARY_FN        # TXT file to output classsify summary (default: *.classify_summary.txt
```

Advanced [(see wiki for details)](https://github.com/PacificBiosciences/cDNA_primer/wiki/RS_IsoSeq-%28v2.3%29-Tutorial-%231.-Getting-full-length-reads)

```
    [--ignore_polyA]            # For sequencing RT products w/no polyAs
    [--min_score MIN_SCORE]
    [--detect_chimera_nfl]      # Cleaner data, but super resource hungry!
```

**Example**

**Note:** assuming all input files in current dir or soft-links:
```bash
pbtranscript.py classify LT1_reads_of_insert.fasta LT1_isoseq_draft.fasta \
-p bc_primers.fasta --cpus 15 --min_seq_len 300 --flnc LT1_isoseq_flnc.fasta \
--nfl LT1_isoseq_nfl.fasta
```


### Part IIIb: Demultiplexing the samples
---

Only need to do if your samples were barcoded, pooled prior to library prep and sequenced together.  This step requires that you added the custom_barcoded file seen in the classify step

##### 1.  Download the 'DemultiplexFlncFasta.py' cloning into home directory
*   `git clone https://github.com/puapinyoying/IsoSeqScripts.git`

##### 2.  Give it run permissions
*   `chmod 755 DemultiplexFlncFasta.py`

##### 3.  For each sample input only the full length classified reads (e.g. LT1_isoseq_flnc.fasta)
*   `DemultiplexFlncFasta.py LT1_isoseq_flnc.fasta`

##### 4.  Demultiplexed samples will be in the `demultiplexed` subfolder in current directory


### Part IV: Cluster
---

This step will take your classified full-length reads and cluster them if their sequences are similar (same isoform).  Then the full-length reads are polished the non-full-length reads to increase read quality using the ICE & QUIVER algorythms

**Note**: For demultiplexed samples, the split sample files needs to be paired with the parent pooled non-full-length fasta file for the cluster/polish step. See explaination in 2. of this section.

##### 1.  Run parameters for cluster
```bash
pbtranscript.py cluster 
    isoseq_flnc.fasta                   # Input full-length non-chimeric reads in fasta format, used for clustering consensus isoforms
    final.consensus.fa                  # Output predicted (unpolished) consensus isoforms in fasta file. 
    --nfl_fa isoseq_nfl.fasta           # Input non-full-length reads in fasta format, used for polishing consensus isoforms
    -d cluster                          # Directory to store temporary and output cluster files.(default: output/)
    --ccs_fofn reads_of_insert.fofn     # A FOFN of ccs.h5 files, which contain quality values of reads of insert.
    --bas_fofn input.fofn               # A FOFN of bax/bas.h5 files (e.g., input.fofn), which contain quality values of raw reads and subreads
    --cDNA_size above3k                 # Estimated cDNA size. {under1k,between1k2k,between2k3k,above3k}  
    --quiver                            # Call quiver to polish consensus isoforms using non-full-length non-chimeric reads of insert.
    --blasr_nproc 15 
    --quiver_nproc 15
    --hq_isoforms_fa HQ_ISOFORMS_FA     # Quiver polished, high quality isoforms in fasta, default: root_dir/output/all_quivered_hq.fa
    --hq_isoforms_fq HQ_ISOFORMS_FQ     # Quiver polished, high quality isoforms in fastq, default: root_dir/output/all_quivered_hq.fq
    --lq_isoforms_fa LQ_ISOFORMS_FA     # Quiver polished, low quality isoforms in fasta,  default: root_dir/output/all_quivered_lq.fa
    --lq_isoforms_fq LQ_ISOFORMS_FQ     # Quiver polished, low quality isoforms in fastq, default: root_dir/output/all_quivered_lq.fq
    # --help
```

**Example command**:
```bash
pbtranscript.py cluster LT1_isoseq_flnc.fasta LT1_final.consensus.fasta --nfl_fa LT1_isoseq_nfl.fasta -d cluster_out --ccs_fofn LT1_reads_of_insert.fofn --bas_fofn LT1_input.fofn --cDNA_size above3k --quiver --blasr_nproc 15 --quiver_nproc 15
```

##### 2.  Optional: Perform if samples were demultiplexed

Demultiplexed example output files from previous steps:
```
*   LT1_isoseq_nfl.fasta            # Non-full-length fasta from classify
*   LT1_isoseq_flnc.fasta           # Full-length fasta from classify
    -   LT1_soleus_flnc.fasta       # barcode0_flnc.fasta (demultiplexed)
    -   LT1_edl_flnc.fasta          # barcode1_flnc.fasta (demultiplexed)
    -   LT1_cardiac_flnc.fasta      # barcode2_flnc.fasta (demultiplexed)
```

pbtranscript.py cluster must be done once per barcoded sample. Notice only the input flnc file is changing in each command. 

**Examples:**
```bash
pbtranscript.py cluster \
LT1_soleus_isoseq_flnc.fasta \
LT1_final.consensus.fasta \
--nfl_fa LT1_isoseq_nfl.fasta \
-d cluster_out \
--ccs_fofn LT1_reads_of_insert.fofn \
--bas_fofn LT1_input.fofn \
--cDNA_size above3k \
--quiver \
--blasr_nproc 15 \
--quiver_nproc 15

pbtranscript.py cluster \
LT1_edl_isoseq_flnc.fasta \
LT1_final.consensus.fasta \
--nfl_fa LT1_isoseq_nfl.fasta \
-d cluster_out \
--ccs_fofn LT1_reads_of_insert.fofn \
--bas_fofn LT1_input.fofn \
--cDNA_size above3k \
--quiver \
--blasr_nproc 15 \
--quiver_nproc 15

pbtranscript.py cluster \
LT1_cardiac_isoseq_flnc.fasta \
LT1_final.consensus.fasta \
--nfl_fa LT1_isoseq_nfl.fasta \
-d cluster_out \
--ccs_fofn LT1_reads_of_insert.fofn \
--bas_fofn LT1_input.fofn \
--cDNA_size above3k \
--quiver \
--blasr_nproc 15 \
--quiver_nproc 15
```

**Note**: The more full-length data going into cluster, the more isoforms you may find and more accurate they will be too. Therefore, it may be useful to run cluster on the multiplexed group too. Might be useful for screening out isoform artifacts, or for finding very rare transcripts. However, there is a balance.  Running cluster with too many full length reads (> 12 cells) is computationally intensive and will take a very long time or out right fail.


### Part V: Setup and alignment with STARlong
---
STAR is a fast and accurate splice aware aligner. STARlong, is an alternate algorithm for handling long reads and is perfect for our situation. See PacBio's official comparison of [GMAP vs STARlong](https://github.com/PacificBiosciences/cDNA_primer/wiki/Bioinfx-study%3A-Optimizing-STAR-aligner-for-Iso-Seq-data)

##### 1.  Download reference genome files and install STAR
*   See requirements, in Part 0 for details.

##### 2.  **Index the genome** 
* Genomes must be indexed by the aligner you plan to use before they can be used. Check the [STAR pdf manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for more detailed options

**Example STAR genome indexing command**:
```bash
STAR --runThreadN 15 --runMode genomeGenerate \
--genomeDir /gencode/mouse/index \
--genomeFastaFiles /gencode/mouse/GRCm38.primary_assembly.genome.fa \
--sjdbGTFfile /gencode/mouse/gencode.vM8.primary_assembly.annotation.gtf \
--sjdbOverhang 100
```

##### 3.  Alignment to genome

**STARlong example command**:
```bash
STARlong --runThreadN 15 \
--genomeDir /data/references/gencode/mouse/mm38 \
--runMode alignReads \
--outSAMattributes NH HI NM MD \
--readNameSeparator space \
--outFilterMultimapScoreRange 1 \
--outFilterMismatchNmax 2000 \
--scoreGapNoncan -20 \
--scoreGapGCAG -4 \
--scoreGapATAC -8 \
--scoreDelOpen -1 \
--scoreDelBase -1 \
--scoreInsOpen -1 \
--scoreInsBase -1 \
--alignEndsType Local \
--seedSearchStartLmax 50 \
--seedPerReadNmax 100000 \
--seedPerWindowNmax 1000 \
--alignTranscriptsPerReadNmax 100000 \
--alignTranscriptsPerWindowNmax 10000 \
--readFilesIn /data/brian/test/all_quivered_hq.100_30_0.99.fastq
```

Output should be an `Aligned.out.sam` file by default


### Part VI: File conversions with SamTools
---
To use matchAnnot, and the collapse isoform step requires the aligned SAM as input. However, the SAM must be sorted first. If you plan to view the data in a genome browser such as IGV, you will also need to convert the SAM to a BAM file which must be sorted and indexed.  I am not sure how to sort a SAM file without first converting it into a BAM file and then converting it back (basic unix sort command may not account for the SAM header...).  

##### 1.  Convert a SAM to BAM (`-b flag to stdout bam`; `-S <input.sam>`)
*   `samtools view -b -S file.sam > file.bam`

##### 2.  Sort BAM file (`-m <intG> gigs of ram`; `-@ <int> threads/cores`)
*   `samtools sort -m 60G -@ 15 file.bam outputPrefix` 

##### 3.  Index the BAM
*   `samtools index sorted.bam`

*   A one-liner (shortcut) for step 1-3
    -   `samtools view -Sb Aligned.out.sam | samtools sort -m 60G -@ 15 - sorted && samtools index sorted.bam`

##### 4.  Convert sorted and indexed BAM back to SAM
*   `samtools view -h sorted.bam > sorted.sam`

**Examples**
```bash
samtools view -Sb LT1_soleus_aligned.sam | samtools sort - LT1_soleus_sorted && samtools index LT1_soleus_sorted.bam
samtools view -h LT1_soleus_sorted.bam > LT1_soleus_sorted.sam
```

### Part VII: Collapse Redundant Isoforms
---

***NOTE***: IsoSeq analysis using exCOVator.py and exPhaser made use of uncollapsed data. We found that collapsing the isoforms results in loss of quantitative ability. Thus we did not end up using part VII, VIII, IV for the project in the manuscript. Please skip directly to [Part X: HTSeq-count](#part-x-htseq-count)

This step requires installation of pbtranscript-TOFU. [Official documentation](https://github.com/PacificBiosciences/cDNA_primer/wiki/tofu-Tutorial-%28optional%29.-Removing-redundant-transcripts)

`pbtranscript.py cluster`, is conservative when it groups consensus reads.  Only exact reads of the same exon structure AND length are merged.  However, some reads may have longer 5' ends. Therefore, the shorter reads can be considered subsets of the longer reads. Displaying reduntant information presents a lot of noise during analysis.  It can help to collapse these shorter redundant reads into one representative (longest) read of the set. 

The script ['collapse_isoforms_by_sam.py'] requires:

*   Sorted SAM file (post alignment with STARlong or GMAP)
*   Matching fasta or fastq from QUIVER step. ('all_quivered_hq.100_30_0.99.fastq')

Before entering any commands, you must load the smrtanalysis enviornmental variables, the smrtshell and the virtual enviornment containing TOFU that you setup.

1. `source /opt/smrtanalysis/current/etc/setup.sh`
2. `/opt/smrtanalysis/smrtcmds/bin/smrtshell`
3. `source /home/smrtcore/vtofu/VENV_TOFU/bin/activate`

Test with `collapse_isoforms_by_sam.py --help`

**Command format:**

*   `collapse_isoforms_by_sam.py --input all_quivered_hq.100_30_0.99.fastq --fq -s Aligned.out.sorted.sam -o Aligned.out`

"The script also has a minimum threshold for alignment accuracy (`default: 0.85`) and coverage (`default: 0.99`), which you can alter using the `--min-coverage` and `--min-identity` options."

**Examples**:
```bash
collapse_isoforms_by_sam.py --input LT1_pooled_hq_quivered.fastq --fq -s LT1_pooled_sorted.sam -o LT1_pooled_collapsed
collapse_isoforms_by_sam.py --input LT1_cardiac_hq_quivered.fastq --fq -s LT1_cardiac_sorted.sam -o LT1_cardiac_collapsed
collapse_isoforms_by_sam.py --input LT1_edl_hq_quivered.fastq --fq -s LT1_edl_sorted.sam -o LT1_edl_collapsed
collapse_isoforms_by_sam.py --input LT1_soleus_hq_quivered.fastq --fq -s LT1_soleus_sorted.sam -o LT1_soleus_collapsed
```

**Output**

*   "The output of the script contains a GFF, a group file, a representative FASTA file, and the list of sequences that were ignored because it was unmapped or below the threshold."


### Part VIII: STARlong and Samtools revisited
---
Do a another round of alignments using STARlong. The same parameters as Part IV, step 4 apply. The only replacement will be the input fastq file from the collapse_by_sam script. Follow the samtools steps in Part V on the new collapsed aligned.sam file, in order to convert the sam into a bam file.  Then sort and index the bam and convert it back into a SAM file for matchAnnot.


### PART IX: MatchAnnot
---
[MatchAnnot](https://github.com/TomSkelly/MatchAnnot) is a useful tool for novel isoform discovery. Each transcript/collapsed-consensus read that mapped to a gene is compared to all the known annotations for that particilar gene and scored based on sequence similarity (e.g. number of exons, position of exons, length of exons, etc...). A score of 0 is means the transcript is completely different from the annotation, and a score of 5 is an exact match.  What we want is a score that is somewhere in between.  Too novel and it may be an artifact, too similar and it may not be interesting. In addition to the scores, the number of supporting full-length reads for the consensus read can be used to lend credibility to the potential novel isoform. The more full-length reads, the less likely the consensus read is an artifact or error of sequencing.  Read more about MatchAnnot at the [wiki](https://github.com/TomSkelly/MatchAnnot/wiki)

**MatchAnnot expects the following inputs**:
```
--gtf               #   Annotation file, in format as described by --format option (Mandatory).
--format            #   Format of annotation file: 'standard', 'alt' or 'pickle' (default: standard).
--clusters          #   cluster_report.csv as produced by IsoSeq (Optional).
<stdin or file>     #   SAM file of IsoSeq transcripts aligned to genomic reference (Mandatory).
--outpickle         #   Output pickle file for clusterView
```
More details about [input](https://github.com/TomSkelly/MatchAnnot/wiki/How-to-Run-matchAnnot)  

The SAM file is output from STARlong 
*   Must be sorted by position first. See part VIII

The MatchAnnot output is very detailed, please read the author's [guide](https://github.com/TomSkelly/MatchAnnot/wiki/How-to-Interpret-matchAnnot-Output)

**Command example:**
```bash
cat /data/brian/pooled/collapsed_sorted.sam | matchAnnot.py \
--cluster /data/brian/cluster/cluster_report.csv \
--gtf /data/references/gencode/mouse/gencode.vM8.primary_assembly.annotation.gtf \
--outpickle matchAnnot_results.pickle > matchAnnot_results.txt
```

### PART X:  HTSeq-count
---
Diving into matchAnnot output can be daunting without some kind of system. I found that one way to efficiently prioritize genes of interest is to rank them by the number of unique isoforms aligned to them first, then chosing biologically relevant genes within the top ranked. I would run htseq-count for all samples, place output into Excel columns for each gene and sort.

**Important Parameters from "htseq-count --help':**
```
-f SAMTYPE, --format=SAMTYPE
                    type of <alignment_file> data, either 'sam' or 'bam'
                    (default: sam)

-s STRANDED, --stranded=STRANDED
                    whether the data is from a strand-specific assay.
                    Specify 'yes', 'no', or 'reverse' (default: yes).
                    'reverse' means 'yes' with reversed strand
                    interpretation

-a MINAQUAL, --minaqual=MINAQUAL
                    skip all reads with alignment quality lower than the
                    given minimum value (default: 10)

-t FEATURETYPE, --type=FEATURETYPE
                    feature type (3rd column in GFF file) to be used, all
                    features of other type are ignored (default, suitable
                    for Ensembl GTF files: exon)

-i IDATTR, --idattr=IDATTR
                    GFF attribute to be used as feature ID (default,
                    suitable for Ensembl GTF files: gene_id)

-m MODE, --mode=MODE  mode to handle reads overlapping more than one feature
                    (choices: union, intersection-strict, intersection-
                    nonempty; default: union)

-o SAMOUT, --samout=SAMOUT
                    write out all SAM alignment records into an output SAM
                    file called SAMOUT, annotating each line with its
                    feature assignment (as an optional field with tag
                    'XF')
```

**Example**: Count reads per feature (I chose to go by gene_name instead of gene_id)
*   `htseq-count -i gene_name sorted.sam mm10_sorted.gtf > gene_count.txt`

For more ways to define the bahavior and in depth explainations of htseq-count see their [official guide](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)


### PART XI: Visualizing Data
---

Two ways I have been using to explore the data has been through the integrated genome viewer (IGV) and clusterView. IGV for general browsing and clusterview for individual genes.

#### Integrated Genome Visualizer (IGV)

##### Setup to view very large genes
1.  Increase memory to at least 4gb (or however much you have available, you can do this by editing the launch script. Open the 'igv.sh' file in linux/unix (or the igb.bat file in windows)

2.  Find the line containing the java command:
    *   `java  -Xmx1200m ...`
    *   Adjust by changing the value in mb `-Xmx1200m` to `-Xmx4000m` for 4 gigs. You can use more memory if available, but be sure to leave at least 2gb for your OS or it may slow down dramatically!

3.  Adjust the [preferences](http://www.broadinstitute.org/igv/Preferences) in IGV once launched to optimize viewing
    *   Specifically under the 'alignments tab', change the 'visibility range' to accomindate the largest gene of interest and uncheck downsampling

##### Setup the annotation files (GTF)

1.  When searching for novel/unannotated transcripts, you should obtain annotation files from multiple databases/sources including the one used for your alignment (e.g. gencode.gtf). IGV automatically includes RefSeq, so focus on other databases.
    *   UCSC Gene Tables allows you to export many annotations in GTF format Their known genes is comprehensive and is useful.
    *   Gencode builds on top of Ensembl's annotations, so Ensembl can be optional

2.  Sort the annotation file by going to  "Tools > Run igvtools"
    Command "Sort"

3.  Then use the same dropdown to 'index' the gtf file

4.  Load up the sorted GTFs and right-click on the track, select expand and then choose set visibility window again, set this to the largest gene of interest (e.g. 4000 for DMD)

**Here are a few tips from the FAQ:**

**Q.** When I launch IGV with the 2GB or 10GB options, I get an error that says, "Could not create the Virtual Machine."   I have plenty of memory on my computer so why is this happening?  
**A.** You need both a 64-bit OS and a 64-bit version of Java.  On many computers, 32-bit Java is installed by default, even if the OS is 64-bit.

**Q.** A gene locus often produces several splicing isoforms, so how do I show different splicing isoforms for a single gene locus in IGV?
**A.** Right-click on the track and select Expand Track, or click the triangle icon () on the left-hand side of the track.  This displays overlapping features in separate rows so that separate isoforms are visible.

If there are additional problems, check the full [IGV FAQ](https://www.broadinstitute.org/igv/FAQ).

##### ClusterView

ClusterView will output an image of the input gene along with annotations lined up horizontally for easy comparisons.  It also removes the introns and replaces them with a vertical redline to separate exons.  This saves space and makes it easier to compare isoforms.  However, I find it may not work well with large genes or with genes that have many unique isoforms. 

***WARNING***: Bug in current version of clusterview (09/24/2015)
When running `clusterView.py` you may run into an error.

```
raise RuntimeError ('no length in name: %s' % cluster.name)
RuntimeError: no length in name: c930/f10p6/1846
```

Reason for this, regular expression error in `cluster_view.py`.[See issue here.](https://github.com/TomSkelly/MatchAnnot/issues/4) However, we can fix this ourselves.

Line **41** in `MatchAnnot/clusterView.py` has an error that can be manually fixed by opening up the python script and changing this:
*   `REGEX_LEN = re.compile ('\|(\d+)$') # cluster length in cluster name`
    
To this:
*   `REGEX_LEN = re.compile ('\/(\d+)$') # cluster length in cluster name`

**Parameters**:
```
Usage: clusterView.py [options]

Options:
  --version          show program's version number and exit
  -h, --help         show this help message and exit
  --gtf=GTF          annotations file, in format specified by --format (optional)
  --format=FORMAT    format of annotation file: standard, alt, pickle (def: standard)
  --matches=MATCHES  pickle file from matchAnnot.py (optional)
  --gene=GENE        gene to plot (required)
  --omit=OMIT        clusters to ignore, e.g., c1234,c2345,c3456
  --show=SHOW        clusters to force shown, even if underpopulated
  --howmany=HOWMANY  how many clusters to plot (def: all)
  --nodups           discard exact duplicate clusters (def: keep all)
  --minlen=MINLEN    minimum length of plotted cluster (def: all)
  --maxlen=MAXLEN    maximum length of plotted cluster (def: all)
  --output=OUTPUT    output plot file name (def: exons.png)
  --flip             reverse plot orientation (def: mRNA 5' on left)
  --yscale=YSCALE    amount by which to scale Y axis (def: 1.0)
  --details=DETAILS  output file name for details in text format (def: no output)
  --fasta=FASTA      output directory for fasta files for chosen clusters (def: no output)
  --title=TITLE      title for top of figure
  --notes=NOTES      file of notations to add to plot (experts only!)

    Note: omit, show, how many parameters can be useful for readability of figure,
        especially if there is a very large artifact or outlier.
```

**Authors's example:**
```
~skellytf/MatchAnnot-master/clusterView.py \
    --gtf /is2/projects/pacbio/scratch/skellytf/data/makeGTF/GENCODE-21/gencode.v21.annotation_with_As_cached_65.pickle \
    --format pickle \
    --matches annotations.pickle \
    --gene KRAS \
    --howmany 100 \
    --omit c12345,c23456 \
    --yscale 0.7 \
    --details cv_KRAS_details_a549_F02_poly_65.txt \
    --output  cv_KRAS_plot_a549_F02_poly_65.png \
    --title "Top 100 transcripts for KRAS/a549, run 2015-05-19a_836/F02_1 (Ensembl GRCh38) with polyAs at 65%"
```

**My Example:**
```
clusterView.py --gtf gencode.vM8.primary_assembly.annotation.gtf \
--format standard --matches M2_pooled_matchAnnot_results.pickle \
--gene Dmd \
--yscale 0.7 \
--details Dmd_pooled.txt \
--output  Dmd_LT1_pooled.png \
--title "Dmd transcripts in soleus, edl, and cardiac muscles (Mouse GRCh38) with polyAs"
```

**How to interpret output**

*   **Green bars** - isoform from annotation (gtf)
*   **Blue bars** - from matchAnnot output (sample data)
*   If IsoSeq Quality scores are:
    -   0 or below 5 AND MatchAnnot transcript score is a perfect 5 = **Magenta bars**; if not perfect 5 = (light) **DodgerBlue bars**.
    -   Above threshold AND transcript score is perfect 5 = **purple bars**; if not perfect 5 = **Blue bars**

