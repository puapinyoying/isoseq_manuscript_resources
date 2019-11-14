## IsoSeqScripts
Custom scripts related to PacBio IsoSeq analysis
- mixture of makefiles, bash scripts and python 2.7 scripts

dependencies:
- python=2.7
- anaconda
- pandas
- tqdm
- biopython
- pip
- pip:
  - argparse

Script  |  Input  | Output |  Requirements |  Description  |
------  |  -----  | ------ | ------------- | ------------  |
s3_demultiplexFlncFasta.py  |  'isoseq_flnc.fasta'  | Folder of barcode fasta files | BioPython, and tqdm python packages | Divides the pooled full-length reads from the classify step into individual fasta files based on barcode/primer.  Requires that you added the custom primer/barcode file during the classify step. If done correctly, sequence descriptions should include the primer field*.   |

##### *Example header from isoseq_flnc.fasta file:
>\>m151120_074440_42136_c100836382550000001823182312311555_s1_p0/98678/5437_73_CCS strand=-;fiveseen=1;polyAseen=1;threeseen=1;fiveend=25;polyAend=5389;threeend=5415;**primer=0**;chimera=0
