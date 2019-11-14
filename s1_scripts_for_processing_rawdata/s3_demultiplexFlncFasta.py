#! /usr/env python

# Script to demultiplex barcoded full-length isoseq reads into their own 
# individual fasta files (barcodes are referred to as primers in the sequence
# files because they are technically custom primers. Will use primer/barcode
# interchangably)

import re
import os
import sys
import glob
import errno
import argparse
import textwrap
from tqdm import *
from Bio import SeqIO  # 

FILE_NAME_REGEXP = r'(.+)\..+'
PRIMER_REGEXP = r"primer\=(\d+)"
STEP3 = 's3_classify'
CELL_DESC = '5-8cell_grps'
GRP_NAME_REGEXP = r'^[\w\-]+_(grp\d+)_pooled_(\d+cell)'

def createFolder(fullDir):
    """Create subfolders for storing file of filenames and scripts"""
    # catch any erorrs with permissions, ignore if dirs already exist
    try:
        # os.makedirs works like 'mkdir -p' by making all parent dirs too
        os.makedirs(fullDir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def genGrpDirList(resultsRoot, ignorePolyA, size):
    """Generates a list of group directories for size fraction and if polyA was
    included or ignored as full-length criteria. Returns the names of the group
    folders and the fullpaths"""
    if ignorePolyA:
        dirPolyA = 'noPolyA'
    else: # Default is to not add the flag
        dirPolyA = 'polyA'

    # generate full path to parent folder (fraction size)
    pathToGrpDirs = os.path.join(resultsRoot, STEP3, dirPolyA, CELL_DESC, size)
    currListing = os.listdir(pathToGrpDirs) # list of all files/folders in dir

    grpDirNames = [] # will use these dir names for naming files
    grpDirPaths = []  # filter for listing for directories only
    for item in currListing:
        currDir = os.path.join(pathToGrpDirs, item)
        if os.path.isdir(currDir):
            grpDirNames.append(item)
            grpDirPaths.append(currDir)

    return grpDirNames, grpDirPaths

def getGrpCellNum(grpDirNames):
    grpNumList = []
    cellNumList = []
    for name in grpDirNames:
        reResults = re.search(GRP_NAME_REGEXP, name)
        grpNum = reResults.group(1)  # grp1, grp2, .. etc
        grpNumList.append(grpNum)

        cellNum = reResults.group(2) # 8cell, 5cell, .. etc
        cellNumList.append(cellNum)

    return grpNumList, cellNumList


def getDirsFastaLists(grpDirPaths):
    flncFastaList = [] # List of flnc fasta input files
    dmpDirList = [] #  For demultiplexed dir paths to house results
    
    for grpDir in grpDirPaths:
        expr = os.path.join(grpDir, '*_isoseq_flnc.fasta')
        globList = glob.glob(expr)
        flncFastaList.append(globList[0])

        currDir = os.path.join(grpDir, 'demultiplexed')
        createFolder(currDir)
        dmpDirList.append(currDir)

    return flncFastaList, dmpDirList


def demultiplexFasta(flncFastaList, dmpDirList, grpNumList, cellNumList, size):
    for i, fastaFile in enumerate(flncFastaList):
        # make some new file names for demultiplexed samples
        bar0_name = "{0}_{1}_soleus_{2}_isoseq_flnc.fasta".format(size, grpNumList[i], cellNumList[i])
        bar1_name = "{0}_{1}_edl_{2}_isoseq_flnc.fasta".format(size, grpNumList[i], cellNumList[i])
        bar2_name = "{0}_{1}_cardiac_{2}_isoseq_flnc.fasta".format(size, grpNumList[i], cellNumList[i])

        # Build full path to write output files
        bar0_path = os.path.join(dmpDirList[i], bar0_name)
        bar1_path = os.path.join(dmpDirList[i], bar1_name)
        bar2_path = os.path.join(dmpDirList[i], bar2_name)

        # Open all 3 barcode files to write to
        bar0_writer = open(bar0_path, "wb")
        bar1_writer = open(bar1_path, "wb")
        bar2_writer = open(bar2_path, "wb")

        # Loop thorugh each sequence header in the input file, if primer # in header
        # matches one of the barcodes, output the whole record to its respective file
        print "Demultiplexing: {0} ...".format(os.path.basename(fastaFile))
        for seq_record in tqdm(SeqIO.parse(fastaFile, "fasta")):
            # seq_record.description gets complete fasta header
            searchObj = re.search(PRIMER_REGEXP, seq_record.description)
            primerNum = searchObj.group(1) # grabs just primer number

            # Output each barcode's sequence records to their own fasta file
            if primerNum == '0':
                SeqIO.write(seq_record, bar0_writer, "fasta")

            elif primerNum == '1':
                SeqIO.write(seq_record, bar1_writer, "fasta")

            elif primerNum == '2':
                SeqIO.write(seq_record, bar2_writer, "fasta")

            else:
                print "That was unexpected! This script can only parse 6 barcodes"
                print seq_record.description
                print "Are you sure barcode{0} exists?".format(primerNum)

        # Close all files
        bar0_writer.closed
        bar1_writer.closed
        bar2_writer.closed
        print "Complete."

def parseInput():
    """Use argparse to handle user input for program"""
    
    # Create a parser object
    parser = argparse.ArgumentParser(
        prog="DemultiplexFlncFasta.py",
        
        formatter_class=argparse.RawDescriptionHelpFormatter,
        
        description="""Custom demultiplex script for the 'isoseq_flnc.fasta' 
        Fastas are output to a new folder called 'demultiplexed' in the 
        'isoseq_flnc.fasta' directory.""")

    # Parse the polyA or the noPolyA groups
    parser.add_argument('--ignorePolyA', help="", action='store_true')
    parser.add_argument("-v", "--version", action="version",
                        version=textwrap.dedent("""\
        %(prog)s
        -----------------------   
        Version:    0.3 
        Updated:    08/29/2016
        By:         Prech Uapinyoying   
        Website:    https://github.com/puapinyoying"""))

    # Essential args
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('--resultsRoot', help="", required=True)
    requiredNamed.add_argument('--size', help="", required=True)
    args = parser.parse_args()
    
    return args

#############
#### Run ####
#############

args = parseInput()

# Generate lists of directories for each group of samples to demultiplex
grpDirNames, grpDirPaths = genGrpDirList(args.resultsRoot, args.ignorePolyA, args.size)

# Get the group's number (grp1, grp2, etc) and the cell number (8cells, 5cells etc)
grpNumList, cellNumList = getGrpCellNum(grpDirNames)

# Get the list of fasta files and demultiplex directories
flncFastaList, dmpDirList = getDirsFastaLists(grpDirPaths)

# Carry out the demultiplexing on all the groups in the sample.
demultiplexFasta(flncFastaList, dmpDirList, grpNumList, cellNumList, args.size)
