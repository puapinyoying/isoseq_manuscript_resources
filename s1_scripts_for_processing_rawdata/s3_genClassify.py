#! /usr/bin/env python

# This script will include a branch point in the classify parameter to ignore
# the polyA signature for full-length read criteria. A full-length read by
# default requires detection of the 5' primer, 3' primer and polyA signature. By
# ignoring the polyA signature, it may be possible to obtain more full-length
# reads that are from the middle of a long transcript (rather than the 3' end).
# However, this increase in sensitivity will likely result in loss of specificity.
# In otherwords, you might get an more artifacts.

import re
import os
import sys
import glob
import errno
import argparse

# Constants
STEP2 = 's2_roi'
STEP3 = 's3_classify'
CELL_DESC = '5-8cell_grps'
FILE_REGEXP = r'^([\w\-]+)\.fasta'

# /lustre/groups/cbi/brianu/projects/isoseq_muscle/results/s3_classify/<polyA or no>/5-8cell_grps/3-6kb/<grpName>

def createFolder(fullDir):
    """Create subfolders for storing file of filenames and scripts"""
    # catch any erorrs with permissions, ignore if dirs already exist
    try:
        # os.makedirs works like 'mkdir -p' by making all parent dirs too
        os.makedirs(fullDir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def getRoiFilesAndStems(resultsRoot, size):
    roiDir = os.path.join(resultsRoot, STEP2, CELL_DESC, size)
    roiExpFasta = os.path.join(roiDir, '*.fasta') # don't need fastq for classify
    roiFastaList = glob.glob(roiExpFasta)

    nameStemList = [] # will use this to make folders for each individual sample
    for f in roiFastaList:
        currFileName = os.path.basename(f)
        currStem = re.search(FILE_REGEXP, currFileName)
        nameStemList.append(currStem.group(1))

    # roiFastaList - full path to each roi fasta file, used in genClassifyScripts
    return roiFastaList, nameStemList # Stem list used in makeOutputDirs


def makeOutputDirs(ignorePolyA, nameStemList, resultsRoot, scriptsRoot, logRoot, size):
    """Generates script, log and group output folders and returns them.  The group
    directories are output as a list of dir strings, the other two are plain strings.
    Directories are decided based on the decision to ignorePolyA or not"""

    if ignorePolyA:
        dirPolyA = 'noPolyA'
    else: # Default is to not add the flag (polyA included is more stringent/recommended)
        dirPolyA = 'polyA'

    scriptsDir = os.path.join(scriptsRoot, STEP3, dirPolyA, CELL_DESC, size)
    logDir = os.path.join(logRoot, STEP3, dirPolyA, CELL_DESC, size)
    resultsDir = os.path.join(resultsRoot, STEP3, dirPolyA, CELL_DESC, size)

    createFolder(scriptsDir)
    createFolder(logDir)

    # Generate directories for each roi group
    grpOutDirList = [] # will return to use when generating scripts
    for name in nameStemList:
        currGrpDir = os.path.join(resultsDir, name)
        createFolder(currGrpDir)
        grpOutDirList.append(currGrpDir)

    return scriptsDir, logDir, grpOutDirList

def genClassifyScripts(ignorePolyA, barcodePath, roiFastaList, nameStemList, 
    scriptsDir, logDir, grpOutDirList):
    """Generate the classify submit scripts for colonial one, there is an option
    to ignore the polyA signature at this point. Also if samples are pooled, 
    the barcode option should be used."""

    # Check for polyA option
    if ignorePolyA:
        polyAFlag = '--ignore_polyA'
    else: # Default is to not add the flag
        polyAFlag = ''

    # Check for barcode option
    if barcodePath is not None:
        barcodeFlag = '-p {0}'.format(barcodePath)
    else:
        barcodeFlag = ''

    for i, roiFasta in enumerate(roiFastaList):
        # Generate the file names for the output files
        draftFasta = "{0}_isoseq_draft.fasta".format(nameStemList[i])
        flncFile = "{0}_isoseq_flnc.fasta".format(nameStemList[i])
        nflFile = "{0}_isoseq_nfl.fasta".format(nameStemList[i])
        logOutFile = "{0}_log.out".format(nameStemList[i])
        logErrFile = "{0}_log.err".format(nameStemList[i])

        # Add the full path to the files
        draftPath = os.path.join(grpOutDirList[i], draftFasta)
        flncPath = os.path.join(grpOutDirList[i], flncFile)
        nflPath = os.path.join(grpOutDirList[i], nflFile)
        logOutPath = os.path.join(logDir, logOutFile)
        logErrPath = os.path.join(logDir, logErrFile)


        output='''#!/bin/sh
# Classify script for {0} {1}

# Specify debug partition
#SBATCH -p defq

# Specify number of nodes 
#SBATCH -N 1

# Specify number of cores (each node has 16 cores)
#SBATCH -n 16

# one hour timelimit:
#SBATCH --time 2-00:00:00

cd {2}

echo "STARTED pbtranscript.py classify: $(date)" >> {8}

# Load SmrtAnalysis software
module load smrt_analysis/2.3.0

# Load SmrtAnalysis enviornmental variables and the smrtShell
source /c1/apps/smrt_analysis/2.3.0/current/etc/setup.sh
/c1/apps/smrt_analysis/2.3.0/smrtcmds/bin/smrtshell

pbtranscript.py classify {3} {4} {5} \
--cpus 15 \
--min_seq_len 300 \
--flnc {6} \
--nfl {7} {1} \
>> {8} 2>> {9}

echo "COMPLETED pbtranscript.py Classify : $(date)" >> {8}'''.format(
		nameStemList[i], polyAFlag, grpOutDirList[i], roiFasta, draftPath,
        barcodeFlag, flncPath, nflPath, logOutPath, logErrPath)
    
        scriptName = '{0}_classify.sh'.format(nameStemList[i]) 
        scriptFullPath = os.path.join(scriptsDir, scriptName)
        with open(scriptFullPath,'wb') as scriptFile:
            scriptFile.write(output)

def writeRunScript(scriptsDir, size):
    # Get list of roi scripts
    classifyScriptPattern = "{0}_*.sh".format(size)
    classifyScriptList = glob.glob(os.path.join(scriptsDir, classifyScriptPattern))

    # Create the run script
    runScriptName = 'submitAll_{0}.sh'.format(size)
    fullRunScriptPath = os.path.join(scriptsDir, runScriptName)
    header = "#! /usr/bin/env bash\n\n"
    with open(fullRunScriptPath, 'wb') as runScript:
        runScript.write(header)
        for script in classifyScriptList:
            runCmd = "sbatch {0}\n".format(script)
            runScript.write(runCmd)

def parseInput():
    parser = argparse.ArgumentParser()
    # Required only for multiplexed/pooled samples (3-6KB & 5-10kb only)
    parser.add_argument('--barcodePath', help="")
    # Only True if the flag is added, otherwise is false
    parser.add_argument('--ignorePolyA', help="", action='store_true')

    # Essential args
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('--resultsRoot', help="", required=True)
    requiredNamed.add_argument('--scriptsRoot', help="", required=True)
    requiredNamed.add_argument('--logRoot', help="", required=True)
    requiredNamed.add_argument('--size', help="", required=True)

    args = parser.parse_args()

    return args

# Parse input args from makefile
args = parseInput()

# get the roi filenames, directories and capture stem names
roiFastaList, nameStemList = getRoiFilesAndStems(args.resultsRoot, args.size)

# Create output directories and return some
scriptsDir, logDir, grpOutDirList = makeOutputDirs(args.ignorePolyA,
    nameStemList, args.resultsRoot, args.scriptsRoot, args.logRoot, args.size)

# Generate the classify scripts
genClassifyScripts(args.ignorePolyA, args.barcodePath, roiFastaList, nameStemList, 
    scriptsDir, logDir, grpOutDirList)

# Write up a run all script
writeRunScript(scriptsDir, args.size)
