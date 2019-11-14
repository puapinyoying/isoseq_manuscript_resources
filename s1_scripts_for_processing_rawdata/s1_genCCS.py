#! /usr/bin/env python

import os
import sys
import glob
import errno
import argparse

# Script will generate submit scripts for generating circular consensus reads
# for each individual cell. Each cell contains 3 bax.h5 files

# Constants
STEP1 = 's1_ccs'
CELL_DESC = 'indv_cell'

def createFolder(fullDir):
    """Create subfolders for storing file of filenames and scripts"""
    # catch any erorrs with permissions, ignore if dirs already exist
    try:
        # os.makedirs works like 'mkdir -p' by making all parent dirs too
        os.makedirs(fullDir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def buildNamesPaths(dataRoot, resultsRoot, size, tissue, cellList, runFolder, shortFolderName):
    """Construct paths and stem names for files"""
    dataPathList = []
    resultsPathList = []
    nameList = [] # more concise name for file

    for cell in cellList:
        # input data paths
        currDataPath = os.path.join(dataRoot, runFolder, cell, 'Analysis_Results')
        dataPathList.append(currDataPath)

        # Output data filename and paths
        # FileName example: 'size_tissue_shortFolderName_runCell.sh'
        # e.g. 3-6kb_pooled_M2F3_A01_1.sh
        currName = '{0}_{1}_{2}_{3}'.format(size, tissue, shortFolderName, cell)
        currResultsPath = os.path.join(resultsRoot, STEP1, CELL_DESC, size, runFolder, cell)
        
        resultsPathList.append(currResultsPath)
        nameList.append(currName)

    return dataPathList, resultsPathList, nameList

def genOutputDirs(resultsPathList, fofnRoot, scriptsRoot, logRoot, size):
    """Generate output directories for results, scripts, logs, and fofn files. 
    Function also returns the scripts, fofn, and log directories."""
    scriptsDir = os.path.join(scriptsRoot, STEP1, CELL_DESC, size)
    fofnDir = os.path.join(fofnRoot, STEP1, CELL_DESC, size)
    logDir = os.path.join(logRoot, STEP1, CELL_DESC, size)
    createFolder(scriptsDir)
    createFolder(fofnDir)
    createFolder(logDir)

    for path in resultsPathList:
        createFolder(path)

    return scriptsDir, fofnDir, logDir


def createFofn(dataPathList, fofnDir, nameList):
    """Generate file of filenames for CCS script input"""
    for j, dataPath in enumerate(dataPathList):
        fofnFileName = nameList[j] + '.fofn'
        # relative path to fofn directory
        fofnPath = os.path.join(fofnDir, fofnFileName) 
        baxPathList = sorted(glob.glob(dataPath + '/*.bax.h5'))

        with open(fofnPath, 'wb') as currFofn:
            currFofn.write('\n'.join(baxPathList)) # place each on new line


def genCcsScripts(nameList, resultsPathList, fofnDir, logDir):
    """Generate the CCS submit scripts for colonial one"""
    for i, name in enumerate(nameList):
        # form full path to current fofn
        fofnFileName = name + '.fofn'
        fofnFilePath = os.path.join(fofnDir, fofnFileName)

        # Log files
        logOutName = name + '_ccs_log.out'
        logErrName = name + '_ccs_log.err'
        logOutPath = os.path.join(logDir, logOutName)
        logErrPath = os.path.join(logDir, logErrName)

        # Results path
        resultsDir = resultsPathList[i]

        # Scripts path
        scriptName = '{0}_ccs.sh'.format(name)
        scriptFullPath = os.path.join(scriptsDir, scriptName)
        

        output='''#!/bin/sh

# Submit script to generate CCS reads for: {0}

# Specify debug partition
#SBATCH -p defq

# Specify number of nodes 
#SBATCH -N 1

# Specify number of cores (each node has 16 cores)
#SBATCH -n 16

# one hour timelimit:
#SBATCH --time 1-00:00:00

# Load SmrtAnalysis software
module load smrt_analysis/2.3.0

# Load SmrtAnalysis enviornmental variables and the smrtShell
source /c1/apps/smrt_analysis/2.3.0/current/etc/setup.sh
/c1/apps/smrt_analysis/2.3.0/smrtcmds/bin/smrtshell

echo "START RUN: $(date)" >> {1}

ConsensusTools.sh CircularConsensus \
--minFullPasses 0 \
--minPredictedAccuracy 75 \
--parameters /c1/apps/smrt_analysis/2.3.0/current/analysis/etc/algorithm_parameters/2015-11/ \
--numThreads 15 \
--fofn {2} \
-o {3} >> {1} 2>> {4}

echo "COMPLETE RUN: $(date)" >> {1}'''.format(name, logOutPath, fofnFilePath,
        resultsDir, logErrPath)
        with open(scriptFullPath,'wb') as scriptFile:
            scriptFile.write(output)

def writeRunScript(scriptsDir, size):
    # Get list of roi scripts
    scriptPattern = "{0}_*.sh".format(size)
    scriptList = sorted(glob.glob(os.path.join(scriptsDir, scriptPattern)))

    # Create the run script
    runScriptName = 'submitAll_{0}.sh'.format(size)
    fullRunScriptPath = os.path.join(scriptsDir, runScriptName)
    header = "#! /usr/bin/env bash\n\n"
    with open(fullRunScriptPath, 'wb') as runScript:
        runScript.write(header)
        for script in scriptList:
            runCmd = "sbatch {0}\n".format(script)
            runScript.write(runCmd)

def parseInput():
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('--cellList', nargs='*', help="", required=True)
    requiredNamed.add_argument('--dataRoot', help="", required=True)
    requiredNamed.add_argument('--resultsRoot', help="", required=True)
    requiredNamed.add_argument('--scriptsRoot', help="", required=True)
    requiredNamed.add_argument('--fofnRoot', help="", required=True)
    requiredNamed.add_argument('--logRoot', help="", required=True)
    requiredNamed.add_argument('--runFolder', help="", required=True)
    requiredNamed.add_argument('--tissue', help="", required=True)
    requiredNamed.add_argument('--size', help="", required=True)
    requiredNamed.add_argument('--shortFolderName', help="", required=True)

    args = parser.parse_args()

    return args

###########
### Run ###
###########

# Parse input from makefile
args = parseInput()

# Build lists of names and paths to use
dataPathList, resultsPathList, nameList = buildNamesPaths(args.dataRoot, 
    args.resultsRoot, args.size, args.tissue, args.cellList, args.runFolder, 
    args.shortFolderName)

# Generate output directories and return some dir variables
scriptsDir, fofnDir, logDir = genOutputDirs(resultsPathList, args.fofnRoot, 
    args.scriptsRoot, args.logRoot, args.size)

# Create the fofns and place in makescript/fofn folder
createFofn(dataPathList, fofnDir, nameList)

# Generate submit scripts for ccs jobs
genCcsScripts(nameList, resultsPathList, fofnDir, logDir)

# Write a bash run all script to make it easier to run entire size fraction
writeRunScript(scriptsDir, args.size)