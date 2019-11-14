#! /usr/bin/env python

import re
import os
import sys
import glob
import errno
import textwrap
import argparse

# Constants
STEP3 = 's3_classify'
STEP4 = 's4_cluster'
CELL_DESC = '5-8cell_grps'
FASTA_NAME_REGEXP = r'^([\w\-]+_grp\d+_\w+_\d+cell)_isoseq_\w+\.fasta'

def createFolder(fullDir):
    """Create subfolders for storing file of filenames and scripts"""
    # catch any erorrs with permissions, ignore if dirs already exist
    try:
        # os.makedirs works like 'mkdir -p' by making all parent dirs too
        os.makedirs(fullDir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

# /lustre/groups/cbi/brianu/projects/isoseq_muscle/results/s3_classify/<polyA or no>/5-8cell_grps/3-6kb/<grpName>

def genInputGrpDirLists(ignorePolyA, resultsRoot, size):
    """We want to build the directories to the grouped classify results (which
    will be where we find the input files),"""
    if ignorePolyA:
        dirPolyA = 'noPolyA'
    else: # Default is to not add the flag (polyA included is more stringent/recommended)
        dirPolyA = 'polyA'

    # generate full path to parent folder (fraction size) for input/classify results
    pathToInputGrpDirs = os.path.join(resultsRoot, STEP3, dirPolyA, CELL_DESC, size)
    currInputListing = os.listdir(pathToInputGrpDirs) # list of all files/folders in dir

    grpInputDirNames = [] # Will extract these dir names for naming files
    grpInputDirPaths = [] # For full path to the directories
    
    for item in currInputListing:
        currInputDir = os.path.join(pathToInputGrpDirs, item)
        if os.path.isdir(currInputDir): # filter for listing for directories only
            grpInputDirNames.append(item)
            grpInputDirPaths.append(currInputDir)

    return grpInputDirNames, grpInputDirPaths

def genFastaDict(grpInputDirPaths, demultiplex):
    """Create a dictionary of fasta files for making the cluster script. 
    fastaDict = {nflFastaPath: [FlncFastaPathList]}"""

    fastaPathDict = {} # dict keys will be nfl fasta : value will be flncFastas

    for i, grpDir in enumerate(grpInputDirPaths):
        flncFastaList = [] # List of flnc fasta input files (will be dict values)

        # non-full-length fastas always in parent group folder
        # There is only one per folder so will use its full path as key
        nflExpr = os.path.join(grpDir, '*_isoseq_nfl.fasta')
        nflGlob = sorted(glob.glob(nflExpr)) # can't risk files being out of order
        nflKey = nflGlob[0] # glob throws stuff into lists even if there is only 1 match
        
        if demultiplex: # use nested demultiplexed folder in path
            dirFlnc = 'demultiplexed'
        else: # Default is to not add the flag; 7-15kb are not pooled. So use grp folder
            dirFlnc = ''

        flncDir = os.path.join(grpDir, dirFlnc) # append 'demultiplexed' if true
        flncExpr = os.path.join(flncDir, '*_isoseq_flnc.fasta')
        flncGlobList = sorted(glob.glob(flncExpr))
        flncFastaList = [] # make new list with appended dirFlnc
        for path in flncGlobList:
            flncFastaList.append(path)

        fastaPathDict[nflKey] = flncFastaList
        # Note: This dictionary will contain one nfl key and one flnc if the group
        # was not demultiplexed and 3 if it was demultiplexed (for each tissue)

    return fastaPathDict # dictionary of fullpath to input files for cluster script

def genNameStemDict(fastaPathDict):
    """The FLNC fasta file names will be the basis for the cluster output folder
    names because it has the size, group, tissue/pooled, and number of cells in
    the group. Will use for making output directories and file output names. 
    Note: Although a bit long and unwieldly the NFL fasta full path will remain 
    the key, which will be useful when creating the script and can use the same 
    key for all dictionaries"""
    # r'^([\w\-]+_grp\d+_\w+_\d+cell)_isoseq_\w+\.fasta'
    nameStemsDict = {}
    for nflPath, flncPathList in fastaPathDict.iteritems():
        flncNameStems = [] # house list of flnc Name stems for building new dict
        for flncPath in flncPathList:
            flncBaseName = os.path.basename(flncPath)
            reResults = re.search(FASTA_NAME_REGEXP, flncBaseName)
            flncNameStem = reResults.group(1)
            flncNameStems.append(flncNameStem)

        nameStemsDict[nflPath] = flncNameStems
      
    return nameStemsDict

def genFofnDict(fastaPathDict, fofnRoot, size):
    """Create a dictionary of bax and ccs fofn fullpaths. It's rather complex, maybe
    there is a simpler solution out there, but this is how it is organized:
    {nflPath: {'ccs': ccs.fofn, 'bax': bax.fofn}}"""

    fofnCcsDir = os.path.join(fofnRoot, STEP4, CELL_DESC, 'ccs', size)
    fofnBaxDir = os.path.join(fofnRoot, STEP4, CELL_DESC, 'bax', size)

    fofnDict = {}
    for nflPath in fastaPathDict:
        nflBaseName = os.path.basename(nflPath)
        reResults = re.search(FASTA_NAME_REGEXP, nflBaseName)
        nflNameStem = reResults.group(1)
        
        baxFofn = nflNameStem + '_bax.fofn'
        ccsFofn = nflNameStem + '_ccs.fofn'
        fullBaxPath = os.path.join(fofnBaxDir, baxFofn)
        fullCcsPath = os.path.join(fofnCcsDir, ccsFofn)

        fofnDict[nflPath] = {'bax': fullBaxPath, 'ccs': fullCcsPath}

    return fofnDict

def makeOutDirs(ignorePolyA, nameStemsDict, resultsRoot, logRoot, scriptsRoot, size):
    """Make output directories for cluster results, logs, and scritps"""
    if ignorePolyA:
        dirPolyA = 'noPolyA'
    else: # Default is to not add the flag (polyA included is more stringent/recommended)
        dirPolyA = 'polyA'

    # make scripts folder
    scriptsDir = os.path.join(scriptsRoot, STEP4, dirPolyA, CELL_DESC, size)
    createFolder(scriptsDir)

    # make log folder
    logDir = os.path.join(logRoot, STEP4, dirPolyA, CELL_DESC, size)
    createFolder(logDir)

    # for cluster output
    clusterOutDict = {} # new dict for cluster directories, also tied to nfl path key
    clusterRoot = os.path.join(resultsRoot, STEP4, dirPolyA, CELL_DESC, size)
    for nflPath, flncNameStemList in nameStemsDict.iteritems():
        clusterOutList = []
        for flncNameStem in flncNameStemList:
            currClusterDir = os.path.join(clusterRoot, flncNameStem)
            clusterOutList.append(currClusterDir)
            createFolder(currClusterDir)
        clusterOutDict[nflPath] = clusterOutList

    return clusterOutDict, scriptsDir, logDir

def genClusterScripts(fastaPathDict, fofnDict, nameStemsDict, clusterOutDict, scriptsDir, logDir):
    """Generate the cluster submit scripts for colonial one."""

    # Loop through fastaPath dictionary, grab {nflPath : [flncPathList]}
    for nflPath, flncPathList in fastaPathDict.iteritems():
        # loop through list of flncPaths
        for i, flncPath in enumerate(flncPathList):
            
            # Output directory for cluster results
            clusterOutDir = clusterOutDict[nflPath][i]

            # Generate the file names and fullpath for the output files
            finalConsFasta = "{0}_final_consensus.fasta".format(nameStemsDict[nflPath][i])
            pathToFinalConsFasta = os.path.join(clusterOutDir, finalConsFasta)

            # unpack fofn paths
            baxFofn = fofnDict[nflPath]['bax']
            ccsFofn = fofnDict[nflPath]['ccs']

            # Add the full path the log, and script
            logOutFile = "{0}_log.out".format(nameStemsDict[nflPath][i])
            logErrFile = "{0}_log.err".format(nameStemsDict[nflPath][i])
            logOutPath = os.path.join(logDir, logOutFile)
            logErrPath = os.path.join(logDir, logErrFile)

            scriptName = '{0}_cluster.sh'.format(nameStemsDict[nflPath][i]) 
            scriptFullPath = os.path.join(scriptsDir, scriptName)


            output='''#!/bin/sh
# Cluster script for {0}

# Specify partition
#SBATCH -p 128gb

# Specify number of nodes 
#SBATCH -N 1

# Specify number of cores (each node has 16 cores)
#SBATCH -n 16

# one hour timelimit:
#SBATCH --time 14-00:00:00

cd {1}

echo "STARTED pbtranscript.py cluster: $(date)" >> {7}

# Load SmrtAnalysis software
module load smrt_analysis/2.3.0

# Load SmrtAnalysis enviornmental variables and the smrtShell
source /c1/apps/smrt_analysis/2.3.0/current/etc/setup.sh
/c1/apps/smrt_analysis/2.3.0/smrtcmds/bin/smrtshell

pbtranscript.py cluster \
{2} \
{3} \
--nfl_fa {4} \
-d cluster_out \
--bas_fofn {5} \
--ccs_fofn {6} \
--cDNA_size above3k \
--quiver \
--blasr_nproc 15 \
--quiver_nproc 15 \
>> {7} 2>> {8}

echo "COMPLETED pbtranscript.py cluster : $(date)" >> {7}'''.format(nameStemsDict[nflPath][i],
            clusterOutDir, flncPath, pathToFinalConsFasta, nflPath, baxFofn, ccsFofn, logOutPath, logErrPath)

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
    """Use argparse to handle user input for program"""
    
    # Create a parser object
    parser = argparse.ArgumentParser(
        prog="s4_genCluster.py",
        
        formatter_class=argparse.RawDescriptionHelpFormatter,
        
        description="""Python script to create submit scripts for cluster step.""")

    # Parse the polyA or the noPolyA groups
    parser.add_argument('--ignorePolyA', help="", action='store_true')
    parser.add_argument('--demultiplex', help="", action='store_true')
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
    requiredNamed.add_argument('--dataRoot', help="", required=True)
    requiredNamed.add_argument('--logRoot', help="", required=True)
    requiredNamed.add_argument('--scriptsRoot', help="", required=True)
    requiredNamed.add_argument('--fofnRoot', help="", required=True)
    requiredNamed.add_argument('--size', help="", required=True)
    args = parser.parse_args()
    
    return args

args = parseInput()

# Generate lists of input directories and paths based on polyA and size fraction
grpInputDirNames, grpInputDirPaths = genInputGrpDirLists(args.ignorePolyA, 
    args.resultsRoot, args.size)

# Build a dictionary of all the nfl and flnc fasta files, special directory for 
# demultiplexed files. Dictionary uses nfl.fasta path as key and flnc.fastas as values
fastaPathDict = genFastaDict(grpInputDirPaths, args.demultiplex)

# RegExp flnc.fasta name stems (grp, #cells, size ect) put into dict using nfl.fasta as the key
nameStemsDict = genNameStemDict(fastaPathDict)

# Make dict of ccs and bax fofn file paths using nfl.fasta path as key again
fofnDict = genFofnDict(fastaPathDict, args.fofnRoot, args.size)

# generate output directories for cluster, scripts and log files. The clusterOutDict
# has the same nfl.fasta key, with the values being cluster dirs
clusterOutDict, scriptsDir, logDir = makeOutDirs(args.ignorePolyA, nameStemsDict, 
    args.resultsRoot, args.logRoot, args.scriptsRoot, args.size)

# Now, generate the cluster scripts with all the created dicts and vars
genClusterScripts(fastaPathDict, fofnDict, nameStemsDict, clusterOutDict, scriptsDir, logDir)

# Finally, create a 'run all' script for the whole size fraction.
writeRunScript(scriptsDir, args.size)