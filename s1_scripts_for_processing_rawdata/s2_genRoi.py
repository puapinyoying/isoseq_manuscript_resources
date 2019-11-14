# %load test_roi.py
import os
import sys
import glob
import errno
import argparse

# This script is to merge several ccs data cells into analysis groups I'll call 
# reads of insert files. Since this script contains information on which cells
# are being combined, I also generate file of filenames for the ccs.h5 and 
# bax.h5 files that are needed for step 4: clustering.

###############################
### CUSTOM FOR THIS DATASET ###
###############################

STEP1 = "s1_ccs"
STEP2 = "s2_roi"
STEP4 = "s4_cluster"
CELL_DESC = '5-8cell_grps'
SIZES = ['3-6kb', '5-10kb', '7-15kb']

# Run folders
M2T1 = "2015-11-02_M2_Loading_Titration1"
M2T2 = "2015-11-06_M2_Loading_Titration2"
M2F5 = "2015-11-17_M2_5-10kb_FullRun"
M2F3 = "2015-11-20_M2_3-6kb_FullRun"
E7 = "2016-06-29_EDL_7-15kb_8cell"
CE7 = "2016-07-07_Cardiac_EDL_7-15kb_16cell"
SE7 = "2016-07-13_Soleus_EDL_7-15kb_16cell"

# Run cells, divide them into groups 8 or less cells
P3A_CELL  = ["A01_1","B01_1","C01_1", "D01_1"]
P3B_CELL  = ["A01_1","B01_1","C01_1", "D01_1"]
P3C1_CELL = ["A01_1", "C01_1", "D01_1", "E01_1", "F01_1", "G01_1", "H01_1"] # B01_1 was bad cell
P3C2_CELL = ["A02_1", "B02_1", "C02_1", "D02_1", "E02_1", "F02_1", "G02_1", "H02_1"]

# 5-10kb, A = M2T1, B = M2T2, C = M2F5
P5A_CELL  = ["E01_1", "F01_1", "G01_1", "H01_1"]
P5B_CELL  = ["E01_1", "F01_1", "G01_1", "H01_1"]
P5C1_CELL = ["A01_1", "B01_1", "C01_1", "D01_1", "E01_1", "F01_1", "G01_1", "H01_1"]
P5C2_CELL = ["A02_1", "B02_1", "C02_1", "D02_1", "E02_1", "F02_1", "G02_1", "H02_1"]

# 7-15kb, S7 = SE7, E7A = E7, E7B = CE7, E7C = SE7, C7 = CE7
# Keep the sample types separate. Some samples will have 5 or 6
S7A1_CELL = ["A01_1", "B01_1", "C01_1", "D01_1", "E01_1", "F01_1"]
S7A2_CELL = ["G01_1", "H01_1", "A02_1", "B02_1", "C02_1", "D02_1"]

E7A_CELL  = ["A01_1", "B01_1", "C01_1", "D01_1", "E01_1", "F01_1", "G01_1", "H01_1"]
E7B_CELL  = ["E02_1", "F02_1", "G02_1", "H02_1"]
E7C_CELL  = ["E02_1", "F02_1", "G02_1", "H02_1"]

C7A1_CELL = ["A01_1", "B01_1", "C01_1", "D01_1", "E01_1", "G01_1"]
C7A2_CELL = ["H01_1", "A02_1", "B02_1", "C02_1", "D02_1"]

def createFolder(fullDir):
    """Create subfolders for storing file of filenames and scripts"""
    # catch any errors with permissions, ignore if dirs already exist
    try:
        # os.makedirs works like 'mkdir -p' by making all parent dirs too
        os.makedirs(fullDir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def genCcsFullPaths(resultsRoot, size, runFolder, runCells):
    """Generate list of directories/fullpaths to CCS data folders"""
    ccsCellPathList = []
    for cell in runCells:
        currPath = os.path.join(resultsRoot, STEP1, 'indv_cell', size, runFolder, cell)
        ccsCellPathList.append(currPath)
    
    return ccsCellPathList

def genBaxFullPaths(dataRoot, runFolder, runCells):
    """Generate list of directories/fullpaths to bax.h5 data folders (this will
    be nessesary for the future cluster step)"""
    baxCellPathList = []
    for cell in runCells:
        currPath = os.path.join(dataRoot, runFolder, cell, "Analysis_Results")
        baxCellPathList.append(currPath)
    
    return baxCellPathList

def makeOutputDirs(dataRoot, resultsRoot, fofnRoot, scriptsRoot, size):
    # build each directory path
    ccsFofnDir = os.path.join(fofnRoot, STEP4, CELL_DESC, 'ccs', size)
    baxFofnDir = os.path.join(fofnRoot, STEP4, CELL_DESC, 'bax', size)
    scriptsDir = os.path.join(scriptsRoot, STEP2, CELL_DESC, size)
    resultsDir = os.path.join(resultsRoot, STEP2, CELL_DESC, size)

    # Create the directories
    createFolder(ccsFofnDir)
    createFolder(baxFofnDir)
    createFolder(resultsDir)
    createFolder(scriptsDir)

    return resultsDir, scriptsDir, baxFofnDir, ccsFofnDir

def genRoiCmds(ccsCellPathList, baxCellPathList, resultsDir, prefixName):
    """Build cat commands to combine fastq and fasta files and get list of h5 files"""
    # will contain * in this list. Will have to expand later
    roiCcsH5List = [] # for ccs.h5 fofn
    roiBaxH5List = [] # for bax.h5 fofn

    # Start of the cat commands
    catRoiFasta = "cat " # leave space after
    catRoiFastq = "cat "
    for i, ccsCell in enumerate(ccsCellPathList):
        # Build middle of cat command
        catRoiFasta += "{0}/*.ccs.fasta ".format(ccsCell) # space after each string
        catRoiFastq += "{0}/*.ccs.fastq ".format(ccsCell)

        ccsH5Path = ccsCell + "/*.ccs.h5" # make expr for globbing ccs.h5 files
        roiCcsH5List.append(ccsH5Path)

        baxH5Path = baxCellPathList[i] + "/*.bax.h5" # make expr for globbing bax.h5 files
        roiBaxH5List.append(baxH5Path)

    # End of the cat commands
    catRoiFasta += "> {0}/{1}.fasta".format(resultsDir, prefixName)
    catRoiFastq += "> {0}/{1}.fastq".format(resultsDir, prefixName)

    return catRoiFasta, catRoiFastq, roiCcsH5List, roiBaxH5List

def writeFofn(prefixName, roiCcsH5List, roiBaxH5List, ccsFofnDir, baxFofnDir):
    """Create the file of filenames for ccs.h5 & bax.h5 files for cluster step"""
    ccsFofnFileName = prefixName + '_ccs.fofn'
    ccsFofnFullPath = os.path.join(ccsFofnDir, ccsFofnFileName)
    with open(ccsFofnFullPath,'wb') as ccsFofnFile:
        for ccsH5 in roiCcsH5List:
            # roiCcsH5List contians *, needs to be expanded for full file names
            expNames = sorted(glob.glob(ccsH5)) # yeilds list of 3 per cell
            for h in expNames:
                ccsFofnFile.write(h + '\n')

    # Bax fofn files
    baxFofnFileName = prefixName + '_bax.fofn'
    baxFofnFullPath = os.path.join(baxFofnDir, baxFofnFileName)
    with open(baxFofnFullPath, 'wb') as baxFofnFile:
        for baxH5 in roiBaxH5List:
            expNames = sorted(glob.glob(baxH5)) # yeilds list of 3 per cell
            for h in expNames:
                baxFofnFile.write(h + '\n')


def writeRoiScript(prefixName, catRoiFasta, catRoiFastq, scriptsDir):
    """ """
    roiScriptName = prefixName + '_roi.sh' 
    header = "#! /usr/bin/env"
    roiScriptFullPath = os.path.join(scriptsDir, roiScriptName)
    with open(roiScriptFullPath, 'wb') as roiScript:
        roiScript.write(header + '\n\n')
        roiScript.write("# Generate fasta roi\n")
        roiScript.write(catRoiFasta + '\n\n')
        roiScript.write("# Generate fastq roi\n")
        roiScript.write(catRoiFastq + '\n\n')


def genFofnAndRoiScript(ccsCellPathList, baxCellPathList, prefixName, resultsDir, scriptsDir):
    """Wrapper script for generating the roi script and fofn's """
    catRoiFasta, catRoiFastq, roiCcsH5List, roiBaxH5List = genRoiCmds(
        ccsCellPathList, baxCellPathList, resultsDir, prefixName) # Get lists

    # Generate roi script and FOFN files
    writeRoiScript(prefixName, catRoiFasta, catRoiFastq, scriptsDir)
    writeFofn(prefixName, roiCcsH5List, roiBaxH5List, ccsFofnDir, baxFofnDir)

def writeRunScript(scriptsDir, size):
    # Get list of roi scripts
    scriptPattern = "{0}_*.sh".format(size)
    scriptList = sorted(glob.glob(os.path.join(scriptsDir, scriptPattern)))

    # Create the run script
    runScriptName = 'runAll_{0}.sh'.format(size)
    fullRunScriptPath = os.path.join(scriptsDir, runScriptName)
    header = "#! /usr/bin/env bash\n\n"
    with open(fullRunScriptPath, 'wb') as runScript:
        runScript.write(header)
        for script in scriptList:
            runCmd = "bash {0}\n".format(script)
            runScript.write(runCmd)

def parseInput():
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('--size', help="", required=True)
    requiredNamed.add_argument('--scriptsRoot', help="", required=True)
    requiredNamed.add_argument('--fofnRoot', help="", required=True)
    requiredNamed.add_argument('--resultsRoot', help="", required=True)
    requiredNamed.add_argument('--dataRoot', help="", required=True)
    args = parser.parse_args()

    return args

#################
### Run Start ###
#################

args = parseInput()

# Build all the directories paths and make the folders
resultsDir, scriptsDir, baxFofnDir, ccsFofnDir = makeOutputDirs(
    args.dataRoot, args.resultsRoot, args.fofnRoot, args.scriptsRoot, 
    args.size)

# we want our reads of insert files to contain 5-8 cells so it won't get bogged
# down when performing the cluster step.
if args.size == SIZES[0]:
    # 3-6kb full paths (starred samples will be combined later)
    # genCcsFullPaths(resultsRoot, size, runFolder, runCells)
    p3a_ccsList  = genCcsFullPaths(args.resultsRoot, SIZES[0], M2T1, P3A_CELL)    # 4 cells
    p3b_ccsList  = genCcsFullPaths(args.resultsRoot, SIZES[0], M2T2, P3B_CELL)    # 4 cells
    p3ab_ccsList = p3a_ccsList + p3b_ccsList                                      # 8 cells*
    p3c1_ccsList = genCcsFullPaths(args.resultsRoot, SIZES[0], M2F3, P3C1_CELL)   # 7 cells
    p3c2_ccsList = genCcsFullPaths(args.resultsRoot, SIZES[0], M2F3, P3C2_CELL)   # 8 cells

    # For the Bax.h5 files
    # genBaxFullPaths(dataRoot, runFolder, runCells)
    p3a_baxList  = genBaxFullPaths(args.dataRoot, M2T1, P3A_CELL)    # 4 cells
    p3b_baxList  = genBaxFullPaths(args.dataRoot, M2T2, P3B_CELL)    # 4 cells
    p3ab_baxList = p3a_baxList + p3b_baxList                         # 8 cells*
    p3c1_baxList = genBaxFullPaths(args.dataRoot, M2F3, P3C1_CELL)   # 7 cells
    p3c2_baxList = genBaxFullPaths(args.dataRoot, M2F3, P3C2_CELL)   # 8 cells
    
    # grp names include the full path to file, samples are multiplexed
    grp1PrefixName = '3-6kb_grp1_pooled_8cell'
    grp2PrefixName = '3-6kb_grp2_pooled_7cell'
    grp3PrefixName = '3-6kb_grp3_pooled_8cell'

    # genFofnAndRoiScript(ccsCellPathList, baxCellPathList, prefixName, scriptsDir)
    genFofnAndRoiScript(p3ab_ccsList, p3ab_baxList, grp1PrefixName, resultsDir, scriptsDir)
    genFofnAndRoiScript(p3c1_ccsList, p3c1_baxList, grp2PrefixName, resultsDir, scriptsDir)
    genFofnAndRoiScript(p3c2_ccsList, p3c2_baxList, grp3PrefixName, resultsDir, scriptsDir)

    writeRunScript(scriptsDir, SIZES[0])
    
elif args.size == SIZES[1]:
    # 5-10kb full paths
    p5a_ccsList  = genCcsFullPaths(args.resultsRoot, SIZES[1], M2T1, P5A_CELL)    # 4 cells
    p5b_ccsList  = genCcsFullPaths(args.resultsRoot, SIZES[1], M2T2, P5B_CELL)    # 4 cells
    p5ab_ccsList = p5a_ccsList + p5b_ccsList                                      # 8 cells*
    p5c1_ccsList = genCcsFullPaths(args.resultsRoot, SIZES[1], M2F5, P5C1_CELL)   # 8 cells
    p5c2_ccsList = genCcsFullPaths(args.resultsRoot, SIZES[1], M2F5, P5C2_CELL)   # 8 cells

    p5a_baxList  = genBaxFullPaths(args.dataRoot, M2T1, P5A_CELL)     # 4 cells
    p5b_baxList  = genBaxFullPaths(args.dataRoot, M2T2, P5B_CELL)    # 4 cells
    p5ab_baxList = p5a_baxList + p5b_baxList                          # 8 cells*
    p5c1_baxList = genBaxFullPaths(args.dataRoot, M2F5, P5C1_CELL)   # 8 cells
    p5c2_baxList = genBaxFullPaths(args.dataRoot, M2F5, P5C2_CELL)    # 8 cells

    # grp names include the full path to file, samples are multiplexed
    grp1PrefixName = '5-10kb_grp1_pooled_8cell'
    grp2PrefixName = '5-10kb_grp2_pooled_8cell'
    grp3PrefixName = '5-10kb_grp3_pooled_8cell'

    # genFofnAndRoiScript(ccsCellPathList, baxCellPathList, prefixName, scriptsDir)
    genFofnAndRoiScript(p5ab_ccsList, p5ab_baxList, grp1PrefixName, resultsDir, scriptsDir)
    genFofnAndRoiScript(p5c1_ccsList, p5c1_baxList, grp2PrefixName, resultsDir, scriptsDir)
    genFofnAndRoiScript(p5c2_ccsList, p5c2_baxList, grp3PrefixName, resultsDir, scriptsDir)

    writeRunScript(scriptsDir, SIZES[1])

elif args.size == SIZES[2]:
    # 7-15kb full paths
    s7a1_ccsList = genCcsFullPaths(args.resultsRoot, SIZES[2], SE7, S7A1_CELL)   # 6 cells
    s7a2_ccsList = genCcsFullPaths(args.resultsRoot, SIZES[2], SE7, S7A2_CELL)   # 6 cells
    e7a_ccsList  = genCcsFullPaths(args.resultsRoot, SIZES[2], E7, E7A_CELL)     # 8 cells   
    e7b_ccsList  = genCcsFullPaths(args.resultsRoot, SIZES[2], CE7, E7B_CELL)    # 4 cells
    e7c_ccsList  = genCcsFullPaths(args.resultsRoot, SIZES[2], SE7, E7C_CELL)    # 4 cells
    e7bc_ccsList = e7b_ccsList + e7c_ccsList                                     # 8 cells*
    c7a1_ccsList = genCcsFullPaths(args.resultsRoot, SIZES[2], CE7, C7A1_CELL)   # 6 cells
    c7a2_ccsList = genCcsFullPaths(args.resultsRoot, SIZES[2], CE7, C7A2_CELL)   # 5 cells

    # Bax.h5 files
    s7a1_baxList = genBaxFullPaths(args.dataRoot, SE7, S7A1_CELL)     # 6 cells
    s7a2_baxList = genBaxFullPaths(args.dataRoot, SE7, S7A2_CELL)     # 6 cells
    e7a_baxList  = genBaxFullPaths(args.dataRoot, E7, E7A_CELL)       # 8 cells   
    e7b_baxList  = genBaxFullPaths(args.dataRoot, CE7, E7B_CELL)      # 4 cells
    e7c_baxList  = genBaxFullPaths(args.dataRoot, SE7, E7C_CELL)      # 4 cells
    e7bc_baxList = e7b_baxList + e7c_baxList                          # 8 cells*
    c7a1_baxList = genBaxFullPaths(args.dataRoot, CE7, C7A1_CELL)     # 6 cells
    c7a2_baxList = genBaxFullPaths(args.dataRoot, CE7, C7A2_CELL)     # 5 cells

    # grp names include the full path to file (sampels not multiplexed)
    grp1PrefixName = '7-15kb_grp1_soleus_6cell'
    grp2PrefixName = '7-15kb_grp2_soleus_6cell'
    grp3PrefixName = '7-15kb_grp3_edl_8cell'
    grp4PrefixName = '7-15kb_grp4_edl_8cell'
    grp5PrefixName = '7-15kb_grp5_cardiac_6cell'
    grp6PrefixName = '7-15kb_grp6_cardiac_5cell'

    # genFofnAndRoiScript(ccsCellPathList, baxCellPathList, prefixName, scriptsDir)
    genFofnAndRoiScript(s7a1_ccsList, s7a1_baxList, grp1PrefixName, resultsDir, scriptsDir)
    genFofnAndRoiScript(s7a2_ccsList, s7a2_baxList, grp2PrefixName, resultsDir, scriptsDir)
    genFofnAndRoiScript(e7a_ccsList,  e7a_baxList,  grp3PrefixName, resultsDir, scriptsDir)
    genFofnAndRoiScript(e7bc_ccsList, e7bc_baxList, grp4PrefixName, resultsDir, scriptsDir)
    genFofnAndRoiScript(c7a1_ccsList, c7a1_baxList, grp5PrefixName, resultsDir, scriptsDir)
    genFofnAndRoiScript(c7a2_ccsList, c7a2_baxList, grp6PrefixName, resultsDir, scriptsDir)

    writeRunScript(scriptsDir, SIZES[2])

else:
    print "ERROR: Unknown fraction size ({0})".format(args.size)
