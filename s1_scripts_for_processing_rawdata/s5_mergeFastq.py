import re
import os
import csv
import sys
import glob
import errno
import operator
import argparse
from Bio import SeqIO
from tqdm import *

# This script merges the cluster reports and fastq files within each size
# fraction and an option to further merge across size fractions.  The way it is 
# done is to alter the cluster names/number of each subsequent file merged after
# the first by adding the final cluster number of the previous file plus 1 to 
# the all the cluster numbers of to the next file. I add the 1 because often the
# first cluster number of each new file starts with cluster 0. (c0)

# E.g.last cluster number of previous file + 1) + current cluster number
# Example: file 1 has clusters from c0 to c5000. So last cluster number + 1 = 5001
# Now the cluster number of each row in the next file is captured and 5001 is 
# added to it.

# File 1 | Merged file
# --------------------
# c0    --> c0
# ...
# c5000 --> c5000

# File2  |
# --------
# c0    --> c5001
# ...
# c6000 --> c11001

# File3  |
# --------
# c0    --> c11002
# ...
# c4000 --> c15002

STEP4 = 's4_cluster'
STEP5 = 's5_mergeFastq'
CELL_DESC = '5-8cell_grps'

# /lustre/groups/cbi/brianu/projects/isoseq_muscle/results/s4_cluster/polyA/5-8cell_grps/7-15kb/7-15kb_grp2_soleus_6cell/cluster_out
REPORT_PATH =  'cluster_out/output/cluster_report.FL_nonFL.csv'
FASTQ_PATH = 'cluster_out/all_quivered_hq.100_30_0.99.fastq'

# groupLists, separated by sample tissue and size fraction
SOL_3KB = ['3-6kb_grp1_soleus_8cell','3-6kb_grp2_soleus_7cell','3-6kb_grp3_soleus_8cell']
EDL_3KB = ['3-6kb_grp1_edl_8cell','3-6kb_grp2_edl_7cell','3-6kb_grp3_edl_8cell']
CAR_3KB = ['3-6kb_grp1_cardiac_8cell','3-6kb_grp2_cardiac_7cell','3-6kb_grp3_cardiac_8cel']

SOL_5KB = ['5-10kb_grp1_soleus_8cell','5-10kb_grp2_soleus_8cell','5-10kb_grp3_soleus_8cell']
EDL_5KB = ['5-10kb_grp1_edl_8cell','5-10kb_grp2_edl_8cell','5-10kb_grp3_edl_8cell']
CAR_5KB = ['5-10kb_grp1_cardiac_8cell','5-10kb_grp2_cardiac_8cell','5-10kb_grp3_cardiac_8cell']

SOL_7KB = ['7-15kb_grp1_soleus_6cell','7-15kb_grp2_soleus_6cell']
EDL_7KB = ['7-15kb_grp3_edl_4cell', '7-15kb_grp5_edl_4cell', '7-15kb_grp5_edl_4cell', '7-15kb_grp6_edl_4cell']
CAR_7KB = ['7-15kb_grp7_cardiac_6cell','7-15kb_grp8_cardiac_5cell']

THREE = [SOL_3KB, EDL_3KB, CAR_3KB]
FIVE =  [SOL_5KB, EDL_5KB, CAR_5KB]
SEVEN = [SOL_7KB, EDL_7KB, CAR_7KB]

# CONSTANTS
DESC_REG = r'c(\d+)(/f\d+p\d+/\d+)( isoform=)c(\d+)(;.+)'
CLUST_REG = r'c(\d+)'
GRP_REG = r'^([\w\-]+)_grp\d+_(\w+)_\d+cell'

def checkFileExist(filePath):
    """Make sure the file exists"""
    if not os.path.isfile(filePath):
        fileName = os.path.basename(filePath)
        print "{0} was not found. Are you sure it was generated?".format(fileName)
        print "Be sure to run '--mergeWithin' before or at the same time as '--mergeAcross'"
        print "Full path to file: {0}".format(filePath)
        sys.exit(1)

def lowerUserList(inputList):
    """Lowercase user inputLists in case there are misspellings. (e.g. 3-6KB)""" 
    # clean list
    loweredList = []
    for item in inputList:
        loweredItem = item.lower()
        loweredList.append(loweredItem)

    return loweredList


def createFolder(fullDir):
    """Create subfolders for storing file of filenames and scripts"""
    # catch any erorrs with permissions, ignore if dirs already exist
    try:
        # os.makedirs works like 'mkdir -p' by making all parent dirs too
        os.makedirs(fullDir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def checkPolyA(ignorePolyA):
    """Check if --ignorePolyA flag was used. If so, change polyA directory to 
    no polyA""" 
    if ignorePolyA:
        dirPolyA = 'noPolyA'
    else: # Default is to not add the flag
        dirPolyA = 'polyA'

    return dirPolyA

def genSymLinks(groupList, resultsRoot, dirPolyA):
    """This function generates symbolic links from the cluster step to the 
    symlinks folder in the s5_mergeFastq directory.  Files in the cluster step
    are not named with any unique descriptors, so this is a chance to label them
    without duplicating the data for the merge steps. Function also returns a 
    tuple list of each cluster_report and fastq pair for each group, the size
    and tissue type of the samples."""

    groupTupList = []
    tissue = ''
    size = ''
    
    for groupName in tqdm(groupList):
        reResults = re.search(GRP_REG, groupName)
        size = reResults.group(1)
        tissue = reResults.group(2)

        srcGroupRoot = os.path.join(resultsRoot, STEP4, dirPolyA, CELL_DESC, size, groupName)
        srcReport = os.path.join(srcGroupRoot, REPORT_PATH)
        srcFastq = os.path.join(srcGroupRoot, FASTQ_PATH)

        # sanity checks
        checkFileExist(srcReport)
        checkFileExist(srcFastq)

        dstGroupRoot = os.path.join(resultsRoot, STEP5, 'clusterLinks', dirPolyA, size)
        
        reportName = "{0}_cluster_report.csv".format(groupName)
        fastqName = "{0}_quivered_hq.fastq".format(groupName)

        dstReport = os.path.join(dstGroupRoot, reportName)
        dstFastq = os.path.join(dstGroupRoot, fastqName)

        groupTupList.append((dstReport, dstFastq))

        # make the destination folder and symbolic links
        createFolder(dstGroupRoot) 

        # Only make symlinks if not made yet 
        if not os.path.isfile(dstReport):
            os.symlink(srcReport, dstReport)

        if not os.path.isfile(dstFastq):
            os.symlink(srcFastq, dstFastq)
        
    return groupTupList, tissue, size


def mergeReport(sampleTupList, mergedReportPath):
    """This function merges all the cluster report files that are part of the
    sample/group tuple list given. It's important to merge the cluster report
    first because it is more inclusive than the fastq files. Meaning, it often 
    has more unique clusters listed in it than the fastq files. It also returns
    a dictionary of the old and new cluster numbers that were altered during the
    merging step for each pair."""

    reportOutFile = open(mergedReportPath, 'wb')
    csvWriter = csv.writer(reportOutFile)
    clustKeyDict = {} # will use to change cluster #s in fastq
    repClustCount = 0 # Tracks final cluster number from previous file (int only)

    for tupIndex, (report, fastq) in enumerate(sampleTupList):
        reportInFile = open(report, 'rU')
        csvReader = csv.reader(reportInFile)
        
        # keep track of current cluster number 
        currClustNum = repClustCount
        currDict = {} # keep dictionary of unique old and new cluster numbers

        for rowIndex, row in enumerate(tqdm(csvReader)):
            if rowIndex == 0:
                csvReader.next() # skip header
            else:
                clustResults = re.search(CLUST_REG, row[0]) # regexp cluster num
                currNum = int(clustResults.group(1)) # conv output from str to int

                if tupIndex == 0: # if first file in sampleTupList
                    csvWriter.writerow(row) # just write directly
                    if currNum not in currDict: # so doesn't overwrite over and over
                        currDict[currNum] = currNum  # key & value are same in first file

                    currClustNum = currNum  # update with CURRENT cluster number

                else: # subsequently, we will continue numbering from repClustCount
                    # by adding prev file's total cluster count (+1)
                    newNum = currNum + repClustCount
                    row[0] = 'c{0}'.format(str(newNum))
                    csvWriter.writerow(row)  # write the new row to file
                    if currNum not in currDict:
                        currDict[currNum] = newNum  # key & val are old vs new now
                    currClustNum = newNum  # update NEW cluster number

        # we want the fastq file to be the main key values
        clustKeyDict[fastq] = currDict

        # Add 1 to repClustCount, most files start with cluster 0. 
        repClustCount = currClustNum + 1
        print "{0} clusters in file {1}".format(repClustCount, tupIndex+1)
        reportInFile.closed
    reportOutFile.closed
    repClustCount = 0 # set it back to zero for next file

    return clustKeyDict


def mergeFastq(sampleTupList, clustKeyDict, mergedFastqPath):
    """Merges the fastq files from the tuple list.  This time cluster number 
    changes are based on the clustKeyDict generated during the merging of the
    cluster reports. That way the numbers agree"""
    fastqWriter = open(mergedFastqPath, 'wb')
    
    for tupIndex, (report, fastq) in enumerate(sampleTupList):
        fastqReader = open(fastq, 'rU') 
        # get cluster dictionary so we can exchange cluster numbers in fastq
        # headers following the cluster report file.
        currDict = clustKeyDict[fastq] 

        # iterate over each sequence record in fastq file.
        for seqIndex, seq_record in enumerate(tqdm(SeqIO.parse(fastqReader, "fastq"))):
             # Capture current cluster name
            descResult = re.search(DESC_REG, seq_record.description) # capture
            currNum = int(descResult.group(1)) # capture current cluster number
            newNum = currDict[currNum] # use dict to get new number

            # create new seq record and replace the original/old records
            newSeqDesc = "c{0}{1}{2}c{0}{3}".format(newNum, descResult.group(2), 
                                        descResult.group(3), descResult.group(5))
            newSeqName = "c{0}{1}".format(newNum, descResult.group(2))
            seq_record.description = newSeqDesc
            seq_record.name = newSeqName
            seq_record.id = newSeqName

            # write the altered sequence record to the merged fastq file
            SeqIO.write(seq_record, fastqWriter, 'fastq') # Write to fq

        fastqReader.closed # close current fastq file
    fastqWriter.closed # close merged fastq file

def printClustKeyDict(clustKeyDict, mergedKeyPath):
    """Print the entire cluster dictionary to a file"""
    with open(mergedKeyPath, 'wb') as mergedKeyOutFile:
        csvWriter = csv.writer(mergedKeyOutFile)
        header = ['fastq_name', 'original', 'new']
        csvWriter.writerow(header)
        # Transform the {dict:{dict}} structure into a (tuple, {dict}) so it 
        # can be printed in order of original fastq file names
        sortedClustKeyTupDict = sorted(clustKeyDict.items(), key=operator.itemgetter(0))
        for fastq, clustDict in sortedClustKeyTupDict:
            for original, new in clustDict.iteritems():
                row = [fastq, original, new]
                csvWriter.writerow(row)


def mergeReportFastqKey(tupList, mergedDir, prefix1, prefix2):
    """Use as wrapper for generating the merged report, fastq and cluster key"""
    createFolder(mergedDir)

    mergedReportName = "{0}_{1}_merged_cluster_report.csv".format(prefix1, prefix2)
    mergedFastqName = "{0}_{1}_merged_quivered_hq.fastq".format(prefix1, prefix2)
    mergedKeyName = "{0}_{1}_merged_cluster_key.csv".format(prefix1, prefix2)
    
    mergedReportPath = os.path.join(mergedDir, mergedReportName)
    mergedFastqPath = os.path.join(mergedDir, mergedFastqName)
    mergedKeyPath = os.path.join(mergedDir, mergedKeyName)

    clustKeyDict = mergeReport(tupList, mergedReportPath)
    mergeFastq(tupList, clustKeyDict, mergedFastqPath)
    printClustKeyDict(clustKeyDict, mergedKeyPath)

def mergeTissueSameSize(groupList, resultsRoot, dirPolyA):
    """For use inside mergeWithinSize()"""
    print "Creating system links for",
    groupTupList, tissue, size = genSymLinks(groupList, resultsRoot, dirPolyA)
    print "{0} samples, in {1} size fraction.".format(tissue,size)
    # Generate merged result path and file names
    mergedDir = os.path.join(resultsRoot, STEP5, 'merged', dirPolyA, size)
    print "Merging {0} samples within the {1} size fraction...".format(tissue,size)
    mergeReportFastqKey(groupTupList, mergedDir, size, tissue)
    print
    



# open output fastq & cluster report files
def mergeWithinSize(mergeWithin, resultsRoot, ignorePolyA):
    """Main function for merging files within a size fraction. For example, if
    there were 10-16 cells ran in the fraction, the data would have been divided
    up into groups of 5-8 cells. (e.g. 3-6kb_grp1_soleus_8cell, 
    3-6kb_grp1_soleus_7cell).  Files for these groups of the same size and 
    tissue type will be merged into a single fastq, and cluster report.  The
    cluster key generated will be for tracing the old cluster names if needed."""
    
    if mergeWithin: # mergeWithin is an arg list of size fractions
        dirPolyA = checkPolyA(ignorePolyA)      
        for fractionName in mergeWithin:
            if fractionName =='3-6kb':
                for groupList in THREE:
                    mergeTissueSameSize(groupList, resultsRoot, dirPolyA)

            elif fractionName =='5-10kb':
                for groupList in FIVE:
                    mergeTissueSameSize(groupList, resultsRoot, dirPolyA)

            elif fractionName =='7-15kb':
                for groupList in SEVEN:
                    mergeTissueSameSize(groupList, resultsRoot, dirPolyA)

            else:
                print "Unknown fraction size: '{0}'".format(fractionName)
                print "Please make sure you typed in the options correctly."
                print "'--mergeWithin' takes arguments: '3-6kb', '5-10kb, 7-15kb"
                print "If you want to specify more than 1 size fraction please use space"
                print "as a delimiter. (e.g. '--mergeWithin 3-6kb 5-10kb 7-15kb')"
                sys.exit(0)

def getMergedWithinTupLists(mergeAcross, resultsRoot, dirPolyA):
    """Use this to obtain lists of full paths to the merged files from the 
    merged within sizes function (first step) depending on which combination of
    sizes are specified in --mergeAcross argument. The paths to the merged cluster 
    report, and fastq files are organized into tuples. Then they are stored into
    lists (tupLists) according to tissues. Tuples of files from different sizes
    are added to the tupLists based on user input from mergeAcross. All 3
    tupLists for each tissue are combined into an allTissueList and returned. 
    
    example: user input '--mergeAcross 3-6kb 5-10kb' 
    -------
    allTissueList = [solTupList, edlTupList, carTupList]
    solTupList = [(sol3kb_report.csv, sol3kb.fastq), (sol5kb_report.csv, sol5kb.fastq)]
    ...
    """
    
    # Directories for withinSize merged files (for source files)
    mergedDir = os.path.join(resultsRoot, STEP5, 'merged', dirPolyA)

    solTupList = []
    edlTupList = []
    carTupList = []

    merged3 = os.path.join(mergedDir, '3-6kb')
    merged5 = os.path.join(mergedDir, '5-10kb')
    merged7 = os.path.join(mergedDir, '7-15kb')

    # size directory list to check will depend on user input into mergeAcross
    sizeDirList = []
    if len(mergeAcross) == 3:
        if ('3-6kb' in mergeAcross) and ('5-10kb' in mergeAcross) and ('7-15kb' in mergeAcross):
            sizeDirList = [merged3, merged5, merged7]
        else:
            print 'Error: check your spelling.'
            print "Your input '--mergeAcross {0} {1} {2}'".format(mergeAcross[0],
                mergeAcross[1],mergeAcross[2])

    elif len(mergeAcross) == 2:
        if ('3-6kb' in mergeAcross) and ('5-10kb' in mergeAcross):
            sizeDirList = [merged3, merged5]
        elif ('3-6kb' in mergeAcross) and ('7-15kb' in mergeAcross):
            sizeDirList = [merged3, merged7]
        elif ('5-10kb' in mergeAcross) and ('7-15kb' in mergeAcross):
            sizeDirList = [merged5, merged7]
        else:
            print 'Error: check your spelling.'
            print "Your input '--mergeAcross {0} {1}'".format(mergeAcross[0], mergeAcross[1])

    elif len(mergeAcross) == 1:
        if mergeAcross[0] == 'all':
            sizeDirList = [merged3, merged5, merged7]
        else:
            print 'Error: check your spelling.'
            print "Your input '--mergeAcross {0}'".format(mergeAcross[0])
            print "Unless it's 'all', no single size allowed for '--mergeAcross'"

    for sizeDir in sizeDirList:
    # expect one match each
        solReport = sorted((os.path.join(sizeDir, '*soleus*report.csv')))
        solFastq = sorted((os.path.join(sizeDir, '*soleus*.fastq')))

        edlReport = sorted((os.path.join(sizeDir, '*edl*report.csv')))
        edlFastq = sorted((os.path.join(sizeDir, '*edl*.fastq')))

        carReport = sorted((os.path.join(sizeDir, '*cardiac*report.csv')))
        carFastq = sorted((os.path.join(sizeDir, '*cardiac*.fastq')))

        globList = [solReport, solFastq, edlReport, edlFastq, carReport, carFastq]
        for file in globList:
            checkFileExist(file[0])  # glob returns a list

        # zip turns two lists into one list of tuples
        solTup = zip(solReport, solFastq)
        edlTup = zip(solReport, solFastq)
        carTup = zip(solReport, solFastq)

        solTupList.append(solTup)
        edlTupList.append(edlTup)
        carTupList.append(carTup)

    allTissueList = [solTupList, edlTupList, carTupList]

    return allTissueList
    
def mergeTissueAcrossSize(mergeAcross, resultsRoot, ignorePolyA):
    """Merges the tissues across sizes.  This function requires that the
    samples have already been merged within each size first."""
    mergeAcross = lowerUserList(mergeAcross)

    dirPolyA = checkPolyA(ignorePolyA)
    # Generate merged result path and file names
    mergedParentDir = os.path.join(resultsRoot, STEP5, 'merged', dirPolyA)

    # lists of tuplelists from withinSize merged files
    # combolists of tuple lists. TupleLists are organized by tissue in the
    # following order: [sol, edl, car], in each tissue tuple list data is 
    # organized as (report, fastq)
    
    mergedDir = ''
    # All 3 fraction sizes
    if ('3-6kb' in mergeAcross) and ('5-10kb' in mergeAcross) and ('7-15kb' in mergeAcross):
        name = '3-6kb_5-10kb_7-15kb'
        mergedDir = os.path.join(mergedParentDir, name)
    
    # 3-6kb and 5-10 only
    elif ('3-6kb' in mergeAcross) and ('5-10kb' in mergeAcross) and ('7-15kb' not in mergeAcross):
        name = '3-6kb_5-10kb'
        mergedDir = os.path.join(mergedParentDir, name)

    # 3-6kb and 7-15 only
    elif ('3-6kb' in mergeAcross) and ('5-10kb' not in mergeAcross) and ('7-15kb' in mergeAcross):
        name = '3-6kb_7-15kb'
        mergedDir = os.path.join(mergedParentDir, name)

    # 5-10kb and 7-15 only
    elif ('3-6kb' not in mergeAcross) and ('5-10kb' in mergeAcross) and ('7-15kb' in mergeAcross):
        name = '5-10kb_7-15kb'
        mergedDir = os.path.join(mergedParentDir, name)

    # Start merging the files
    tissueDict = {0:'soleus', 1:'edl', 2:'cardiac'} # for prefix1 of mergeReportFastqKey()
    allTissueList = getMergedWithinTupLists(mergeAcross, resultsRoot, dirPolyA)
    for tissueIndex, tupList in enumerate(allTissueList):
        print "Merging {0} files across sizes...".format(tissueDict[tissueIndex])
        mergeReportFastqKey(tupList, mergedDir, tissueDict[tissueIndex], name)


### RUN ###

# Setup argparse
parser = argparse.ArgumentParser()
parser.add_argument('--mergeWithin', help="", nargs='+')
parser.add_argument('--mergeAcross', help="", nargs='+')

# Only True if the flag is added, otherwise is false
parser.add_argument('--ignorePolyA', help="", action='store_true')

# Essential args
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('--resultsRoot', help="", required=True)

args = parser.parse_args()

if not (args.mergeWithin or args.mergeAcross):
    parser.error('No action requested, add --mergeWithin, --mergeAcross or both')

if args.mergeWithin:
    mergeWithinSize(args.mergeWithin, args.resultsRoot, args.ignorePolyA)

if args.mergeAcross:
    mergeTissueAcrossSize(args.mergeAcross, args.resultsRoot, args.ignorePolyA)

print 'Done!'