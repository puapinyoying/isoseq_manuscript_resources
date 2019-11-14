#! /bin/python3

import os
import re
import sys
import HTSeq
import argparse
import textwrap
import itertools
import pandas as pd
from tqdm import tqdm
from pandas import Series, DataFrame
from collections import defaultdict

# Regular expression patterns for parsing SAM alignment read names (mainly use NOT_COLLAPSED)
COLLAPSED_REGEXP = r'[PB\.\d]+\|[\w\:\(\)\-\+]+\|c\d+\/f(\d+)p(\d+)\/\d+'
NOT_COLLAPSED_REGEXP = r'c\d+\/f(\d+)p(\d+)\/\d+'

##############################
### Data Loading functions ###
##############################

def check_file_exist(file_path):
    """Make sure the file exists"""
    if not os.path.isfile(file_path):
        file_name = os.path.basename(file_path)
        print(f"{file_name} was not found. Check path")
        print(f"Full path to file: {file_path}")
        sys.exit(1)

def load_gff(gff_path):
    check_file_exist(gff_path)
    gff_reader = HTSeq.GFF_Reader(gff_path, end_included=True)
    
    return gff_reader

def load_bed(bed_path):
    check_file_exist(bed_path)
    bed_reader = HTSeq.BED_Reader(bed_path) # also call it a sam reader for ease

    return bed_reader

def load_bam(bam_path):
    check_file_exist(bam_path)
    bam_reader = HTSeq.BAM_Reader(bam_path)

    return bam_reader

def get_filename(bamPath):
    '''Gets the filename from each bamfile without the path or the extension. Can use as default sample names.'''
    path, filenameExt = os.path.split(bamPath)
    filename, ext = os.path.splitext(filenameExt)
    
    return filename

def load_bam_list(bamPathList):
    '''Loads all the bams and returns a list of bamReader objects. Also returns a list of filenames'''
    bamReaderList = []
    bamFilenameList = []
    
    for bamPath in bamPathList:
        # Load bamReader objects from baths
        bamReader = load_bam(bamPath)
        bamReaderList.append(bamReader)
        
        # Grab filenames
        filename = get_filename(bamPath)
        bamFilenameList.append(filename)
        
    return bamReaderList, bamFilenameList


#########################################
### Read sorting based on input exons ###
#########################################


def get_bed_list(bedReader):
    '''Obtain a list of genomic intervals from the bed file and return number of bed objects '''
    bedList = []
    for bed in bedReader:
        bedList.append(bed)
    
    bedNum = len(bedList)
    
    return bedList, bedNum

def get_interval_range(bedList):
    '''Calculate the minimum start coordinate and maximum end coordinate of all the exons 
    from the input bed. Variable is called intervalRange which is a genomic interval'''
    smallestStart = 0
    largestEnd = 0
    chrom = ''
    strand = ''

    for i, bed in enumerate(bedList):
        if i == 0:
            smallestStart = bed.iv.start
            largestEnd = bed.iv.end
            chrom = bed.iv.chrom
            strand = bed.iv.strand
        else:
            if bed.iv.start < smallestStart:
                smallestStart = bed.iv.start
            if bed.iv.end > largestEnd:
                largestEnd = bed.iv.end
    
    intervalRange = HTSeq.GenomicInterval(chrom, smallestStart, largestEnd, strand)

    return intervalRange

def read_name_pattern(collapsed):
    '''Choose regular expression patters for format of pacbio read names'''
    if collapsed:
        name_pattern = COLLAPSED_REGEXP
    else:
        name_pattern = NOT_COLLAPSED_REGEXP
        
    return name_pattern

def grab_fullLength_count(name_pattern, readAlnObj):
    '''Gets the full-length read count from the pacbio read name'''
    result = re.search(name_pattern, readAlnObj.read.name)
    fullLengthCount = int(result.group(1))
    
    return fullLengthCount

def sort_reads(bamReader, intervalRange):
    '''Only select reads that contain the full range of all exons provided in bed file. 
    This sorting is rather strict because the read must contain the interval completely.
    If it overlaps the interval in a partial way it will be put in the
    readsOverlapInterval pile (total reads). Also keeps counts for stats.''' 
    # Contains does not mean all exons match, they just contain the coordinates of the intervalRange
    readsContainInterval = [] # note name changes to lists
    readsOverlapInterval = [] # This is less string
    
    # Keep Track for stats
    readStatsDict = {}
    readsContainInterval_clusterCount = 0
    readsContainInterval_flCount = 0
    readsOverlapInterval_clusterCount = 0
    readsOverlapInterval_flCount = 0
    allReads_clusterCount = 0
    allReads_flCounts = 0
    
    for read in tqdm(bamReader):
        flCount = grab_fullLength_count(NOT_COLLAPSED_REGEXP, read)
        if read.iv.contains(intervalRange):
            readsContainInterval.append(read)
            readsContainInterval_clusterCount += 1
            readsContainInterval_flCount += flCount
            
        if read.iv.overlaps(intervalRange):
            readsOverlapInterval.append(read)
            readsOverlapInterval_clusterCount += 1
            readsOverlapInterval_flCount += flCount
        
        allReads_clusterCount += 1
        allReads_flCounts += flCount
    
    # Fill readStatsDict
    readStatsDict['readsContainInterval_clusterCount'] = readsContainInterval_clusterCount
    readStatsDict['readsContainInterval_flCount']      = readsContainInterval_flCount
    readStatsDict['readsOverlapInterval_clusterCount'] = readsOverlapInterval_clusterCount
    readStatsDict['readsOverlapInterval_flCount']      = readsOverlapInterval_flCount
    readStatsDict['allReads_clusterCount']             = allReads_clusterCount
    readStatsDict['allReads_flCounts']                 = allReads_flCounts
            
    return readsContainInterval, readsOverlapInterval, readStatsDict


def sort_reads_from_bamReaderList(bamReaderList, bamFilenameList, intervalRange, addReadStart=0, addReadEnd=0):
    '''Sort the reads of the whole bamReaderList and return a dictionary of lists 
    for each sample that contain the intervalRange and reads that do not. This interval 
    can be adjusted using the addReadStart, addReadEnd. '''
    # Adjust intervalRange for sorting reads
    readIntervalRange = intervalRange.copy() # make a copy called readIntervalRange so not to clobber original range
    readIntervalRange.start += addReadStart
    readIntervalRange.end += addReadEnd
    
    # master dictionary of samples containing another dictionary of sorted sample reads (containsExons vs not)
    readDict = defaultdict(dict) 
    
    for filename, bamReader in zip(bamFilenameList, bamReaderList):
        # Sort current bam's reads
        readsContainInterval, readsOverlapInterval, readStatsDict = sort_reads(bamReader, readIntervalRange)
        # Add it to the master dictionary
        readDict[filename] = {'readsContainInterval': readsContainInterval, 
                              'readsOverlapInterval': readsOverlapInterval,
                              'readStatsDict': readStatsDict}
        
    return readDict, readIntervalRange



##################################################
### Determine splicing patterns of annotations ###
##################################################

# Need to redesign script to remove the generation of a boolian matrix. It is
# an unnessisary step that doesn't scale well. A better method would be to 
# directly store the patterns found in the annotation and bam reads

def create_boolMatrix(bedList, bedNum):
    '''Create a binary matrix of all possible patterns of exons and convert to boolean matrix. 
    BoolMatrix is used determining splicing patterns of transcripts and reads'''
    binMatrix = [] # create a binary matrix as reference for the boolean matrix

    # Fill the binary matrix with binary numbers zero to exonNum-1
    # I call it a matrix, but its really a python list of lists
    for i in tqdm(range(2**bedNum)): # bedNum is the count of number of exons in bed file
        binNum = f'{i:0{bedNum}b}'
        binMatrix.append(binNum)
    
    # Create an boolean matrix of the same size as the binary one, but all 'False' entries
    boolMatrix=[]
    for i in tqdm(range(2**bedNum)):
        l = []
        for j in range(bedNum):
            l.append(False)
        boolMatrix.append(l) # Create boolean matrix full of False 

    # Now change the entries to match the binary matrix 1 = True, 0 = False
    for i in tqdm(range(2**bedNum)):
        for j in range(bedNum):
            if binMatrix[i][j] == '0':
                boolMatrix[i][j] = False
            elif binMatrix[i][j] == '1':
                boolMatrix[i][j] = True
    
    return boolMatrix # Lets return the boolMatrix too so we can use it for the transcript dictionary

def sort_trans(gffReader, intervalRange, addTransStart=0, addTransEnd=0):
    # Only select transcripts that contain the full range of all exons profided 
    # (start of ftransIntervalRangest exon and end of last exon
    transContainInterval = [] # contain not mean all exons match, just included in the interval
    transOverlapInterval = []
    
    # Adjust the interval range
    transIntervalRange = intervalRange.copy() # make a copy
    transIntervalRange.start += addTransStart
    transIntervalRange.end += addTransEnd
    
    print('Sorting transcripts...')
    for annot in tqdm(gffReader):
        if annot.type == 'transcript': # only look at annotated transcripts this time
            if annot.iv.contains(transIntervalRange):
                transContainInterval.append(annot)
            if annot.iv.overlaps(transIntervalRange):
                transOverlapInterval.append(annot)
            
    return transContainInterval, transOverlapInterval, transIntervalRange


def group_trans_exons(transContainInterval, gffReader):

    transDict = defaultdict(list)

    print('Grouping exons by transcript...')
    for annot in gffReader:
        if annot.type == 'exon': # Only look at exons
            for trans in transContainInterval:
                if trans.attr['transcript_id'] == annot.attr['transcript_id']:
                    transDict[trans.attr['transcript_id']].append(annot)
    
    return transDict


def check_trans_patterns(transDict, bedList, bedNum, boolMatrix):    
    
    # Dict similar to the readsBoolCountDict, except instead of fl_read counts
    # It will contain a list of transcript ids that share that splice pattern
    transBoolDict = {} 
    
    print('Checking transcript patterns...')
    for i in tqdm(range(2**bedNum)): # bedNum = number of exons in the bed file
        # Turn the list into a tuple, and assign it a '' string as default
        transBoolDict[tuple(boolMatrix[i])] = '' 
        
    # Now fill up the annotated ones with the transcript Ids
    for transId, transExonList in tqdm(transDict.items()):
        matchList = [False] * bedNum # Store matched exons (Default False)
        for i, bed in enumerate(bedList):
            for transExon in transExonList:
                if bed.iv.overlaps(transExon.iv):
                    matchList[i] = True

        matchTuple = tuple(matchList)
        if transBoolDict[matchTuple] == '':
            transBoolDict[matchTuple] = f'{transId}' # the first, simple assignment
        else:
            transBoolDict[matchTuple] += f', {transId}' # otherwise add a comma and space
        
#     # Fill up the non-annotated patterns
#     for boolKey, transIdList in transBoolDict.items():
#         if not transIdList:
#             transBoolDict[boolKey] += 'NaN'
    
    return transBoolDict


def gen_boolTransDf(transBoolDict, boolMatrix, bedList):
    boolTransDf = pd.DataFrame(boolMatrix) # turn boolMatrix into a dataframe
    
    # Get names of bed/exons from bedlist and assign as column names
    bedNames = []
    for bed in bedList:
        bedNames.append(bed.name)
    boolTransDf.columns = bedNames
    
    for i, (patternKey, transIdList) in enumerate(transBoolDict.items()):
        boolTransDf.at[i, 'transcript_ids'] = transIdList
    
    return boolTransDf


#######################################################
### Determine and quantify splice patterns in reads ###
#######################################################


def rm_indels_in_cigar(cigar_obj):
    """Pacbio data has many long-reads with insertion and deletion (InDel) mistakes. 
    In the sam alignment file, we want to get coordinates of whole exons in order to 
    determine if they are annotated or not.  However, InDels split up exon coordinates 
    ('M') in the cigar string creating what looks like multiple short exons between 
    introns when parsing with HTSeq.  Since we are only interested in the matched 
    coordinates, this function will merge all matched coordinates between introns ('N's)
    and return a list of genomic interval objects corresponding to the full exons"""
    new_cigar = []
    whole_exon = {'chr':'', 'start':'', 'end':'', 'strand':''}
    first_match_of_exon = True
    for c in cigar_obj:
        if (c.type == "M"):
            if first_match_of_exon: # add new coordinate information
                whole_exon['chr'] = c.ref_iv.chrom
                whole_exon['start'] = c.ref_iv.start
                whole_exon['end'] = c.ref_iv.end
                whole_exon['strand'] = c.ref_iv.strand
                first_match_of_exon = False
            else: # change only the end coordinate of exon (skips InDels)
                whole_exon['end'] = c.ref_iv.end

        if c.type == 'N': # This is an intron, signals end of exon
            new_giv = HTSeq.GenomicInterval(whole_exon['chr'], whole_exon['start'], 
                                            whole_exon['end'], whole_exon['strand'])
            new_cigar.append(new_giv)
            first_match_of_exon = True # reset to true for next exon

    # If not a first match, that means the cigar string likely had a final match with no
    # indels and/or ended on an intron, soft-clip or hard-clip 
    if first_match_of_exon == False: # Let's add the final match to the new cigar list
            new_giv = HTSeq.GenomicInterval(whole_exon['chr'], whole_exon['start'], 
                                            whole_exon['end'], whole_exon['strand'])
            new_cigar.append(new_giv)
            
    return new_cigar


def gen_empty_readBoolDict(boolMatrix, bedNum):
    readBoolDict = {}
    for i in range(2**bedNum): # bedNum = number of exons in the bed file
        # Turn the list into a tuple, and assign it a zero count as default
        readBoolDict[tuple(boolMatrix[i])] = 0 
    
    return readBoolDict

def count_bed_patterns(readDict, boolMatrix, bedList, bedNum):
    
    sampleBoolDict = {}
    for sampleKey, sampleDict in tqdm(readDict.items()):
        readsContainInterval = sampleDict['readsContainInterval']
        
        readBoolDict = gen_empty_readBoolDict(boolMatrix, bedNum)
        for read in tqdm(readsContainInterval):
            newCigar = rm_indels_in_cigar(read.cigar)
            matchList = [False] * bedNum # Store matched exons (Default False)
            for cig in newCigar:
                for i, bed in enumerate(bedList):
                    if bed.iv.overlaps(cig):
                        matchList[i] = True

            matchTuple = tuple(matchList)
            flCount = grab_fullLength_count(NOT_COLLAPSED_REGEXP, read)
            readBoolDict[matchTuple] += int(flCount)
        
        sampleBoolDict[sampleKey] = readBoolDict
    
    return sampleBoolDict


###########################################################################
### Combine and filter annotation and read data into a single dataframe ###
###########################################################################

def combine_trans_reads(boolTransDf, sampleBoolDict):
    
    sampColNameList = [] # keep track of these so we can use them to filter our dataframe
    for sampleKey, readBoolDict in tqdm(sampleBoolDict.items()):
        for i, (patternKey, flCount) in enumerate(readBoolDict.items()):
            sampColName = f'{sampleKey}_flCount'
            sampColNameList.append(sampColName)
            boolTransDf.at[i, sampColName] = flCount
    
    sampColNameSet = set(sampColNameList) # keep only uniques
    totalDf = boolTransDf # rename to avoid confusion
    
    return totalDf, sampColNameSet


def filter_df(totalDf, sampColNameSet):
    '''Function to filter the totalDf to include only rows that either have a transcript id or sample reads'''
    # Main query string that we will use for filtering
    queryString = "transcript_ids != ''"
    
    for sampColName in sampColNameSet:
        sampString = f"or {sampColName} > 0"
        queryString += sampString
    
    finalDf = totalDf.query(queryString)
    finalDf.index.name = 'pattern_num'
    
    return finalDf


######################################
### Read and transcript statistics ###
######################################
        

def calc_read_stats(readDict, intervalRange, readIntervalRange):
    '''Calculates sorted reads from each sample and returns a dataframe.'''
    dfDict = {}
    for sampleKey, sampleDict in tqdm(readDict.items()):
        # Reads that contain all exons
        readStatsDict = sampleDict['readStatsDict']
        containCluster = readStatsDict['readsContainInterval_clusterCount']
        containFl      = readStatsDict['readsContainInterval_flCount']
        overlapCluster = readStatsDict['readsOverlapInterval_clusterCount']
        overlapFl      = readStatsDict['readsOverlapInterval_flCount']
        allCluster     = readStatsDict['allReads_clusterCount']
        allFl          = readStatsDict['allReads_flCounts']
        
        dfDict[sampleKey] = {('reads_contain_interval','cluster_count'): containCluster,
                            ('reads_contain_interval','fl_count'): containFl,
                            ('reads_overlap_interval','cluster_count'): overlapCluster,
                            ('reads_overlap_interval','fl_count'): overlapFl,
                            ('all_reads', 'cluster_count'): allCluster,
                            ('all_reads', 'fl_count'): allFl }
        
        readStatsDf = pd.DataFrame(dfDict).T
        readStatsDf['original_interval_range'] = intervalRange
        readStatsDf['read_interval_range'] = readIntervalRange
        readStatsDf.index.name = 'samples'
        
    return readStatsDf

def trans_stats(transContainInterval, transOverlapInterval, intervalRange, transIntervalRange):
    containList = []
    overlapList = []
    containStr = ''
    overlapStr = ''
    droppedStr = ''
    
    for annot in transContainInterval:
        transId = annot.attr['transcript_id']
        containList.append(transId)
        if containStr == '':
            containStr = transId
        else:
            containStr += f', {transId}'
    
    for annot in transOverlapInterval:
        transId = annot.attr['transcript_id']
        overlapList.append(annot.attr['transcript_id'])
        if overlapStr == '':
            overlapStr = transId
        else:
            overlapStr += f', {transId}'
        
    # Get what transcripts that intersect, or ones that overlap, but not contain
    droppedTrans = list(set(overlapList)- set(containList))
    for transId in droppedTrans:
        if droppedStr == '':
            droppedStr = transId
        else:
            droppedStr += f', {transId}'
    
    transStatsDf = pd.DataFrame([containStr, overlapStr, droppedStr, intervalRange, transIntervalRange], 
              index=['transcripts_contain_interval_range', 'transcripts_overlap_interval_range',
                    'transcript_dropped', 'original_interval_range', 'transcript_interval_range'])
    
    transStatsDf.index.name = 'statistic'
    transStatsDf.columns = ['detail']

    return transStatsDf


##########################
### Output data to csv ###
##########################

def output_data(outputPrefix, finalDf, readStatsDf, transStatsDf):
    finalDf.to_csv(f'{outputPrefix}_main_analysis.csv')
    readStatsDf.to_csv(f'{outputPrefix}_read_statistics.csv')
    transStatsDf.to_csv(f'{outputPrefix}_annotation_statistics.csv')


##########################
### Parsing user input ###
##########################


def parse_input():
    """Use argparse to handle user input for program"""
    
    # Create a parser object
    parser = argparse.ArgumentParser(
        prog='exPhaser.py',
        
        formatter_class=argparse.RawDescriptionHelpFormatter,
        
        description="""Help
        
        Usage:

        exPhaser.py --bam <sample2.bam sample2.bam ... > --gtf <gene.gtf> --bed <gene_exons.bed> -outputPrefix <prefix> 

        Optional parameters (--addTransStart <#> --addTransEnd <#> --addReadStart <#> --addReadEnd <#>)


        Description:
        ExPhaser is a script that helps phases exons to determine if neighboring
        splicing events are happening on the same transcript or not. When using 
        an 'exon-based' analysis approach this phasing data is lost. Therefore,
        we can use this program to select key (cassette) exons that could help
        define unique transcripts in the reads and determine which annotations
        they belong to with partial transcript information.

        For example, if we look at the percent spliced in of two exons
        independently, we might see results like below:


        | ==A== |                | ==B== |
        | ==A== |                | ===== |
        | ===== |                | ==B== |
        | ===== |                | ===== |
        |_______|                |_______|
        |  50%  |                |  50%  |

        Observing the percent spliced in of Exon A and Exon B, it's 50%
        for both. This is rather ambiguous because its possible that when phased
        the splice pattern may look like below:

        ===========A=================B================
        ===========A==================================
        =============================B================
        ==============================================

        With phasing information, we can see that in the context of neighboring
        exons that there are four unique splice transcripts. ExPhaser will be
        able to help determine the annotation of each splice pattern if it exists
        and quantify the number of transcripts in the provided sample reads.
        """)

    parser.add_argument('--addTransStart', 
        help="Adjust intervalRange start coordinate for selecting transcrpt annotations", 
        type=int, default=0)

    parser.add_argument('--addTransEnd', 
        help="Adjust intervalRange end coordinate for selecting transcrpt annotations", 
        type=int, default=0)

    parser.add_argument('--addReadStart', 
        help="Adjust intervalRange start coordinate for selecting reads", 
        type=int, default=0)

    parser.add_argument('--addReadEnd', 
        help="Adjust intervalRange end coordinate for selecting reads", 
        type=int, default=0)

    parser.add_argument("-v", "--version", action="version",
                        version=textwrap.dedent("""\
        %(prog)s
        -----------------------   
        Version:    0.1 
        Updated:    02/01/2019
        By:         Prech Uapinyoying   
        Website:    https://github.com/puapinyoying"""))

    # Essential args
    requiredNamed = parser.add_argument_group('Required named arguments:')
    
    requiredNamed.add_argument('--bams', nargs='+',
        help="IsoSeq sample.bam file(s) to analyze, uncollapsed data (not processed through ToFu", required=True)
    
    requiredNamed.add_argument('--gtf',  
        help="""An annotation (GTF) file containing gene, transcript and exon features.
        It works best if you use a gtf from a single gene.
        e.g. grep 'gene_id' gencode_mm10.gtf > gene.gtf""", required=True)
    
    requiredNamed.add_argument('--bed', 
        help="""UCSC style tab delimited bed file contianing genomic coordinates for exons to analyze. 
        The exPhaser requires the bed file to have five columns of data: chromosome, start, end, name, 
        score and strand. Easiest way to make a bed file is to use USCS Table browser to export exon 
        features as bed and manually chose exons you would like to include in the analysis. The name 
        column in the bed file can be changed if desired. It is used as column headers for the csv 
        output.""", required=True)
    
    requiredNamed.add_argument('-o','--outputPrefix', 
        help="Specify prefix name for output files", required=True)

    args = parser.parse_args()
    
    return args


#############
#### Run ####
#############


if __name__ == "__main__":
    args = parse_input()

    # Parse user input
    gtfPath = args.gtf
    bedPath = args.bed
    bamPathList = args.bams
    outputPrefix = args.outputPrefix
    addTransStart = args.addTransStart
    addTransEnd = args.addTransEnd
    addReadStart = args.addReadStart
    addReadEnd = args.addReadEnd

    # Load data
    gffReader = load_gff(gtfPath)
    bedReader = load_bed(bedPath)
    bamReaderList, bamFilenameList = load_bam_list(bamPathList)

    # Calculate interval range based on bed exons
    bedList, bedNum = get_bed_list(bedReader)
    intervalRange = get_interval_range(bedList)

    # Sort for reads that contain and overlap the interval range
    print('Processing sample bam files...')
    readDict, readIntervalRange = sort_reads_from_bamReaderList(
        bamReaderList, bamFilenameList, intervalRange, addReadStart, addReadEnd)

    # Generate a true/false matrix of all possible splice patterns
    print("Generating matrix of all possible splice patterns...")
    boolMatrix = create_boolMatrix(bedList, bedNum)

    # Select transcritps that fall in the interval range and determine splice 
    # patterns of transcript annotations
    transContainInterval, transOverlapInterval, transIntervalRange = sort_trans(
        gffReader, intervalRange, addTransStart, addTransEnd)
    transDict = group_trans_exons(transContainInterval, gffReader)
    transBoolDict = check_trans_patterns(transDict, bedList, bedNum, boolMatrix)
    boolTransDf = gen_boolTransDf(transBoolDict, boolMatrix, bedList)

    # Counting bed/exon splice patterns in reads
    print('Counting splice patterns seen in reads...')
    sampleBoolDict = count_bed_patterns(readDict, boolMatrix, bedList, bedNum)

    # Combining splice pattern, transcript_ids and read counts together in a 
    # pandas dataframe & filter
    totalDf, sampColNameSet = combine_trans_reads(boolTransDf, sampleBoolDict)
    finalDf = filter_df(totalDf, sampColNameSet)

    # Calculating Stats
    print('Calculating statistics...')
    readStatsDf = calc_read_stats(readDict, intervalRange, readIntervalRange)
    transStatsDf = trans_stats(transContainInterval, transOverlapInterval, 
        intervalRange, transIntervalRange)

    # Output all the data to csv
    output_data(outputPrefix, finalDf, readStatsDf, transStatsDf)
    print('Complete!')