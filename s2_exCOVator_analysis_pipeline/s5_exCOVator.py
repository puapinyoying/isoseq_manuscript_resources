#! /usr/bin/env python

from __future__ import division
import os
import re
import sys
import csv
import HTSeq
import argparse
import textwrap
import pyensembl
import itertools
import pandas as pd
from tqdm import tqdm
from operator import attrgetter
from pandas import Series, DataFrame
from collections import Counter, defaultdict
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

plt.style.use('fivethirtyeight')

# Regular expression patterns for parsing SAM alignment read names
COLLAPSED_REGEXP = r'[PB\.\d]+\|[\w\:\(\)\-\+]+\|c\d+\/f(\d+)p(\d+)\/\d+'
NOT_COLLAPSED_REGEXP = r'c\d+\/f(\d+)p(\d+)\/\d+'

def check_file_exist(file_path):
    """Make sure the file exists"""
    if not os.path.isfile(file_path):
        file_name = os.path.basename(file_path)
        print "{0} was not found. Check path".format(file_name)
        print "Full path to file: {0}".format(file_path)
        sys.exit(1)

def load_ensembl_data(release, species):
    """Load the ensembl data so we can convert gene_ids into gene_names for more
    intuitive output"""
    try:
        ensembl_data = pyensembl.EnsemblRelease(release=release, 
            species=species)

    except ValueError as error_details:
        print error_details
        print "Make sure you have the pyensembl library and annotation release \
you want to use downloaded and properly installed. e.g.:"
        print "$ pip install pyensembl"
        print "$ pyensembl install --release 85 --species mouse"
        sys.exit(1)

    return ensembl_data

def get_gene_name(gene_id, ensembl_data):
    """Get the gene name from the ensembl gene id (gencode same as ensembl)"""
    base_id = gene_id.split('.')[0]
    return ensembl_data.gene_name_of_gene_id(base_id)


def load_sam(sam_path):
    check_file_exist(sam_path)
    sam_reader = HTSeq.SAM_Reader(sam_path)

    return sam_reader


def load_bam(bam_path):
    check_file_exist(bam_path)
    sam_reader = HTSeq.BAM_Reader(bam_path) # also call it a sam reader for ease

    return sam_reader


def load_gff(gff_path):
    check_file_exist(gff_path)
    gff_reader = HTSeq.GFF_Reader(gff_path, end_included=True)
    return gff_reader

def get_file_name(file_path):
    base_name = os.path.basename(file_path)
    file_name, file_extention = os.path.splitext(base_name)

    return file_name

def read_input_files(args):
    """Reads in the BAM/SAM file and GFF file and returns reader objects. The
    alignments will be put into a list of sam_reader objects, while the
    gff_reader is alone"""
    sam_reader_dict = {}

    if args.bam:
        for bam in args.bam:
            sam_reader = load_bam(bam)
            file_name = get_file_name(bam)
            sam_reader_dict[file_name] = sam_reader

    elif args.sam:
        for sam in args.sam:
            sam_reader = load_sam(sam)
            file_name = get_file_name(sam)
            sam_reader_dict[file_name] = sam_reader

    else:
        print "No BAM or SAM file given."
        sys.exit(1)

    if args.gff:
        gff_reader = load_gff(args.gff)

    return sam_reader_dict, gff_reader


def create_annot_dicts(gff_reader, ensemble_data):
    """dexseq_prepare_annotation.py"""
    gene_annot_dict = defaultdict(list)
    exon_annot_dict = defaultdict(list)
    
    #print "Creating GFF dictionary...",
    for feature in tqdm(gff_reader, "Creating GFF dictionary"):
        gene_name = get_gene_name(feature.attr["gene_id"], ensemble_data)
        if feature.type == 'aggregate_gene':
            gene_annot_dict[gene_name].append(feature)
        elif feature.type == 'exonic_part':
            exon_annot_dict[gene_name].append(feature)
        else:
            print "what is {0}".format(feature.type)

    return gene_annot_dict, exon_annot_dict


def create_sam_aln_dict(sam_reader_dict, gene_annot_dict):
    """Generate a dictionary to organize all the sam_alignments by gene_name.
    This will help when checking for coverage and whenlooking for unannotated 
    exons. Could also be used to check whole gene coverage too, but seems less
    useful because you could easily do it using a one liner with HTSeq_count.py"""
    sam_dict = defaultdict(lambda: defaultdict(set))

    for sam_name, sam_reader in sam_reader_dict.iteritems():
        #gene_aln_count = Counter() # count all reads alned to gene
    
        for aln in tqdm(sam_reader, "Parsing SAM alignments for '{0}'".format(sam_name)):
            if aln.aligned: # Make sure that the read is aligned to genome  
                for gene_name, feature in gene_annot_dict.iteritems():
                    if aln.iv.overlaps(feature[0].iv):
                        sam_dict[sam_name][gene_name].add(aln)
                        #gene_aln_count[gene_name] += 1
    
    return sam_dict #, gene_aln_count


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


def collapse_unannotated_exons(stranded, unannot_exon_dict):
    """Collapses the list of unannotated exon genomic itervals. Similar to
    what dexseq_prepare_annotation.py does, it takes all overlapping exons and
    collapses them down into a single unique coordinate that contains the rest.
    returns a set of unqiue genomic interval objects"""
     
    ga_dict = {} # genomic array dictionary
    for gene_name, exon_list in unannot_exon_dict.iteritems():
        # Create a GenomicArray object to get largest exons
        ga = HTSeq.GenomicArray('auto', stranded=stranded)
        
        for exon_iv in exon_list:
            ga[exon_iv] = 1 # real exons are set to 1, inbetween are 0's
        
        ga_dict[gene_name] = ga

    # create new defaultdict set to hold unique exon coordinates
    collapsed_unannot_dict = defaultdict(set)
    
    for gene_name, ga in ga_dict.iteritems():      
        for exon_iv, num in ga.steps():
            if num == 1:
                collapsed_unannot_dict[gene_name].add(exon_iv)
    
    return collapsed_unannot_dict


def find_unannotated(sam_dict, exon_annot_dict):
    # many alignments may have same unannotated exon, use 'set' to remove dups
    unannotated_exon_dict = defaultdict(set)
    
    stranded = True
    for sample_name, aln_dict in tqdm(sam_dict.iteritems(), 
        "Finding unannotated exons in alignment files"):
        for gene_name, sam_alignments in aln_dict.iteritems():
            for alignment in sam_alignments:
                # check strandedness for collapse function
                if alignment.iv.strand == '.': 
                    stranded = False
                aln_exon_list = rm_indels_in_cigar(alignment.cigar)
                for aln_exon in aln_exon_list:
                    match = False
                    for annot_exon_feature in exon_annot_dict[gene_name]:
                        if annot_exon_feature.iv.overlaps(aln_exon):
                            match = True
                            break
                    if match == False:
                        unannotated_exon_dict[gene_name].add(aln_exon)
            
    collapsed_unannot_dict = collapse_unannotated_exons(stranded, 
        unannotated_exon_dict)
        
    return collapsed_unannot_dict


def read_name_pattern(collapsed):
    '''Small script to choose regular expression patters for format of pacbio 
    read names'''
    if collapsed:
        name_pattern = COLLAPSED_REGEXP
    else:
        name_pattern = NOT_COLLAPSED_REGEXP
        
    return name_pattern


def get_unannot_coverage(sam_dict, unannot_exon_dict, collapsed):
    name_pattern = read_name_pattern(collapsed)
      
    unannot_cov_dict = {}
    for sam_name, aln_dict in sam_dict.iteritems():
        gene_cov_dict = {}
        for gene_name, unannot_exon_list in tqdm(unannot_exon_dict.iteritems(), 
                "Calculating unannotated exon coverage for {0}".format(sam_name)):
            
            # Create counters for unannotated exons, will use genomic interval as key
            unannot_ctotal = Counter()
            unannot_cmatch = Counter() 
            unannot_ftotal = Counter()
            unannot_fmatch = Counter()
            
            # if no alignments mapped to gene for this file
            if aln_dict[gene_name] == set([]):
                for unannot_exon in unannot_exon_list:       
                    unannot_ctotal[unannot_exon] += 0
                    unannot_ftotal[unannot_exon] += 0
                    unannot_cmatch[unannot_exon] += 0
                    unannot_fmatch[unannot_exon] += 0
                    
                gene_cov_dict[gene_name] = {'unannot_ctotal' : unannot_ctotal,
                                            'unannot_cmatch' : unannot_cmatch,
                                            'unannot_ftotal' : unannot_ftotal,
                                            'unannot_fmatch' : unannot_fmatch }

            else: # otherwise check individual exon alignments   
                for alignment in aln_dict[gene_name]:
                    # grab number of full-length reads from read.name
                    result = re.search(name_pattern, alignment.read.name)
                    num_FL_reads = int(result.group(1))

                    # Indels break up the exon coordinates, remove them
                    aln_exons_rm_indels = rm_indels_in_cigar(alignment.cigar)

                    for unannot_exon in unannot_exon_list:
                    # check total read coverage the unannotated exon region (does not need to match)                
                        if alignment.iv.overlaps(unannot_exon): 
                            unannot_ctotal[unannot_exon] += 1
                            unannot_ftotal[unannot_exon] += num_FL_reads

                            for aln_exon in aln_exons_rm_indels: 
                                # Check if the read sequence matches the exon
                                if aln_exon.overlaps(unannot_exon):
                                    unannot_cmatch[unannot_exon] += 1
                                    unannot_fmatch[unannot_exon] += num_FL_reads
                                    break # stop after first match
                        else:
                            unannot_ctotal[unannot_exon] += 0
                            unannot_ftotal[unannot_exon] += 0
                            unannot_cmatch[unannot_exon] += 0
                            unannot_fmatch[unannot_exon] += 0

                gene_cov_dict[gene_name] = {'unannot_ctotal' : unannot_ctotal,
                                            'unannot_cmatch' : unannot_cmatch,
                                            'unannot_ftotal' : unannot_ftotal,
                                            'unannot_fmatch' : unannot_fmatch }

        unannot_cov_dict[sam_name] = gene_cov_dict

    return unannot_cov_dict

    
def get_annot_coverage(sam_dict, exon_annot_dict, collapsed):
    """Obtains coverage data for each of the individual exonic coordinates given
    in the gff file. The function will calculate and return a dictionary object
    that contains the gene_name as the key and a dictionary of 4 counter 
    objects as values. """

    name_pattern = read_name_pattern(collapsed)
    
    # Dictionary of sam_name: coverage_dict
    annot_cov_dict = {}

    for sam_name, gene_aln_dict in sam_dict.iteritems():

        gene_cov_dict = {}
        # first layer, for each gene get the list of all annot exon features
        for gene_name, exon_feature_list in tqdm(exon_annot_dict.iteritems(), 
                "Calculating annotated exon coverage for {0}".format(sam_name)):

            # Create counters for annotated exons will use feature as key
            ctotal = Counter() #  cluster total counts for specific exon
            cmatch = Counter() # these reads need to actually match
            ftotal = Counter() # break cluster into full_length reads
            fmatch = Counter()
            
            # if no alignments mapped to gene for this file
            if gene_aln_dict[gene_name] == set([]):
                for exon_feature in exon_feature_list:       
                    ctotal[exon_feature] += 0
                    ftotal[exon_feature] += 0
                    cmatch[exon_feature] += 0
                    fmatch[exon_feature] += 0
                    
                gene_cov_dict[gene_name] = {'ctotal': ctotal, 
                                            'cmatch': cmatch,
                                            'ftotal': ftotal,
                                            'fmatch': fmatch}

            else: # otherwise check individual exon alignments
                for alignment in gene_aln_dict[gene_name]:
                    for exon_feature in exon_feature_list:
                        # Get number of full-length reads from read.name
                        result = re.search(name_pattern, alignment.read.name)
                        num_FL_reads = int(result.group(1))

                        # Indels break up the exon coordinates, remove them
                        aln_exons_rm_indels = rm_indels_in_cigar(alignment.cigar)

                         # get total coverage of annotated exon
                        if alignment.iv.overlaps(exon_feature.iv):
                            ctotal[exon_feature] += 1
                            ftotal[exon_feature] += num_FL_reads

                            for aln_exon in aln_exons_rm_indels: # exons from sam alignment
                                if aln_exon.overlaps(exon_feature.iv):
                                    cmatch[exon_feature] += 1
                                    fmatch[exon_feature] += num_FL_reads
                                    break # quit counting after first match
                        else:
                            ctotal[exon_feature] += 0
                            ftotal[exon_feature] += 0
                            cmatch[exon_feature] += 0
                            fmatch[exon_feature] += 0

                gene_cov_dict[gene_name] = {'ctotal': ctotal, 
                                            'cmatch': cmatch,
                                            'ftotal': ftotal,
                                            'fmatch': fmatch}

        annot_cov_dict[sam_name] = gene_cov_dict

    return annot_cov_dict


def merge_coverage_dicts(annot_cov_dict, unannot_cov_dict, output_prefix):
    """Function used to merge the annotated and unannotated coverage 
    dictionaries. It is the first step used to make it easier to convert the  
    count data into a pandas data frame. Outputs list of genes with no
    unannotated exons found (in the set of provided samples)"""
    merged_cov_dict = defaultdict(dict)
    no_unannot = []
    for sam_name, gene_cov_dict in annot_cov_dict.iteritems():
        for gene_name, annot_count_dict in gene_cov_dict.iteritems():
            try:
                unannot_count_dict = unannot_cov_dict[sam_name][gene_name]
                merged_count_dict = annot_count_dict.copy()        
                merged_count_dict.update(unannot_count_dict)
                merged_cov_dict[sam_name][gene_name] = merged_count_dict
                
            except KeyError as error_details:
                # Create an empty dummy entry, so dictionaries will be balanced
                unannot_fmatch = Counter()
                unannot_cmatch = Counter()
                unannot_ftotal = Counter()
                unannot_ctotal = Counter()
                unannot_count_dict = {'unannot_fmatch': unannot_fmatch,
                                      'unannot_cmatch': unannot_cmatch,
                                      'unannot_ftotal': unannot_ftotal,
                                      'unannot_ctotal': unannot_ctotal}
                merged_count_dict = annot_count_dict.copy()        
                merged_count_dict.update(unannot_count_dict)
                merged_cov_dict[sam_name][gene_name] = merged_count_dict

                no_unannot.append(gene_name)
                #print "No unannotated exons found for gene: '{0}'".format(gene_name)

    with open(output_prefix + '_no_unannot.csv', 'wb') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"')
        for gene_name in no_unannot:
            csv_writer.writerow([gene_name])

    print "No unannotated exons found for {0} genes".format(len(no_unannot))
                
    return merged_cov_dict    


def sort_unnannot_exon_list(unannot_exon_list):
    unannot_sorted_list = sorted(unannot_exon_list, key=attrgetter('start'))
    
    return unannot_sorted_list


def prep_dict_for_df(merged_cov_dict):
    sam_results_dict = {} # upper level dictionary using sample name as key
    for sam_name, coverage_dict in merged_cov_dict.iteritems():
        # dictionary for holding coverage results organized by gene
        results_dict = defaultdict(dict)

        for gene_name, all_count_dict in coverage_dict.iteritems():

            annot_feature_list = all_count_dict['ctotal'].keys()
            unannot_exon_list = all_count_dict['unannot_ctotal'].keys()

            # Run through annotated exons first
            for annot_feature in annot_feature_list:
                annotated = True
                exon_num = annot_feature.attr['exonic_part_number']
                chrom = annot_feature.iv.chrom
                start = annot_feature.iv.start
                end = annot_feature.iv.end
                location = '{0}:{1}-{2}'.format(chrom, start, end)
                strand = annot_feature.iv.strand
                length = annot_feature.iv.length

                ctotal  = all_count_dict['ctotal'][annot_feature]
                cmatch = all_count_dict['cmatch'][annot_feature]
                ftotal  = all_count_dict['ftotal'][annot_feature]
                fmatch = all_count_dict['fmatch'][annot_feature]
                
                if ctotal == 0: # prevent div by zero errors
                    cratio = 'NA'
                else:
                    cratio = cmatch / ctotal
                    
                if ftotal == 0:
                    fratio = 'NA'
                else:
                    fratio = fmatch / ftotal

                results_dict[gene_name][exon_num] = { 'annotated': annotated,
                    'location': location, 'chrom': chrom, 'start': start, 
                    'end': end, 'strand': strand, 'length': length, 
                    'ctotal': ctotal, 'cmatch': cmatch, 'cratio': cratio,
                    'ftotal': ftotal, 'fmatch': fmatch, 'fratio': fratio }
                
            # now add the unannotated exons
            # sort the unannot exons by start so exon number will be in order
            sorted_unannot_list = sort_unnannot_exon_list(unannot_exon_list)
            exon_count = 1
            for unannot_exon in sorted_unannot_list:
                annotated = False
                # exonic part numbering will continue -- after annotated ones
                # so that the indexing is simpler for ploting
                exon_num = 'u-{:03d}'.format(exon_count)
                chrom = unannot_exon.chrom
                start = unannot_exon.start
                end = unannot_exon.end
                location = '{0}:{1}-{2}'.format(chrom, start, end)
                strand = unannot_exon.strand
                length = unannot_exon.length

                ctotal  = all_count_dict['unannot_ctotal'][unannot_exon]
                cmatch = all_count_dict['unannot_cmatch'][unannot_exon]
                ftotal  = all_count_dict['unannot_ftotal'][unannot_exon]
                fmatch = all_count_dict['unannot_fmatch'][unannot_exon]
                
                if ctotal == 0: # prevent div by zero errors
                    cratio = 'NaN'
                else:
                    cratio = cmatch / ctotal
                    
                if ftotal == 0:
                    fratio = 'NaN'
                else:
                    fratio = fmatch / ftotal

                results_dict[gene_name][exon_num] = { 'annotated': annotated,
                    'location': location, 'chrom': chrom, 'start': start, 
                    'end': end, 'strand': strand, 'length': length, 
                    'ctotal': ctotal, 'cmatch': cmatch, 'cratio': cratio,
                    'ftotal': ftotal, 'fmatch': fmatch, 'fratio': fratio }

                exon_count += 1

        sam_results_dict[sam_name] = results_dict

    return sam_results_dict


def results_dict_to_df(sam_results_dict):
    
    df_list = [] # list of dataframes for building the large one
    final_df_list = [] # sorted and sliced, ready to merge
    
    # count and iterate though the results dictionary
    for count, (sam_name, results_dict) in enumerate(sam_results_dict.iteritems()):
        # format the dictionary into a more intuitive dataframe structure 
        plain_df = DataFrame.from_dict({(gene, exon_num): results_dict[gene][exon_num] 
                           for gene in results_dict.keys()
                           for exon_num in results_dict[gene].keys()},
                        orient='index')

        # sort the columns
        sorted_df = plain_df.ix[:,['location', 'chrom', 'start', 'end', 'strand', 
                                   'length', 'annotated','ctotal', 'ftotal', 'cmatch', 
                                   'fmatch', 'cratio', 'fratio']]

        # sort the columns
        sorted_df = plain_df.ix[:,['location', 'chrom', 'start', 'end', 'strand', 
                                   'length', 'annotated','ctotal', 'cmatch', 
                                   'cratio','ftotal', 'fmatch',  'fratio']]
        
        # Add an extra header level to describe the sample name
        named_df = pd.concat({sam_name: pd.DataFrame(sorted_df)}, axis=1)
        df_list.append(named_df)
    
    # sort the data frames by sample names
    sorted_df_list = sorted(df_list, key=attrgetter('columns.levels'))
    for count, named_df in enumerate(sorted_df_list):
        if count > 0:
            sliced_df = named_df.ix[:,7:13]
            final_df_list.append(sliced_df)            
        else:
            final_df_list.append(named_df)
    
    # merge df together
    merged_df = pd.concat(final_df_list,axis=1)
    merged_df.index.names = ['gene_name', 'exonic_part_num']
    merged_df.columns.names = ['sample', 'exon_info']

    return merged_df


def get_col_indexes(sample_names, ratio_type):
    num_samples = len(sample_names)
    
    if ratio_type == 'fratio':
        ratio_col = 12 # col index fratio (full-length ratio) of first sample
        total_col = 10
    elif ratio_type == 'cratio':
        ratio_col = 9
        total_col = 7
        
    ratio_col_indexes = [ratio_col]
    total_col_indexes = [total_col]
    
    if num_samples > 1:
        for i in range(0, num_samples-1):
            ratio_col += 6
            total_col  += 6
            ratio_col_indexes.append(ratio_col)
            total_col_indexes.append(total_col)
    
    return ratio_col_indexes, total_col_indexes


def get_fig_params(ratio_df2):
    exon_num = len(ratio_df2)
    if exon_num < 20:
        xticks_size = 10
        fig_width = 7

    elif exon_num >= 20 and exon_num < 50:
        xticks_size = 10
        fig_width = 8

    elif exon_num >= 50 and exon_num < 100:
        xticks_size = 9
        fig_width = 10

    elif exon_num >= 100 and exon_num < 150:
        xticks_size = 9
        fig_width = 15

    elif exon_num >= 150 and exon_num < 200:
        xticks_size = 8
        fig_width = 20

    elif exon_num >= 200 and exon_num < 300:
        xticks_size = 7
        fig_width = 25

    else:
        xticks_size = 6
        fig_width = 40
        
    return xticks_size, fig_width


def plot_ratio_by_gene(merged_df, ratio_type='fratio', output_name='no_name'):
    
    # Set graph variables based on ratio_type user supplies
    if ratio_type == 'fratio':
        total_type = 'ftotal'
        line_ylabel = "FL Match / Total"
        bar_ylabel = "Total FL Reads"
    elif ratio_type == 'cratio':
        total_type = 'ctotal'
        line_ylabel = "Custer Match / Total"
        bar_ylabel = "Total Cluster Reads"
     
    with PdfPages(output_name + '.pdf') as pdf:
        for gene_name in tqdm(merged_df.index.levels[0], "Generating gene plots"):  

            # build a sub dataframe containing full-length or cluster data
            sample_names = merged_df.columns.levels[0]
            num_samp =len(sample_names)     
            ratio_col_indexes, total_col_indexes = get_col_indexes(
                sample_names, ratio_type)
            
            # sort so unannotated exons will be in order by start coordinate 
            sort_start_df = merged_df.ix[gene_name].sort_values(
                by=[(sample_names[0],'start')],ascending=[1])

            #print sort_start_df
            # grab the columns for data of interest
            ratio_df = sort_start_df.ix[:,ratio_col_indexes] # full-length or cluster ratio
            total_df  = sort_start_df.ix[:,total_col_indexes] # total reads, for bottom graph

            # reset the index, so we can make this new exon order perminant
            ratio_df2 = ratio_df.reset_index()
            #print ratio_df2
            total_df2 = total_df.reset_index()
            #print total_df2

            xticks_size, fig_width = get_fig_params(ratio_df2)
            colors = ["#f1595f","#599ad3","#79c36a","#f9a65a","#9e66ab","#727272",
                      "#cd7058","#d77fb3"]

            fig, axarr = plt.subplots(2, sharex=True)
            fig.set_size_inches(fig_width,5)
            
            prev_total_df = 0
            for num, sample in enumerate(sample_names):
                samp_ratio_df = ratio_df2[sample]
                samp_total_df = total_df2[sample]

                # Make sure to remove any NaN values in the fratio column
                samp_ratio_df2 = pd.to_numeric(samp_ratio_df[ratio_type], errors='coerce')
                samp_total_df2 = pd.to_numeric(samp_total_df[total_type], errors='coerce')

                # Make the top line graph
                axarr[0].set_ylim([-0.2,1.2])
                axarr[0].plot(samp_ratio_df2.index, samp_ratio_df2, linestyle='-', marker='s',
                              linewidth=3,alpha=0.5, color=colors[num], label=sample)
                axarr[0].set_title(gene_name)
                axarr[0].set_ylabel(line_ylabel, fontsize=10)

                # Make the bottom barplot
                axarr[1].set_xlabel("Exonic Part Number", fontsize=11)
                axarr[1].set_ylabel(bar_ylabel, fontsize=10)
                if num == 0:
                    axarr[1].bar(samp_total_df2.index, samp_total_df2, align='center', alpha=0.5, 
                                 color=colors[num], label=sample)
    
                else:
                    axarr[1].bar(samp_total_df2.index, samp_total_df2, align='center',  alpha=0.5,
                                 color=colors[num],label = sample, bottom = prev_total_df)
                    
                legend = axarr[1].legend(title='Samples', loc='best', shadow=True, prop={'size':8})
                plt.setp(legend.get_title(),fontsize='xx-small')

                prev_total_df += samp_total_df2

            plt.xticks(ratio_df2.index, ratio_df2['exonic_part_num'], rotation='vertical',size=xticks_size)
            plt.setp(axarr[1].get_yticklabels(),fontsize=11)
            plt.setp(axarr[0].get_yticklabels(),fontsize=11)
            plt.xlim([-1,len(ratio_df2.index)])
            plt.tight_layout()
            pdf.savefig(fig)  # saves the current figure into a pdf page
            plt.close()


def df_to_csv(df, output_prefix):
    output_file = output_prefix + '.csv'
    df.to_csv(output_file, sep=',', encoding='utf-8')


def parse_input():
    """Use argparse to handle user input for program"""
    
    # Create a parser object
    parser = argparse.ArgumentParser(
        prog='exonCoverage.py',
        
        formatter_class=argparse.RawDescriptionHelpFormatter,
        
        description="""Script that calculates the coverage for each exon of the
        genes in the given gff file. The number of matched alignments are
        normalized by the number of total alignments that cover each specific 
        exon. Please use 'dexseq_prepare_annotation script on the gencode gtf 
        file first. Preferably, filter the original gtf file for a few genes of 
        interest or it may take hours.

        The script will be able to look for unannotated exons that show up in 
        sam alignments too.""")

    parser.add_argument('-c', '--collapsed', help="", action='store_true')
    parser.add_argument('--ratio_type', help="", default='fratio')

    parser.add_argument("-v", "--version", action="version",
                        version=textwrap.dedent("""\
        %(prog)s
        -----------------------   
        Version:    0.2 
        Updated:    01/20/2017
        By:         Prech Uapinyoying   
        Website:    https://github.com/puapinyoying"""))

    # Essential args
    requiredNamed = parser.add_argument_group('Required named arguments:')
    requiredNamed.add_argument('--gff', 
        help="Dexseq processed gff file", required=True)
    requiredNamed.add_argument('--species',  
        help="Specify annotation species (e.g. mouse or human)", required=True)
    requiredNamed.add_argument('--release', 
        help="Specify ensembl release number", type=int, required=True)
    requiredNamed.add_argument('-o','--output_prefix', 
        help="Specify prefix name for output file", required=True)


    # One or the other
    sam_or_bam = parser.add_mutually_exclusive_group(required=True)
    sam_or_bam.add_argument('--bam', nargs='+', help="")
    sam_or_bam.add_argument('--sam', nargs='+', help="")

    args = parser.parse_args()
    
    return args

#############
#### Run ####
#############

if __name__ == "__main__":
    args = parse_input()

    # Load ensembl data to get gene_names from gene_names
    #print "Loading Ensembl data.."
    ensembl_data = load_ensembl_data(args.release, args.species)

    # Load input files and generate reader objects for accessing them
    #print "\nLoading alignment and gff file..."
    sam_reader_dict, gff_reader = read_input_files(args)

    # Create an annotation dictionary of each whole gene and the exons
    #print "Creating annotation dictionaries: "
    gene_annot_dict, exon_annot_dict = create_annot_dicts(gff_reader, ensembl_data)
    
    # Create a sam alignment dictionary and a dictionary of count data for all
    # alignments that map to each gene
    sam_dict = create_sam_aln_dict(sam_reader_dict, gene_annot_dict)
    
    # Create a list of genomic_intervals for unannotated exons (not found in
    # gff, but in the sam alignments
    unannot_exon_dict = find_unannotated(sam_dict, exon_annot_dict)

    # Get exon coverage for calculating the ratios
    unannot_cov_dict = get_unannot_coverage(sam_dict, unannot_exon_dict, args.collapsed)
    annot_cov_dict = get_annot_coverage(sam_dict, exon_annot_dict, args.collapsed)

    # Merge the annotated and unannotated coverage dictionaries
    merged_cov_dict = merge_coverage_dicts(annot_cov_dict, unannot_cov_dict, 
        args.output_prefix)
    
    # Create a new dictionary with calculated results that is a bit easier
    # to convert into pandas DataFrames
    sam_results_dict = prep_dict_for_df(merged_cov_dict)

    # Convert the results dictionary into a large pandas DataFrame with all samples
    merged_df = results_dict_to_df(sam_results_dict)

    # Write the DataFrame of results to csv
    df_to_csv(merged_df, args.output_prefix)

    # generate plots for each gene and output to pdf file
    plot_ratio_by_gene(merged_df, args.ratio_type, args.output_prefix)
    
    print "Output results to: {0}.csv and {0}.pdf".format(args.output_prefix)


