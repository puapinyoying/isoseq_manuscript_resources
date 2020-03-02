from __future__ import division
import pandas as pd
import argparse
import textwrap
import os
import sys
import csv
from tqdm import tqdm

def check_file_exist(file_path):
    """Make sure the file exists"""
    if not os.path.isfile(file_path):
        file_name = os.path.basename(file_path)
        print "{0} was not found. Check path".format(file_name)
        print "Full path to file: {0}".format(file_path)
        sys.exit(1)

def import_df_from_csv(file_name):
# Import the csv and turn it into the exact same multi-index, multi-column, dataframe as merge_df
    df = pd.read_csv(file_name,  index_col=[0,1], 
                 header=[0,1],tupleize_cols=False, skipinitialspace=True)
    return df

def get_genes_and_samples(df):
    gene_list = df.index.get_level_values(0).unique()
    sample_list = df.columns.get_level_values(0).unique()
    
    return gene_list, sample_list

def filter_exons(df, gene_list, sample_list, min_cov, cov_readtype, ratio_diff):
    """Filter exons by minimum total coverage and ratio of differential expression between samples"""
    cov_filter_df_list = [] # Store the filtered gene dataframe to create a new large data frame with results
    full_filter_df_list = []
    
    for gene in tqdm(gene_list):
        df1 = df.loc[gene] # grab a slice of the data frame with a single gene
    
        if len(sample_list) == 3: # If there are three samples in csv file
            s1 = sample_list[0] # simplify the sample names into sort variables
            s2 = sample_list[1]
            s3 = sample_list[2]
            
            # Filter exons by  coverage first
            f1 = df1[(df1[s1][cov_readtype] >= min_cov) | \
                    (df1[s2][cov_readtype] >= min_cov) | \
                    (df1[s3][cov_readtype] >= min_cov)]
            
            # Rm exons that are 100% PSI in all 3 samples
            f2 = f1[(f1[s1]['fratio'] < 1) | (f1[s2]['fratio'] < 1) | (f1[s3]['fratio'] < 1)]
            
            m1 = f2[s1]['ftotal'].mean() # Calc the mean of total fl read coverage for each exon
            m2 = f2[s2]['ftotal'].mean()
            m3 = f2[s3]['ftotal'].mean()
            
            # Compare ratios (0-1) for each exon in gene and compare ratio by a given threshold
            if m1 <= 0 and m2 > 0 and m3 > 0:
                # Get the absolute value of the difference & compare to ratio_diff
                f3 = f2[(abs(f2[s2]['fratio'] - f2[s3]['fratio']) >= ratio_diff)]

            elif m2 <= 0 and m1 > 0 and m3 > 0:
                f3 = f2[(abs(f2[s1]['fratio'] - f2[s3]['fratio']) >= ratio_diff)]

            elif m3 <= 0 and m1 > 0 and m2 > 0:
                f3 = f2[(abs(f2[s1]['fratio'] - f2[s2]['fratio']) >= ratio_diff)]

            else: 
                f3 = f2[ (abs(f2[s1]['fratio'] - f2[s2]['fratio']) >= ratio_diff) | \
                         (abs(f2[s1]['fratio'] - f2[s3]['fratio']) >= ratio_diff) | \
                         (abs(f2[s2]['fratio'] - f2[s3]['fratio']) >= ratio_diff) ]
        
        # Its much more straight forward with two samples
        elif len(sample_list) == 2:
            s1 = sample_list[0]
            s2 = sample_list[1]
            f1 = df1[(df1[s1][cov_readtype] >= min_cov) | \
                     (df1[s2][cov_readtype] >= min_cov)]
            f2 = f1[(f1[s1]['fratio'] < 1) | (f1[s2]['fratio'] < 1)]
            m1 = f2[s1]['ftotal'].mean()
            m2 = f2[s2]['ftotal'].mean()
            f3 = f2[(abs(f2[s1]['fratio'] - f2[s2]['fratio']) >= ratio_diff)]
            
        else:
            print "Script does not support > 3 or less than 2 samples"
            sys.exit(1)
            
        if len(f3) > 0: # If its not empty
            f1_copy = f1.copy()
            f1_copy['gene_name'] = gene
            f1_copy.set_index('gene_name', append=True,  inplace=True) # Set it as the index
            
            f3_copy = f3.copy() # make a copy of the slice
            f3_copy['gene_name'] = gene # Put the gene name back on the dataframe
            f3_copy.set_index('gene_name', append=True,  inplace=True) # Set it as the index
            
            # Reorder it to original position in front
            cov_filter_gene_df = f1_copy.reorder_levels(['gene_name', 'exonic_part_num'])
            full_filter_gene_df = f3_copy.reorder_levels(['gene_name', 'exonic_part_num'])
            
            # Sort the dfs so that unannotated exons will be listed in order by start coordinate
            sorted_cov_filter_df = cov_filter_gene_df.sort_values([(sample_list[0], 'start')], ascending=True)
            sorted_full_filter_df = full_filter_gene_df.sort_values([(sample_list[0], 'start')], ascending=True)
            cov_filter_df_list.append(sorted_cov_filter_df)
            full_filter_df_list.append(sorted_full_filter_df)
    
    cov_filter_df = pd.concat(cov_filter_df_list)
    full_filter_df = pd.concat(full_filter_df_list)
    
    return cov_filter_df, full_filter_df

def output_gene_list_from_df(df, output_prefix):  
    output_name = output_prefix + '_gene_list.csv'
    with open(output_name, 'wb') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"')
        gene_list = df.index.get_level_values(0).unique()
        for gene in gene_list:
            csv_writer.writerow([gene])

def df_to_csv(df, output_prefix):
    output_file = output_prefix + '.csv'
    df.to_csv(output_file, sep=',', encoding='utf-8')


def parse_input():
    """Use argparse to handle user input for program"""
    # Create a parser object
    parser = argparse.ArgumentParser(
        prog='s6_filterDiffUsedExons.py',
        
        formatter_class=argparse.RawDescriptionHelpFormatter,
        
        description="""Script to filter for exons by minimum coverage and ones 
        that are differentially used. Options to filter by consensus reads (default)
        or full-length reads), filter by minimum number of reads, and minimum 
        difference between samples as a ratio (e.g. ratio_diff=0.10

        example:""")
    
    parser.add_argument('--min_cov', 
        help="Minimum threshold for each exon's coverage across all samples",
        default=30, type=int)

    parser.add_argument('--cov_readtype',
        help="ctotal = consensus reads, ftotal = full-length reads",
        default='ctotal', type=str)

    parser.add_argument('--ratio_diff', 
        help="Minimum threshold (ratio) for difference in exon usage between 2 or more samples",
        default=0.20, type=float)

    parser.add_argument("-v", "--version", action="version",
                        version=textwrap.dedent("""\
        %(prog)s
        -----------------------   
        Version:    0.2 
        Updated:    2/27/2020
        By:         Prech Uapinyoying   
        Website:    https://github.com/puapinyoying"""))

    # Essential args
    requiredNamed = parser.add_argument_group('Required named arguments:')
    requiredNamed.add_argument('-i', '--input', 
        help="exCOVator output csv file", required=True)
    requiredNamed.add_argument('-o','--output_prefix', 
        help="Prefix for output file", required=True)

    args = parser.parse_args()
    
    return args

if __name__ == '__main__':
    args = parse_input()

    # Run the program
    check_file_exist(args.input)
    df = import_df_from_csv(args.input)
    gene_list, sample_list = get_genes_and_samples(df)
    cov_filter_df, full_filter_df = filter_exons(df, gene_list, sample_list, 
        args.min_cov, args.cov_readtype, args.ratio_diff)

    output_gene_list_from_df(cov_filter_df, '{0}_{1}x_{2}_filtered'.format(
        args.output_prefix, args.min_cov, args.cov_readtype))
    output_gene_list_from_df(full_filter_df, '{0}_{1}x_{2}_{3}_ratio_diff_filtered'.format(
        args.output_prefix, args.min_cov, args.cov_readtype, args.ratio_diff))
    df_to_csv(cov_filter_df, '{0}_{1}x_{2}_filtered_data'.format(
        args.output_prefix, args.min_cov, args.cov_readtype))
    df_to_csv(full_filter_df, '{0}_{1}x_{2}_{3}_ratio_diff_filtered_data'.format(
        args.output_prefix, args.min_cov, args.cov_readtype, args.ratio_diff))