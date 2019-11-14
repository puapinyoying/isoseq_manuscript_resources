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

def filter_exons(df, gene_list, sample_list, ratio_diff):
    # Store the filtered gene dataframe to create a new large data frame with results
    filtered_df_list = []
    
    for gene in tqdm(gene_list):
        df1 = df.loc[gene] # grab a slice of the data frame with a single gene
        
        # If there are three samples in csv file
        if len(sample_list) == 3:
            # simplify the sample names into sort variables
            s1 = sample_list[0]
            s2 = sample_list[1]
            s3 = sample_list[2]
            
            # Remove any exons that are fully matched in all 3 samples
            f1 = df1[(df1[s1]['fratio'] < 1) | \
                     (df1[s2]['fratio'] < 1) | \
                     (df1[s3]['fratio'] < 1)]
            
            # Calculate the mean of the exon total coverage of each full-length read column
            # These averages exclude the fully covered exons filtered out above (f1)
            m1 = f1[s1]['ftotal'].mean()
            m2 = f1[s2]['ftotal'].mean()
            m3 = f1[s3]['ftotal'].mean()
            
            # Second filter step is to compare ratios for each exon and see if they differ
            # from each other by a specific threshold (ratio_diff between 0 and 1) 
            # If one of the samples have consistantly low average coverage across all exons 
            if m1 <= 0 and m2 > 0 and m3 > 0:
                # Get the absolute value of the difference & compare to ratio_diff
                f2 = f1[(abs(f1[s2]['fratio'] - f1[s3]['fratio']) >= ratio_diff)]

            elif m2 <= 0 and m1 > 0 and m3 > 0:
                f2 = f1[(abs(f1[s1]['fratio'] - f1[s3]['fratio']) >= ratio_diff)]

            elif m3 <= 0 and m1 > 0 and m2 > 0:
                f2 = f1[(abs(f1[s1]['fratio'] - f1[s2]['fratio']) >= ratio_diff)]

            else: 
                f2 = f1[ (abs(f1[s1]['fratio'] - f1[s2]['fratio']) >= ratio_diff) | \
                         (abs(f1[s1]['fratio'] - f1[s3]['fratio']) >= ratio_diff) | \
                         (abs(f1[s2]['fratio'] - f1[s3]['fratio']) >= ratio_diff) ]
        
        # Its much more straight forward with two samples
        elif len(sample_list) == 2:
            s1 = sample_list[0]
            s2 = sample_list[1]
            f1 = df1[(df1[s1]['fratio'] < 1) | (df1[s2]['fratio'] < 1)]
            
            m1 = f1[s1]['ftotal'].mean()
            m2 = f1[s2]['ftotal'].mean()
            f2 = f1[(abs(f1[s1]['fratio'] - f1[s2]['fratio']) >= ratio_diff)]
                    
        else:
            print "Script does not support > 3 or less than 2 samples"
            sys.exit(1)
            
        if len(f2) > 0: # If its not empty
            ## For lenient filter
            f2_copy = f2.copy() # make a copy of the slice
            f2_copy['gene_name'] = gene # Put the gene name back on the dataframe
            f2_copy.set_index('gene_name', append=True,  inplace=True) # Set it as the index
            
            # Reorder it to original position in front
            filtered_gene_df = f2_copy.reorder_levels(['gene_name', 'exonic_part_num'])
            
            # Sort the df so that unannotated exons will be listed in order by start coordinate
            sorted_filtered_df = filtered_gene_df.sort_values([(sample_list[0], 'start')], ascending=True)
            filtered_df_list.append(sorted_filtered_df)
          
    filtered_df = pd.concat(filtered_df_list)
    
    return filtered_df

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
        
        description="""Script to filter for exons differentially used. Option to
        provide minimum difference between samples as a ratio (e.g. ratio_diff=0.10

        example:""")
    
    parser.add_argument('--ratio_diff', help="", default=0.10, type=float)

    parser.add_argument("-v", "--version", action="version",
                        version=textwrap.dedent("""\
        %(prog)s
        -----------------------   
        Version:    0.1 
        Updated:    11/13/2019
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
    filtered_df = filter_exons(df, gene_list, sample_list, args.ratio_diff)
    output_gene_list_from_df(filtered_df, args.output_prefix)
    df_to_csv(filtered_df, args.output_prefix)
