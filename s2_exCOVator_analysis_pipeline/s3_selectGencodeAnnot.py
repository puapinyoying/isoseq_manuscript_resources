import os
import sys
import csv
import HTSeq
import itertools
import argparse
import textwrap
from tqdm import tqdm

def check_file_exist(file_path):
    """Make sure the file exists"""
    if not os.path.isfile(file_path):
        file_name = os.path.basename(file_path)
        print "{0} was not found. Check path".format(file_name)
        print "Full path to file: {0}".format(file_path)
        sys.exit(1)

        
def rm_quotes(word):
    if word.startswith('"') and word.endswith('"'):
        word = word[1:-1]
    return word
    

def get_gene_list(gene_file):
    check_file_exist(gene_file)
    gene_list = []

    with open(gene_file, 'rb') as csvfile:
        gene_reader = csv.reader(csvfile, delimiter='\t')
        for row in gene_reader:
            if len(row) > 1:
                print "The gene file: {0} has more than one column."
                sys.exit(1)
            else:
                gene_name = rm_quotes(row[0])
                gene_list.append(gene_name)

    return gene_list


def select_annot_by_gene(input_gtf, output_gtf, gene_list):
    check_file_exist(input_gtf)
    with open(input_gtf, 'rb') as csv_file:
        # quotechar='|' prevents double quoting and writes out same way
        gtf_reader = csv.reader(csv_file, delimiter='\t', quotechar='|')

        with open(output_gtf, 'wb') as gtf_file:
            gtf_writer = csv.writer(gtf_file, delimiter='\t', quotechar='|') 

            for row in tqdm(gtf_reader, "Sifting through {0}".format(input_gtf)):
                if row[0].startswith('#'):
                    continue # skip over header rows

                # column 8 has gtf attributes delimited by ';'
                attr_list = row[8].split(';')
                for item in attr_list:
                    if item: # skip empty columns
                        pair = item.split() # split each item by space
                        if pair[0] == 'gene_name':
                            # names are double quoted and must be removed
                            name = rm_quotes(pair[1]) 
                            if name in gene_list:
                                gtf_writer.writerow(row)

    print "Complete! Selected annotations output to {0}".format(output_gtf)


def parse_input():
    """Use argparse to handle user input for program"""
    # Create a parser object
    parser = argparse.ArgumentParser(
        prog='selectGencodeAnnot.py',
        
        formatter_class=argparse.RawDescriptionHelpFormatter,
        
        description="""Script to select specified gene annotations from a 
        Gencode (GTF) file and place them into a new gtf file.""")

    # Parse the polyA or the noPolyA groups
    parser.add_argument('-c', '--collapsed', help="", action='store_true')
    parser.add_argument("-v", "--version", action="version",
                        version=textwrap.dedent("""\
        %(prog)s
        -----------------------   
        Version:    0.1 
        Updated:    01/13/2017
        By:         Prech Uapinyoying   
        Website:    https://github.com/puapinyoying"""))

    # Essential args
    requiredNamed = parser.add_argument_group('Required named arguments:')
    requiredNamed.add_argument('-i', '--input_gtf', 
        help="Gencode annotation GTF file", required=True)
    requiredNamed.add_argument('-g','--gene_file', 
        help="List of gene names listed in a single column", required=True)
    requiredNamed.add_argument('-o','--output_gtf', 
        help="The selected gene annotation output gtf", required=True)

    args = parser.parse_args()
    
    return args

#############
#### Run ####
#############

if __name__ == "__main__":
    args = parse_input()

    gene_list = get_gene_list(args.gene_file)

    select_annot_by_gene(args.input_gtf, args.output_gtf, gene_list)

