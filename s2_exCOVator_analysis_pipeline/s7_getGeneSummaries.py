#! /usr/bin/env python

import os
import sys
import csv
import time
import mygene
import textwrap
import argparse
import itertools
from Bio import Entrez


def check_file_exist(file_path):
    """Make sure the file exists"""
    if not os.path.isfile(file_path):
        file_name = os.path.basename(file_path)
        print "{0} was not found. Check path".format(file_name)
        print "Full path to file: {0}".format(file_path)
        sys.exit(1)


def rm_quotes(word):
    if word.startswith('"') and word.endswith('"') or \
    word.startswith("'") and word.endswith("'"):
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

# The above three funcitons are used often in several of my scripts, see if I 
# can modify the scripts so it can import from one another another 

def get_mygene_ids(gene_symbols, animal):
    """Function that obtains NCBI gene IDs for all of the gene, symbols
    provided. Returns a list of ids found and ids not found. Mygene.info 
    databases do not have a limit to how many inquires you make using their
    servers.  It also automatically batches the requests for you if you ask for
    more than a thousand genes.  This way we can off load the search for gene
    ids here and use the epost function in Entrez to make a single request on
    NCBI. Keeps us safe from the request limit."""
    
    mg = mygene.MyGeneInfo()

    # By choosing all three species, chances are good that you will find a gene
    # id for all of them.
    results = mg.querymany(gene_symbols, scopes='symbol', fields='entrezgene', 
        species=animal) 

    sym_not_found = []
    found_ids = []
    
    # to check if there are multiple entries for each symbol (human, mouse, zeb)
    # If there is, pick the first one.
    prev_symbol = ''
    for info_dict in results:
        if 'notfound' in info_dict:
            sym_not_found.append(info_dict['query'])

        else:
            # First id in series for each symbol is usually human, if human not 
            #mavail, then falls back on mouse , then zebrafish
            # So we pick the first one. Human has the most detailed summaries
            if info_dict['query'] == prev_symbol:
                continue

            try: # returns int; convert to string
                gene_id = str(info_dict['entrezgene']) 
                prev_symbol = info_dict['query']
                found_ids.append(gene_id)

            except KeyError as error_details:
                print "Could not find {0} for {1}".format(error_details, info_dict['query'])
                print info_dict
                sym_not_found.append(info_dict['query'])
                continue

    return found_ids, sym_not_found


# def get_ncbi_ids(gene_list, animal):
#     """This is only here as a backup, NCBI Entrez has the tendency to fall back
#     on orthologs if it can't identify a gene symbol.  This is great, but
#     bombarding NCBI with hundreds or thousands of searches is very bad. It can
#     get you banned.  So this function should only be used for a small set of 
#     gene symbols, for example for ones that were not found using mygene"""

#     gene_count = len(gene_list)
#     id_list = [] # store gene_ids that will be found
    
#     batch_size = 100 # size of queries to do at a time, NCBI has limits
#     for start in range(0, gene_count, batch_size):
#         end = min(count, start+batch_size)
#         print "going to download record {0} to {1}".format(start+1, end)
        
#         for index in range(start, end):
#             gene = gene_list[index]
#             # Search for the gene from the provided gene_name, use search history
#             search_string = gene+"[Gene] AND "+animal+"[Organism]"
#             search_handle = Entrez.esearch(db="gene", term=search_string)
#             search_results = Entrez.read(search_handle)
#             search_handle.close() # close the handle or it will keep using the first gene

#             result_ids = search_results["IdList"]

#             if len(result_ids) > 1:
#                 print "Warning, found {0} gene id's for '{1}'".format(len(result_ids),gene)
#                 print "Going to use the first entry in the list."
#                 gene_id = result_ids[0]
#             elif len(result_ids) < 1:
#                 print "WARNING, could not find gene id for '{0}'".format(gene)
#                 continue
#             else:
#                 gene_id = result_ids[0]

#             id_list.append(gene_id)
        
#     return id_list

def get_gene_info(id_list):
    """With a list of NCBI gene ID numbers (as strings), get summary data for
    each gene.  NCBI returns a complex dictionary of information."""

    try:
        # With the gene_id, grab the gene summaries
        epost_handle = Entrez.epost("gene",id=",".join(id_list))
        epost_record = Entrez.read(epost_handle)

    except RuntimeError as error_report:
        print "Error, trying to grab gene information from NCBI Entrez"
        print error_report

    # Close the first handle
    epost_handle.close()
    
    # Obtain the epost histroy of all the gene_id searches from the first record 
    webEnv = epost_record["WebEnv"]
    queryKey = epost_record["QueryKey"]

    # In case there is a server error or downtime, make 3 attempts before quit
    attempt = 0
    while attempt < 3:
        attempt += 1
        try: # use the history to obtain summaries of all the genes searched
            esummary_handle = Entrez.esummary(db="gene", webenv=webEnv, 
                query_key = queryKey, rettype='summary', retmode="xml")

        except HTTPError as err:
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 3" % attempt)
                time.sleep(15)
            else:
                raise
            
    genes_info_dict = Entrez.read(esummary_handle)
    esummary_handle.close()
    
    return genes_info_dict


def parse_gene_info(genes_info_dict, id_list):
    """Parses the genes info dictionary returned by NCBI esummary for the gene
    symbol, full name and summary. Input the same gene id_list used to generate
    the genes_info_dict. (E.g. found_ids)"""
    count = len(id_list)
    
    summary_list = []
    for index in range(0, count):
        symbol = genes_info_dict['DocumentSummarySet']['DocumentSummary'][index]['NomenclatureSymbol']
        name = genes_info_dict['DocumentSummarySet']['DocumentSummary'][index]['NomenclatureName']
        summary = genes_info_dict['DocumentSummarySet']['DocumentSummary'][index]['Summary']
        row = [symbol, name, summary]
        summary_list.append(row)
        
    return summary_list


def output_summary(summary_list, output_prefix):
    """Outputs the summary of genes into a CSV file. Takes a list of lists"""

    output_name = output_prefix + '_gene_summaries.csv'
    with open(output_name, 'wb') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"')
        csv_writer.writerows(summary_list)
        print "\n1. Gene summaries were output to '{0}'".format(output_name)


def output_not_found(sym_not_found, output_prefix):
    count = len(sym_not_found)
    if count == 0:
        print "2. Gene ids were found for all provided gene symbols."
    else:
        output_name = output_prefix + '_ids_not_found.csv'
        with open(output_name, 'wb') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"')
            for gene in sym_not_found:
                csv_writer.writerow([gene])

        print "2. {0} gene id(s) could not be found in the mygene.info database.".format(count)
        print "Genes that had no corresponding ids wereoutput to '{0}'".format(output_name)
        print "You can try changing the animal/species or adding additional ones."


def parse_input():
    """Use argparse to handle user input for program"""
    
    # Create a parser object
    parser = argparse.ArgumentParser(
        prog='getGeneSummaries.py',
        
        formatter_class=argparse.RawDescriptionHelpFormatter,
        
        description="""Script that takes a csv of gene symbols and returns a csv
        of gene summaries from NCBI.  Useful for quickly screening the gene for
        its function.""")

    #parser.add_argument('-c', '--collapsed', help="", action='store_true')
    parser.add_argument('-a', '--animal', 
        help="""Animal/Species to look for gene ids for. Takes both common names 
        and scientific names. If you want to provide multiple animals, delimit
        using comma (no spaces). E.g. "Human,Mouse". Defaults to Human, Mouse
        and Zebrafish""", default='Human,Mouse,Zebrafish')

    parser.add_argument("-v", "--version", action="version",
                        version=textwrap.dedent("""\
        %(prog)s
        -----------------------   
        Version:    0.1 
        Updated:    01/25/2017
        By:         Prech Uapinyoying   
        Website:    https://github.com/puapinyoying"""))

    # Essential args
    requiredNamed = parser.add_argument_group('Required named arguments:')
    requiredNamed.add_argument('-g', '--gene_file', 
        help="CSV file containing gene symbols formatted in a single column", 
        required=True)
    requiredNamed.add_argument('-o','--output_prefix', 
        help="Specify prefix name for output file", required=True)

    requiredNamed.add_argument('-e','--email', 
        help="Specify email address for NCBI. **Important***", required=True)

    args = parser.parse_args()
    
    return args

#############
#### Run ####
#############

if __name__ == "__main__":
    args = parse_input()

    # Always supply your email to Entrez incase they need to contact you about
    # your usage, or you may get banned for abuses
    Entrez.email = args.email

    # Get list of gene symbols from the provided csv file
    gene_symbols = get_gene_list(args.gene_file)

    # Get lists of gene ids that were found and symbols of ones that were not
    found_ids, sym_not_found = get_mygene_ids(gene_symbols, args.animal)

    # Get the gene summaries from NCBI
    genes_info_dict = get_gene_info(found_ids)

    # Parse the genes_info dict
    summary_list = parse_gene_info(genes_info_dict, found_ids)

    # Output the summaries to a csv file
    output_summary(summary_list, args.output_prefix)
    output_not_found(sym_not_found, args.output_prefix)

    print "Complete!"





