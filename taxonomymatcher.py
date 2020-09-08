#!/usr/bin/python

import pandas as pd 
import sys, getopt
import argparse

'''
Usage: taxonomymatcher.py -i <kallisto_abundance_file.tsv> -o <output.tsv>
taxonomymatcher.py compiles Kallisto output with PhyloDB taxonomy and annotation databases.
'''

def parse_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', metavar="INPUT_FILE",
                        help='Path to kallisto/salmon pseudocount input file.',
                        required=True, dest='input_file', type=str)
    parser.add_argument('-o', '--outfile', metavar="OUTPUT_FILE",
                        help='Path to write merged taxonomy/annotation file',
                        required=True, dest='out_file', type=str)
    parser.add_argument('-a', '--annotations', metavar="ANNOTATIONS_FILE",
                        help='Path to phylodb annotations text file.',
                        required=True, dest='annot_file', type=str)  
    parser.add_argument('-t', '--taxonomies', metavar="TAXONOMY_FILE",
                        help='Path to input phylodb taxonomy file',
                        required=True, dest='tax_file', type=str)   
    args = parser.parse_args()

    return args

def get_databases(kallisto_abundance_file, annotations, taxonomies):
    '''
    Loads files from databases and input file, preps for merging
    '''
    abundances = pd.read_csv(kallisto_abundance_file, sep='\t')

    gene_database = pd.read_csv(annotations, sep='\t', \
        names=["gene_id", "idstring", "strain_name", "gene_name"])
    

    taxonomy_database = pd.read_csv(taxonomies, sep='\t')


    return abundances, gene_database, taxonomy_database

def data_merger(abundances, gene_database, taxonomy_database):
    '''
    Merges data using joins, cleans up unnecessary columns and 0 TPM rows
    '''

    abundances = abundances.merge(gene_database, left_on = abundances["target_id"], right_on = gene_database["gene_id"])
    
    abundances = abundances.merge(taxonomy_database, left_on = abundances["strain_name"], right_on = taxonomy_database["strain_name"])
    
    abundances = abundances.drop(['strain_name_y', 'idstring', 'gene_id', 'peptide_count'], axis=1)

    #abundances = abundances.drop(abundances[abundances['tpm'] == 0].index) # Can edit in/out depending on if you want to filter

    return abundances

def output_writer(merged_dataframe, outputfile):
    merged_dataframe.to_csv(path_or_buf=outputfile, index=False)

def main():
    
    user_args = parse_arguments()
    
    inputfile = user_args.input_file
    outputfile = user_args.out_file
    annot_file = user_args.annot_file
    tax_file = user_args.tax_file
        

    abundances, gene_database, taxonomies = get_databases(inputfile, annot_file, tax_file)

    merged_table = data_merger(abundances, gene_database, taxonomies)
    output_writer(merged_table, outputfile)

if __name__ == "__main__":
    main()
