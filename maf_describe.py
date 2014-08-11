#!/usr/bin/env python3
"""
Read a MAF file, and describe it. 
"""
import pandas as pd
import argparse
import sys

MAF_COLUMN_NAMES = [
    'Hugo_Symbol',
    'Entrez_Gene_Id',
    'Center',
    'NCBI_Build',
    'Chromosome',
    'Start_Position',
    'End_Position',
    'Strand',
    'Variant_Classification',
    'Variant_Type',
    'Reference_Allele',
    'Tumor_Seq_Allele1',
    'Tumor_Seq_Allele2',
    'dbSNP_RS',
    'dbSNP_Val_Status',
    'Tumor_Sample_Barcode',
    'Matched_Norm_Sample_Barcode',
    'Match_Norm_Seq_Allele1',
    'Match_Norm_Seq_Allele2'
]

def run(args):
    maf_df = get_maf_df(args)
    allele_ref = maf_df['Reference_Allele']
    allele_1 = maf_df['Tumor_Seq_Allele1']
    allele_2 = maf_df['Tumor_Seq_Allele2']
    allele_1_vars = maf_df[(allele_ref == allele_2) & (allele_ref != allele_1)]
    allele_2_vars = maf_df[(allele_ref == allele_1) & (allele_ref != allele_2)]
    allele_both_vars = maf_df[(allele_ref != allele_1) & (allele_ref != allele_2)]
    allele_none_vars = maf_df[(allele_ref == allele_1) & (allele_ref == allele_2)]
    print("Number of lines: %d" % len(maf_df))
    print("Allele 1 is different from the ref and 2 is the same: %d" % len(allele_1_vars))
    print("Allele 2 is different from the ref and 2 is the same: %d" % len(allele_2_vars))
    print("Both alleles are different from the ref: %d" % len(allele_both_vars))
    print("Both alleles are the same as the ref: %d" % len(allele_none_vars))

def get_maf_df(args):
    filename = args.filename[0]
    df = pd.read_csv(filename, sep='\t')
    with open(filename) as fd:
        lines_to_skip = 0
        while next(fd).startswith('#'):
            lines_to_skip += 1
        return pd.read_csv(
                filename,
                skiprows=lines_to_skip,
                sep='\t',
                usecols = range(len(MAF_COLUMN_NAMES)),
                names=MAF_COLUMN_NAMES)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("filename",
            metavar="filename", nargs=1,
            help="a MAF file")
    args = parser.parse_args()
    run(args)
