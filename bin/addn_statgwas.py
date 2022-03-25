#!/usr/bin/env python3
import os
import math
import argparse
import numpy as np
import pandas as pd
import sys

''' 
generate frequency and N for GxE 
'''


def parseArguments():
    parser = argparse.ArgumentParser(description='merge file statistics with file frequence generate by plink')
    parser.add_argument('--file_stat',type=str,required=True, help="association files")
    parser.add_argument('--file_freq',type=str,help="pheno names")
    parser.add_argument('--out', type=str,help="",default="out")
    parser.add_argument('--gwas_chr', type=str,help="",default='chr')
    parser.add_argument('--gwas_ps', type=str,help="",default='ps')
    parser.add_argument('--gwas_rs', type=str,help="",default='rs')
    args = parser.parse_args()
    return args

args = parseArguments()

datastat= pd.read_csv(args.file_stat,delim_whitespace=True)
freqAll = pd.read_csv(args.file_freq,delim_whitespace=True)
freqAll = freqAll[['CHR','SNP', 'NCHROBS']]
freqAll['NCHROBS']=freqAll['NCHROBS']/2
freqAll['CHR']=freqAll['CHR'].replace(' ','').astype(str)
freqAll['SNP']=freqAll['SNP'].replace(' ','')
freqAll=freqAll.rename(index=str, columns={"MAF": "Freq_All", "NCHROBS": "N", "CHR":args.gwas_chr, "SNP":args.gwas_rs})
datastat[args.gwas_chr]=datastat[args.gwas_chr].astype(str)
datastatAll=pd.merge(datastat,freqAll,how="left",on=[args.gwas_chr, args.gwas_rs])
datastatAll.to_csv(args.out, sep='\t', na_rep='NA', header=True, index=False, mode='w')


