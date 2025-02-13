#!/usr/bin/env python3
import os
import argparse
import sys

''' 
generate frequency and N for GxE 
'''

#chr	rs	ps	n_miss	allele1	allele0	af	beta	se	logl_H1	l_remle	l_mle	p_wald	p_lrt	p_score	SNP	N
def parseArguments():
    parser = argparse.ArgumentParser(description='merge file statistics with file frequence generate by plink')
    parser.add_argument('--file_stat',type=str,required=True, help="association files")
    parser.add_argument('--file_filelistinfo',type=str,help="pheno names")
    parser.add_argument('--out', type=str,help="",default="out")
    #parser.add_argument('--gwas_chr', type=str,help="",default='chr')
    #parser.add_argument('--gwas_ps', type=str,help="",default='ps')
    args = parser.parse_args()
    return args

args = parseArguments()

## read 
listfile=[x.replace('\n','') for x in open(args.file_filelistinfo).readlines()]
dicinfo={}
print("--------- begin read info file ------")
for File in listfile :
  liref=open(File) 
  print("file "+File+" analyse") 
  head=liref.readline()
  for line in liref:
     spl=line.replace('\n','').split()
     chro=spl[0]
     chro=chro.replace('chr','').replace('X', '23')
     if chro not in dicinfo :
        dicinfo[chro]={}
     dicinfo[chro][spl[1]]='\t'+spl[2]+"\t"+spl[3]+"\n"
  liref.close() 

print("--------- end read info file ------")
writestat=open(args.out, 'w')
writestaterror=open(args.out+'.error', 'w')
readstat=open(args.file_stat)
header=readstat.readline().replace('\n','')
spl_head=header.split()
#writestat.write(header+"\tTyped\tINFO\n")
writestat.write("\t".join(['VARIANTID','RSID','CHR','BP','EFFECT_ALLELE','OTHER_ALLELE','N','EAF','BETA','SE','P','R2','INFO'])+'\n')
for line in readstat :
  linei=line.replace('\n','')
  (chro,rs,ps,n_miss,allele1,allele0,af,beta,se,logl_H1,l_remle,l_mle,p_wald,p_lrt,p_score,SNP,N)=linei.split()
  chro=chro.replace('chr','')
  id1=chro+"_"+ps+"_"+allele1+"_"+allele0
  newl="\t".join([id1,SNP,chro,ps,allele1,allele0,N,af,beta,se,p_wald]) 
  try :
     writestat.write(newl+dicinfo[chro][ps])
  except :
     writestaterror.write(newl+'\tNA\tNA\n')
     writestat.write(newl+'\tNA\tNA\n')

