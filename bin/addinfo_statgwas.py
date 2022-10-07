#!/usr/bin/env python3
import os
import argparse
import sys

''' 
generate frequency and N for GxE 
'''


def parseArguments():
    parser = argparse.ArgumentParser(description='merge file statistics with file frequence generate by plink')
    parser.add_argument('--file_stat',type=str,required=True, help="association files")
    parser.add_argument('--file_filelistinfo',type=str,help="pheno names")
    parser.add_argument('--out', type=str,help="",default="out")
    parser.add_argument('--gwas_chr', type=str,help="",default='chr')
    parser.add_argument('--gwas_ps', type=str,help="",default='ps')
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
     if spl[0] not in dicinfo :
        dicinfo[spl[0]]={}
     dicinfo[spl[0]][spl[1]]='\t'+spl[2]+"\t"+spl[3]+"\n"
  liref.close() 

print("--------- end read info file ------")
writestat=open(args.out, 'w')
writestaterror=open(args.out+'.error', 'w')
readstat=open(args.file_stat)
header=readstat.readline().replace('\n','')
spl_head=header.split()
poschr=spl_head.index(args.gwas_chr)
posbp=spl_head.index(args.gwas_ps)
writestat.write(header+"\tTyped\tINFO\n")
for line in readstat :
  linei=line.replace('\n','')
  spl=linei.split()
  try :
     writestat.write(linei+dicinfo[spl[poschr]][spl[posbp]])
  except :
     print(spl[poschr]+"\t"+spl[posbp])
     writestaterror.write(linei+'\tNA\n')
     writestat.write(linei+'\tNA\n')

