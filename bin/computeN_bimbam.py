#!/usr/bin/env python3


import sys
import argparse

def computed_pheno(data, pheno, cov_list, bimbam_ind):
  readd=open(data)
  header=readd.readline().replace('\n','').split()
  listpospheno=[header.index(pheno)] 
  if len(cov_list) >0:
   listpospheno+=[header.index(x) for x in cov_list.split(',')]
  CmtInd=3
  listindind=[]
  for line in readd:
     spl=line.replace('\n', '').split()
     CmtNa=0
     for pospheno in listpospheno :
        if spl[pospheno].upper()=='NA' or spl[pospheno]=='-9':
          CmtNa+=1 
     if CmtNa==0:
        listindind.append(spl[0]+'_'+spl[1])
     CmtInd+=1
  readd.close()
  Cmt=3
  infoind=[]
  for x in open(bimbam_ind) :
    ind=x.replace('\n', '')
    if ind in listindind :
      infoind.append(Cmt)
    Cmt+=1
  return infoind
def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--bimbam',type=str,required=False)
    parser.add_argument('--bimbam_ind',type=str,required=False)
    parser.add_argument('--data',type=str,required=True,help="File with phenotype and covariate data")
    parser.add_argument('--cov_list', type=str,help="comma separated list of covariates",default="")
    parser.add_argument('--pheno',type=str,required=True,help="comma separated list of  pheno column")
    parser.add_argument('--out', type=str,help="format output : 1:Gemma, 2:boltlmm, 3:FastLmm, 4:gcta 5: gemma bimbam", required=True)
    args = parser.parse_args()
    return args




args = parseArguments()

posindgood =computed_pheno(args.data, args.pheno, args.cov_list, args.bimbam_ind)
readbimbam=open(args.bimbam)
## computed 

#Agincourt/format/bimbam/filt_chr10.vcf_10.bimbam
#.:10:114926 G A

writestat=open(args.out, 'w')
for lbimbam in readbimbam  :
   splbimbam=lbimbam.replace('\n','').split()  
   for posind in posindgood:
       test=float(splbimbam[posind])      
       if test<0 or test>2 :
         print(test)
   chainebimbam=splbimbam[0]+'\t'+splbimbam[1]+'\t'+splbimbam[2]



