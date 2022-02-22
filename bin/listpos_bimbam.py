#!/usr/bin/env python3
import sys
import argparse

def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--bimbam',type=str,required=False)
    parser.add_argument('--listbimbam',type=str,required=False)
    parser.add_argument('--filepos',type=str,required=False)
    parser.add_argument('--exclude_chr',type=str,required=False,help="File with phenotype and covariate data", default=-1)
    parser.add_argument('--include_chr',type=str,required=False,help="File with phenotype and covariate data")
    parser.add_argument('--out',type=str,required=True,help="File with phenotype and covariate data")

    args = parser.parse_args()
    return args


def extract_listpos(filepos, chro):
  listinfo={}
  lirepossave=open(filepos)
  for line in lirepossave :
   spll=line.split()
   if chro != spll[0] :
     if  spll[0] not in listinfo : 
       listinfo[spll[0]]=set([])
     listinfo[spll[0]].add(spll[1])
  lirepossave.close()
  return listinfo
 
## 
args = parseArguments()
filepos=args.filepos
chro_exclude=args.exclude_chr
chro_include=args.include_chr
filebimbamout=args.out

balisefilepos=False
if filepos :
  listpossave=extract_listpos(filepos,chro_exclude)
  balisefilepos=True


filebimbam=args.bimbam

writebimbam=open(filebimbamout, 'w')
if args.bimbam :
   listbimbam=[args.bimbam]
elif args.listbimbam :
   listbimbam=[x for x in args.listbimbam.split(',') if x[-7::]=='.bimbam']
else :
     print('first column of bim bam file must be rs:chr:pos')
     sys.exit('args..bimbam or args.listbimbam must be initialise')
     sys.exit(2)

for filebimbam in listbimbam:
 readbimbam=open(filebimbam)
 for linebimbam in readbimbam :
   infobimbam=linebimbamspl=linebimbam.split()[0].split(':')
   if(len(infobimbam)!=3) :
     print(':'.join(infobimbam))
     print('first column of bim bam file must be rs:chr:pos')
     sys.exit(2)
   if balisefilepos and (infobimbam[1] in listpossave) and (infobimbam[2] in listpossave[infobimbam[1]]):
     writebimbam.write(linebimbam) 
   elif chro_include and (infobimbam[1] == chro_include) :
     writebimbam.write(linebimbam) 
 readbimbam.close()

writebimbam.close()
