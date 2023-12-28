#!/usr/bin/env python3
import sys
import os
import gzip
import argparse
def splitgz(x) :
  return x.decode("utf-8").split()

def split(x) :
  return x.split()

def openf(File) :
  print(File)
  def is_gz_file(filepath):
      with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'
  def checkexists(path_to_file) :
    if exists(path_to_file)==False :
     sys.exist('file '+ path_to_file+' doesn t exist')
  balisegz=False
  if is_gz_file(File) :
    readf=gzip.open(File)
    spl=splitgz
    tpw='wt'
  else :
    readf=open(File)
    spl=split
    tpw='w'
  return (readf, spl, tpw)


def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--bimbam',type=str,required=False)
    parser.add_argument('--listbimbam',type=str,required=False)
    parser.add_argument('--filepos',type=str,required=False)
    parser.add_argument('--exclude_chr',type=str,required=False,help="File with phenotype and covariate data")
    parser.add_argument('--include_chr',type=str,required=False,help="File with phenotype and covariate data")
    parser.add_argument('--annotation',type=str,required=True,help="File with phenotype and covariate data")
    parser.add_argument('--out',type=str,required=True,help="File with phenotype and covariate data")

    args = parser.parse_args()
    return args


def extract_listpos_chroexcl(filepos, chro):
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


def extract_listpos_chroincl(filepos, chro):
  listinfo={}
  lirepossave=open(filepos)
  for line in lirepossave :
   spll=line.split()
   if chro == spll[0] :
     if  spll[0] not in listinfo :
       listinfo[spll[0]]=set([])
     listinfo[spll[0]].add(spll[1])
  lirepossave.close()
  return listinfo

def extract_listpos(filepos):
  listinfo={}
  lirepossave=open(filepos)
  for line in lirepossave :
   spll=line.split()
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
  if chro_exclude  :
    listpossave=extract_listpos_chroexcl(filepos,chro_exclude)
  elif chro_include :
    listpossave=extract_listpos_chroincl(filepos,chro_include)
  else :   
    listpossave=extract_listpos(filepos)
  balisefilepos=True




filebimbam=args.bimbam

writeannotation=open(args.annotation, 'w')
if args.bimbam :
   listbimbam=[args.bimbam]
elif args.listbimbam :
   listbimbam=[x for x in args.listbimbam.split(',') if x.lower().endswith(('.bimbam', '.bimbam.gz'))]
else :
     print('first column of bim bam file must be rs:chr:pos')
     sys.exit('args..bimbam or args.listbimbam must be initialise')
     sys.exit(2)

if listbimbam[0].endswith(('.bimbam.gz')) :
  writebimbam=gzip.open(filebimbamout, 'w')
else :
  writebimbam=gzip.open(filebimbamout, 'wt')

CmtSnp=0
for filebimbam in listbimbam:
 (readbimbam,split, tpw)=openf(filebimbam)
 for linebimbam in readbimbam :
   linebimbamspl=split(linebimbam)
   infobimbam=linebimbamspl[0].split(':')
   if(len(infobimbam)!=3) :
     print(':'.join(infobimbam))
     print('first column of bim bam file must be rs:chr:pos')
     sys.exit(2)
   if balisefilepos and (infobimbam[1] in listpossave) and (infobimbam[2] in listpossave[infobimbam[1]]):
     writebimbam.write(linebimbam) 
     writeannotation.write(linebimbamspl[0]+", "+infobimbam[1]+", "+infobimbam[2]+'\n') 
     CmtSnp+=1
   elif chro_include and (infobimbam[1] == chro_include) :
     writeannotation.write(linebimbamspl[0]+", "+infobimbam[1]+", "+infobimbam[2]+'\n') 
     writebimbam.write(linebimbam) 
     CmtSnp+=1
if CmtSnp==0 :
 print('no SNPs used for relatness matrix\nexit\n')
 os._exit(2)

readbimbam.close()
writebimbam.close()
writeannotation.close()
