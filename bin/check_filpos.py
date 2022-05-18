#!/usr/bin/env python3
import sys
''' 
take in input file of rs or position and format as :
  * chr pos pos rsid
  * if 1 columm will open bim file and extract positions 
'''

nbcol=-1
listpos=sys.argv[1]
bim=sys.argv[2]
out=sys.argv[3]
readpos=open(listpos)
writepos=open(out, 'w')
firstline=readpos.readline()
firstline=firstline.replace('\n','').split()
nbcol=len(firstline)
## case of 4 => nothing to do
cmt=1
if nbcol==4 :
  writepos.write("\t".join(firstline)+'\n')
  for line in readpos :
    line=line.replace('\n','').split()
    if nbcol!=len(line):
       print("row "+str(cmt)+" doesn't contain "+str(ncol)+"in file "+listpos+"\nexit\t")
       sys.exit(2)
    writepos.write("\t".join(line)+'\n')
elif nbcol==2 :
  writepos.write(firstline[0]+"\t"+firstline[1]+"\t"+firstline[1]+"\t"+firstline[0]+":"+firstline[1]+'\n')
  for line in readpos :
    line=line.replace('\n','').split()
    if nbcol!=len(line):
       print("row "+str(cmt)+" doesn't contain "+str(ncol)+"in file "+listpos+"\nexit\t")
       sys.exit(2)
    writepos.write(line[0]+"\t"+line[1]+"\t"+line[1]+"\t"+line[0]+":"+line[1]+'\n')
elif nbcol==3 and firstline[2].isnumeric()==False:
  writepos.write(firstline[0]+"\t"+firstline[1]+"\t"+firstline[1]+"\t"+firstline[2]+'\n')
  for line in readpos :
    line=line.replace('\n','').split()
    if nbcol!=len(line):
       print("row "+str(cmt)+" doesn't contain "+str(ncol)+"in file "+listpos+"\nexit\t")
       sys.exit(2)
    writepos.write(line[0]+"\t"+line[1]+"\t"+line[1]+"\t"+line[2]+'\n')
elif nbcol==3 and firstline[2].isnumeric():
  writepos.write(firstline[0]+"\t"+firstline[1]+"\t"+firstline[2]+"\t"+firstline[0]+":"+firstline[1]+":"+firstline[2]+'\n')
  for line in readpos :
    line=line.replace('\n','').split()
    if nbcol!=len(line):
       print("row "+str(cmt)+" doesn't contain "+str(ncol)+"in file "+listpos+"\nexit\t")
       sys.exit(2)
    writepos.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[0]+":"+line[1]+":"+line[2]+'\n')
elif nbcol==1 :
  ## consider as rs read in bimfile
  listrs=set(firstline)
  ## read line
  print('read positions for rs')
  for line in readpos :  
    line=line.replace('\n','').split()
    if nbcol!=len(line):
       print("row "+str(cmt)+" doesn't contain "+str(ncol)+"in file "+listpos+"\nexit\t")
       sys.exit(2)
    listrs.add(line[0])
  print('nb positions to search :',len(listrs))
  print('open bim to search chr and pos')
  lirebim=open(bim)
  for lbim in lirebim :
       line=lbim.replace('\n','').split() 
       if line[1] in listrs :
         writepos.write(line[0]+"\t"+line[3]+"\t"+line[3]+"\t"+line[1]+'\n')
         listrs.discard(line[1])
else :
  print("column number not good : "+nbcol)

writepos.close()
readpos.close()
