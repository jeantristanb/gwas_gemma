#!/usr/bin/env python3
import sys

filebimbamread=sys.argv[1]
filebimbamwrite=sys.argv[2]

readdatai=open(filebimbamread)
writedata=open(filebimbamwrite, 'w')

writedata.write(readdatai.readline())

for line in readdatai:
  spl=line.replace('\n','').split()
  info=spl[1].split(':')
  spl[0]=info[1]
  spl[1]=info[0]
  spl[2]=info[2]
  writedata.write("\t".join(spl)+'\n')

writedata.close()
readdatai.close()
