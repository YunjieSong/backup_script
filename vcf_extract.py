#!/usr/bin/python
#_*_ coding:UTF-8 _*_

from optparse import OptionParser
import collections

__author__="SongYunjie"

usage="usage: python %prog -i arg1"
parser=OptionParser(usage)
parser.add_option("-i","--input",dest="inputfile",help="input vcf file")
(options,args)=parser.parse_args()

infile=options.inputfile
# print infile
result=[]
with open(infile, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        else:
            infor=line.strip().split("\t")[7].split(";")
            for i,s in enumerate(infor):
                if "culprit" in s:
                    result.append(infor[i].split("=")[1])
print result
print len(result)
print collections.Counter(result)