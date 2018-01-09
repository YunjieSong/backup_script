#!/usr/bin/python
# _*_ coding:UTF-8 _*_

__author__="SongYunjie"

from optparse import OptionParser
parser=OptionParser()
parser.add_option("--inputfile",dest="infile",help="the input file",metavar="File")
parser.add_option("--outputfile",dest="outfile",help="the output file",metavar="File")
(options,args)=parser.parse_args()

infile=options.infile
outfile=options.outfile

dict={}
try:
    with open(infile,"r") as f:
        for line in f:
            (key, value) = line.strip().split(':')
            dict[key]=value
            #print line
except IOError as ioerr:
    print "your inputfile %s is not exist !" %(infile)


try:
    with open(outfile,"w") as o:
        for key in dict:
            o.write("%s:%s\n"%(key,int(dict[key])+10))
except IOError as ioerr:
    print "Can not touch the output file %s"%(outfile)


print dict