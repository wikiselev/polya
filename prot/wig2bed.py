#! /usr/bin/env python
# coding: utf-8

# TODO: none
#
# Converting 'wig.formated' into BED format (with cluster shape changing). 
# This script converts 'wig.formated' files into BED files obtained using other 
# scripts. Original data is from CLIP experiments from the paper: 
# 
# Martin, G., Gruber, A. R., Keller, W. & Zavolan, M. 
# Genome-wide Analysis of Pre-mRNA 3� End Processing Reveals a Decisive Role of 
# Human Cleavage Factor I in the Regulation of 3� UTR Length. 
# Cell Rep 1, 753�763 (2012).
#
# Author: Vladimir Kiselev

import sys
import os

if len(sys.argv) == 1:
    print "usage: python wig2bed_new.py wigfile bedfile"
else:
    wigfile = sys.argv[1]
    outfile = sys.argv[2]

f = open(wigfile,'r') # open an input .wig file
f2 = open(outfile,'w') # open an output file
sepstring = 'fixedStep' # "sepstring" is used to separate cluster from .wig file 
ALL = f.read() # read an input file
ALL = ALL.split(sepstring)[1:] # split the input file into clusters by "sepstring"
for all in ALL: # for each cluster
    strand = wigfile.split('_')[1] # define the strand
    desline = all.split(os.linesep, 1)[0]
    chr = desline.split()[0].split('=')[1] # define the chromosome
    clust_start = desline.split()[1].split('=')[1]
    length = all.split(os.linesep, 1)[1]
    length = len(length.split('\n')) - 1
    clust_end = int(clust_start) + length - 1
    print >>f2, '\t'.join([chr, clust_start, str(clust_end), ".", ".", strand])
f.close()
f2.close()
