"""
Get alternative splicing events. For each splice-site (SD or SA), all paired splice sites are collected. 
The most abundant is classified as typical, whereas the remaining pairs are either 'short' or 'long' relative to the canonical. 

usage:
python get_alternative_splicing.py -f *.spliced.bed > alternative_splicing.bed

The spliced.bed files are the output from get_spliced_reads.py

"""

__author__ = "Thomas Hansen (tbh@mbg.au.dk)"
__lab__    = "ncRNAlab"
__version__ = "1.0.0"



import sys, os
import glob
from collections import defaultdict
import argparse
import pysam
import numpy as np
import gtf2

dSS = defaultdict (list)
aSS = defaultdict (int)

dAnnot = defaultdict (str)

parser = argparse.ArgumentParser(description='Get alternative splicing.')
parser.add_argument("-f","--files", nargs="+", dest="files", default=[], type=str, help="junction files (*.spliced.bed)")
args = parser.parse_args()


files = []

for x in args.files:
   files.append (x)

files.sort()

   
for filename in files:

    print ("Reading...", filename, file=sys.stderr)
    
    f = open(filename, 'r')
   
    for line in f:
     
        cols = line.split ('\t')

        if len (cols) < 6:
            continue

        if not cols[1].isdigit() or not cols[2].isdigit():            
            continue

        chrom, start, end, name, strand = cols[0], cols[1], cols[2], cols[3], cols[5]
        nreads = int(cols[4])

        (sd, sa) = (start, end) if strand == "+" else (end, start)

        dSS[chrom, sd, strand, filename, "SD"].append ({"ss":sa, "score":nreads})
        dSS[chrom, sa, strand, filename, "SA"].append ({"ss":sd, "score":nreads})
        
        aSS[chrom, sa, strand, "SA"] = 1
        aSS[chrom, sd, strand, "SD"] = 1
        


output =  ["chrom", "start", "end", "name", "score", "strand", "ss_type", "alternative"] + files

print ("\t".join (output))


for (chrom, ss_pos, strand, ss_type) in aSS:
                
    dSUM = defaultdict (int)
    dF = defaultdict (int)

    for f in files:
    
        for ss_pair in dSS[(chrom, ss_pos, strand, f, ss_type)]:
            
            dSUM[ss_pair['ss']] += ss_pair['score']
            dF[ss_pair['ss'], f] += ss_pair['score']

    
    sorted_tuple = sorted(dSUM.items(), key=lambda x: x[1], reverse=True)[0]

    sorted_ss, sorted_count = sorted_tuple
    
    typical_length = abs (int(sorted_ss) - int(ss_pos))

    for sds in dSUM:
    
        sum_counts = dSUM[sds]
    
        ss_length = abs (int(sds) - int(ss_pos))

        rel = ss_length - typical_length

        alternative = "Typical"
        if ss_length == typical_length:
            alternative = "Typical"
        elif ss_length < typical_length:
            alternative = "Short"
        else:
            alternative = "Long"

        output = [chrom, str(min (int(sds), int(ss_pos))), str(max(int(sds), int(ss_pos))), ".", str(sum_counts), strand, ss_type, str(rel)]

        for f in files:

            counts = dF[sds, f]

            output = output + [str (counts)]
        
        print ("\t".join (output))



