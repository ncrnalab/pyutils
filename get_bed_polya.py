#!/usr/bin/env python
import sys
import pysam
import argparse
import numpy as np
import re
from utils import rev_comp

from get_library_type_func import get_library_type_coord



parser = argparse.ArgumentParser(description='Get polyA sequence for each locus')

parser.add_argument('-b', "--bam", dest='bam', type=str, nargs="+", default=[], help="input bamfile")
parser.add_argument('-s', "--stranded", dest='stranded', type=str, default="stranded", help="library type")
parser.add_argument('-g', "--genome", dest='genome', type=str, help="Genome fasta file")

args = parser.parse_args()


fastafile = pysam.FastaFile (args.genome)


for line in sys.stdin:

    cols = line.rstrip('\r\n').split("\t")
        
    if not cols[1].isdigit() or not cols[2].isdigit():
        print ("\t".join (cols + ["polya", "astretch"]))
        continue

    chrom, start, end, strand = cols[0],cols[1], cols[2], cols[5]


    if strand == ".":

        nfw, nre = get_library_type_coord (args.bam, chrom, start, end, "+" if args.stranded == "stranded" else "-")

        strand  = "+" if nfw > nre else "-"
        cols[5] = strand

    polya = 0
    astretch = 0
   
    try:

        # get PAS signal
        seq = fastafile.fetch (chrom, int(start), int(end))


        if strand == "-":
            seq = rev_comp (seq)

        seq = seq.upper()
       
        s = re.search ("A[A|T|U][T|U]AAA", seq)
        if s:
            polya=1

        # get downstream A-richness

        seq = fastafile.fetch (chrom, int(start), int(end) + 50) if strand == "+" else rev_comp (fastafile.fetch (chrom, int(start)-50, int(end)))

        w = 15

        for s in range (0, len (seq)-w):

            wseq = seq[s:s+w]

            nA = wseq.count ("A")

            astretch = max (nA, astretch)

        
    except ValueError:
        polya = 0
        astretch = 0
   


    print ("\t".join (cols + [str (polya), str(astretch)]))

