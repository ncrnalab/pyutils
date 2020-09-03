"""
Get all spliced reads from bam.

usage:
python get_spliced_reads.py -f alignments.bam -g annotation.gtf -fa genome.fa > alignments.spliced.bed

genome.fa must be indexed: samtools faidx genome.dk
annotation.gtf preferable from gencode, may also work using other sources.

"""

__author__ = "Thomas Hansen (tbh@mbg.au.dk)"
__lab__    = "ncRNAlab"
__version__ = "1.0.0"




import sys
import os
import pysam
import argparse
from collections import defaultdict
from utils import rev_comp
import gtf2
import numpy as np



parser = argparse.ArgumentParser(description='Get spliced reads from bam')

parser.add_argument("-f", "--file", dest="file", type=str,help="bamfiles")
parser.add_argument("-g", "--gtf",  dest="gtf", type=str, help="gtf")
parser.add_argument("-fa", "--fasta", dest="fasta", type=str, help="fasta")
parser.add_argument('-s', "--stranded", dest='stranded', type=str, default="reverse", help="Stranded: yes, reverse or no")
parser.add_argument('-t', "--overhang", dest='oh', type=int, default=8, help="Minimun number of aligned reads on each exon") 
parser.add_argument("-mi", "--max_intron", dest='max_intron', type=int, default=500000, help="Max intron length")
parser.add_argument("--head", dest='head', type=int, default=-1, help="Only use the n first reads in bam")


args = parser.parse_args()



lb_type = args.stranded 

print ("library type:", lb_type, file=sys.stderr)

print ("Reading annotation", args.gtf, file=sys.stderr)

gtf = gtf2.gtf2 (args.gtf)

print ("Assemble all known splice sites", file=sys.stderr)

dAnnot = defaultdict (str)

for ts in gtf.get_all_tids ():

    gene_id, chrom, strand = gtf.get_gene_id (ts), gtf.get_chrom (ts), gtf.get_strand (ts)
    intron_s, intron_e = gtf.get_introns (ts)

    for (start, end) in zip (intron_s, intron_e):

        dAnnot[(chrom, start)] = (gene_id, strand)
        dAnnot[(chrom, end)] = (gene_id, strand)


print ("Detecting splicing from bam:", args.file, file=sys.stderr)

dJunc = defaultdict (int)
        
samfile = pysam.Samfile(args.file, "rb")

i = 0                
for read in samfile.fetch():

    if read.mapping_quality  <= 3:
        continue

    i += 1

    if i > args.head and args.head != -1:
        break

    cigar = read.cigarstring        
            
    is_reverse = np.count_nonzero([read.is_reverse, read.is_read2, lb_type == "reverse"]) #ltype == "fr_secondstrand"])

    strand = "-" if is_reverse % 2 == 1 else "+"

    if lb_type == "unstranded":
        strand = "."


    pos = read.pos
    
    match_length_up = 0
    intron_coord = (0,0,0,0)
    
    for ctype, clength in read.cigartuples:
                        
        if ctype == 3: #intron

            intron_coord = (read.reference_name, pos, pos+clength, strand)
            #c.append ((pos, pos+clength)) 
            pos += clength
            

        elif ctype == 0: #match

            if match_length_up != 0 and intron_coord != (0,0,0,0):

                if match_length_up >= args.oh and clength >= args.oh:
                    #print intron_coord
                    dJunc[(intron_coord)] += 1

            match_length_up = clength
        
            pos += clength

        else:
            pos += clength


print ("Outputting confident introns as junctions.bed-like file", file=sys.stderr)

ff = pysam.Fastafile (args.fasta)

ij = 1

print ('track name="' + args.file.split(".")[0] + ".juncs" + '" graphType=junctions description="spliced reads"')

for (chrom, start, end, strand) in dJunc: # intron coordinates

    # check splice sites

    name_array = []

    output_strand = strand  # USE STRAND FROM BAM (not strand from annotation)

    annot_strand_array = []

    if (chrom, start) in dAnnot:
        name, annot_strand = dAnnot[(chrom, start)]
        name_array.append (name)
        annot_strand_array.append (annot_strand)

    if (chrom, end) in dAnnot:
        name, annot_strand = dAnnot[(chrom, end)]
        name_array.append (name)
        annot_strand_array.append (annot_strand)


    annot_strand = "."

    if len(annot_strand_array) >= 1 and list(set(annot_strand_array)) == 1:  #At least one splice sites found (on same strand if two)
        annot_strand = annot_strand_array[0]

    annotated = False

    if len (name_array) >= 2 and len (list(set(name_array))) == 1: # Both sites annotated in same gene_name
        annotated = True


    if not annotated or strand == ".":

        if abs (start-end) > args.max_intron: 
            continue

        sd = ff.fetch (chrom, start, start+2) if strand != "-" else rev_comp (ff.fetch (chrom, end-2, end))
        sa = ff.fetch (chrom, end-2, end) if strand != "-" else rev_comp (ff.fetch (chrom, start, start+2))

        if strand == ".":

            if (sd == "GT" and sa == "AG"):
                #print >>sys.stderr, "Found splice site (+)", sa, sd, chrom, str(start), str(end), str(dJunc[(chrom, start, end, strand)])
                output_strand = "+"

            elif (sd == "CT" and sa == "AC"):
                #print >>sys.stderr, "Found splice site (-)", sa, sd, chrom, str(start), str(end), str(dJunc[(chrom, start, end, strand)])
                output_strand = "-"

            else:
                #print >>sys.stderr, "unstranded, no splice sites", sa, sd, chrom, str(start), str(end), str(dJunc[(chrom, start, end, strand)])
                continue

        elif sd != "GT" or sa != "AG":
            #print >>sys.stderr, "stranded, no splice sites", sa, sd, chrom, str(start), str(end), str(dJunc[(chrom, start, end, strand)])
            continue


    if output_strand != "." and annot_strand != "." and output_strand != annot_strand:
        print (chrom, start, end, dJunc[(chrom, start, end, strand)], "- STRAND DISCREPANCY... found:", output_strand, ", annotated as_", annot_strand, ". Using annotated strand", file=sys.stderr)
    
    name_tag = "|".join (list(set(name_array))) if len(name_array) > 0 else "junc_" + str(ij)

    print ("\t".join ([chrom, str(start), str(end), name_tag, str(dJunc[(chrom, start, end, strand)]), output_strand, str(start), str(end), "100,100,100"]))

    ij += 1





    