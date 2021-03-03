"""
Get all the flanking linear spliced reads for circRNAs to determine circ2linear ratios etc. 
This includes linear spliced reads for zero-count circRNAs which is typically not provided by the circRNA detection pipeline.

usage:
python get_flanking_spliced_reads.py -f *.spliced.bed < circRNA.bed > circRNA.linear_spliced_reads.bed

The spliced.bed files are the output from get_spliced_reads.py

"""

__author__ = "Thomas Hansen (tbh@mbg.au.dk)"
__lab__    = "ncRNAlab"
__version__ = "1.0.0"



import sys, os
import glob
from collections import defaultdict
import argparse
#import mycigar
import pysam

dPos = defaultdict (list)
dJunc = defaultdict (int)
dDist = defaultdict (int)
dCirc = defaultdict (list)


parser = argparse.ArgumentParser(description='Get flanking linear spliced reads for circRNAs')

parser.add_argument("-f","--files", nargs="+", dest="files", default=[], type=str, help="Spliced reads files (*.spliced.bed)")
parser.add_argument("-s","--sep_ss", dest="sep_ss", action = 'store_true', help="Separate sa and sd")


args = parser.parse_args()


header = []

print ("reading stdin...", file=sys.stderr)


# read stdin file
for line in sys.stdin:

    cols = line.strip().split ('\t')

    if cols[0][0] == "#" or not cols[1].isdigit () or not cols[2].isdigit ():
        header = cols
        continue
    
    junction = cols[3]      
    junction = junction[junction.find ('/')+1:]
        
    chrom, start, end, strand = cols[0], int(cols[1]), int(cols[2]), cols[5]
    
    dCirc[(chrom, start, ".")].append (end)
    dCirc[(chrom, end, ".")].append (start)
        
    dPos[(chrom, start, end, strand)] = cols
    


# read spliced.bed files
for i, jfile in enumerate (args.files):
            
    print ("reading", jfile, "-", i+1, "of", len (args.files), file=sys.stderr)

    f = open(jfile, 'r')
    
    for line in f:

        if line[:5] == "track":
            continue
    
        cols = line.strip().split ('\t')
            
        if len(cols) < 6:
            continue
        
        chrom = cols[0]
        start = int (cols[1]) #
        end = int (cols[2])      
        strand = "." #cols[5]
                    

        if len (cols) > 10:            
            overhangs = cols[10].split (',')  # tophat-created junction.bed
        else:
            overhangs = [0,0] 

        intron_start = start + int(overhangs[0])
        intron_end = end - int(overhangs[1])

        
        if (chrom, intron_start, strand) in dCirc:
        
            junction =  dJunc[(chrom, intron_start, strand, jfile)]
            junction = junction + int(cols[4])
            dJunc[(chrom, intron_start, strand, jfile)] = junction

                                
            
        if (chrom, intron_end, strand) in dCirc:
        
            junction =  dJunc[(chrom, intron_end, strand, jfile)]
            junction = junction + int(cols[4])
            dJunc[(chrom, intron_end, strand, jfile)] = junction
                                    


# output

print ("outputting...", file=sys.stderr)


if args.sep_ss:
    output = header + ["total_SA", "total_SD"]
else:
    output = header + ["total_LIN"]
      
for f in args.files:

    colname = os.path.basename(f)
   
    if args.sep_ss:
        output += [colname+"_JUNCS_SA", colname+"_JUNCS_SD"]

    else:
        output += [colname+"_JUNCS_LIN"]
            
print ("\t".join (output))

for (chrom, start, end, strand) in dPos:
               
    output_data = []
    
    total_SA = 0
    total_SD = 0
    total_LIN = 0
    
                     
    cols = dPos[(chrom, start, end, strand)]       


    for jfile in args.files:
    

        SA = dJunc[(chrom, start, ".", jfile)]         
        SD = dJunc[(chrom, end, ".", jfile)]
        
        if strand == "-":
            SA, SD = SD, SA  #NEW!?
    
        total_SA = total_SA + SA
        total_SD = total_SD + SD

        total_LIN = total_LIN + SA + SD

        if args.sep_ss:
            output_data += [str(SA), str(SD)]
        else:   
            output_data += [str(int(SA)+int(SD))]

    if args.sep_ss:
        output = cols + [str(total_SA), str (total_SD)] + output_data
    else:    
        output = cols + [str(total_LIN)] + output_data
  
    print ("\t".join (output))


