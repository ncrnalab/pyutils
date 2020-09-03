"""
Annotate based on gtf-annotation

usage:
python annotate_circ.py -a <path/to/gtf-file> < circRNAs.bed > circRNAs.annotated.bed

"""

__author__ = "Thomas Hansen (tbh@mbg.au.dk)"
__lab__    = "ncRNAlab"
__version__ = "1.0.0"

import argparse
import gtf2
import time
import sys
import statistics
from collections import defaultdict

parser = argparse.ArgumentParser(description='circRNAs annotation script')
parser.add_argument('-a', "--annot", dest='annot', type=str, help="gtf file")
args = parser.parse_args()

start = time.time()

gtf = gtf2.gtf2 (args.annot)

iline = 0

add_cols = ["name", "gids", "ntids", "tids", 
            "mature_length", "exons", "exon_no", \
            "flank_intron", "intron_up", "intron_down", \
            "host_exons", "host_length", "host_orf_length", "dist5", "dist3", \
            "comments", "circ_type"]


for line in sys.stdin:

    cols = line.rstrip('\r\n').split("\t")

    output = defaultdict (list)
        
    if not cols[1].isdigit () or not cols[2].isdigit ():
        
        print ("\t".join (cols + add_cols))
        continue

    iline += 1

    # convert to 1-based start coordinate
    chrom, start, end, strand = cols[0], int(cols[1])+1, int(cols[2]), cols[5]

    tids_s = gtf.get_tids_from_start (chrom, start)
    tids_e = gtf.get_tids_from_end (chrom, end)
    
    itids = [tid for tid in tids_s if tid in tids_e]

    if len (itids) == 0:
        itids = list(set(tids_s+tids_e))
        
    if len (itids) == 0:
        output['comments'].append ("No annotated SS")
            
    
    output['tids'] = itids
    output['ntids'].append (len (itids))

    name = []

    for ts in itids:

        if gtf.get_strand(ts) != strand:
            continue 
        
        output['gids'].append (gtf.get_gene_id (ts))
        output['name'].append (gtf.get_name (ts))

        # get gene coordinates

        #chrom,  strand = gtf.get_chrom(ts), gtf.get_strand(ts) 
        exon_s, exon_e = gtf.get_exons (ts)

        gene_start = min (exon_s)
        gene_end = max (exon_e)

        #get start and stop codons

        cds_s, cds_e = gtf.get_startcodon(ts), gtf.get_stopcodon(ts)

        nexons = len (exon_s)

        output['host_exons'].append (nexons)
        output['dist5'].append (abs (gene_start - start) if strand == "+" else abs(gene_end - end))
        output['dist3'].append (abs (gene_start - start) if strand == "-" else abs(gene_end - end))
        
        hlength, olength = 0,0

        for i in range (nexons):

            hlength += abs (int (exon_e[i]) - int (exon_s[i]))

            if cds_s != cds_e and int (exon_e[i]) > cds_s and int(exon_s[i]) < cds_e:
    
                olength += min (int (exon_e[i]), cds_e) - max (int (exon_s[i]), cds_s)


        output['host_length'].append (hlength)
        output['host_orf_length'].append (olength)

        r_annot = ""
        
        if start in exon_s or end in exon_e:
            
            if cds_s == cds_e:
                r_annot = "non-coding"
            elif (start < cds_s and end > cds_s and strand == "+") or (start < cds_e and end > cds_e and strand == "-"):
                r_annot = "5utr-cds"
            elif (start < cds_e and end > cds_e and strand == "+") or (start < cds_s and end > cds_s and strand == "-"):
                r_annot = "cds-3utr"
            elif start > cds_s and end < cds_e:
                r_annot = "cds"
            elif (start < cds_s and end < cds_s and strand == "+") or (start > cds_e and end > cds_e and strand == "-"):
                r_annot = "5utr"
            elif (start > cds_e and end > cds_e and strand == "+") or (start < cds_s and end < cds_s and strand == "-"):
                r_annot = "3utr"
        
        output['old_annot'].append (r_annot)
        
        if start in exon_s and end in exon_e:  # Both splice sites found
            
            ies = exon_s.index (start)
            iee = exon_e.index (end)
            
            ei_s_string = ",".join ([str(s) for s in exon_s[ies:(iee+1)]])
            ei_e_string = ",".join ([str(s) for s in exon_e[ies:(iee+1)]])
           
                
            # exon-numbers in circRNA
            eno = "-".join ([str (nexons-s+1) if strand == "-" else str(s) for s in range (ies+1, iee+2)])
            output['exon_no'].append (eno)

            # number of exons in circRNA
            ne = abs (iee - ies) + 1 # How many exons                
            output['exons'].append (str(ne))
            
            # mature_length
            ml = 0
            for i in range (ies, iee+1):
                ml += int (exon_e[i]) - int (exon_s[i])
            
            output['mature_length'].append (str(ml))
            
            # flank intron
            if ies > 0 and iee < len (exon_s) - 1:                
                output['intron_up'].append   ((int(exon_s[ies]) - int(exon_e[ies-1])) if strand == "+" else (int(exon_s[iee+1]) - int(exon_e[iee])))
                output['intron_down'].append ((int(exon_s[ies]) - int(exon_e[ies-1])) if strand == "-" else (int(exon_s[iee+1]) - int(exon_e[iee])))

            elif ies > 0:                
                output['intron_up'].append (int(exon_s[ies]) - int(exon_e[ies-1])) if strand == "+" else output['intron_down'].append (int(exon_s[ies]) - int(exon_e[ies-1]))
            elif iee < len (exon_s) - 1:
                output['intron_up'].append (int(exon_s[iee+1]) - int(exon_e[iee])) if strand == "-" else output['intron_down'].append (int(exon_s[iee+1]) - int(exon_e[iee]))

            

        elif start in exon_s:  # Upstream ss found
            
            ies = exon_s.index (start)
            
            if strand == "+" and ies > 0:
                output['intron_up'].append (int(exon_s[ies]) - int(exon_e[ies-1]))
            elif ies > 0:
                output['intron_down'].append (int(exon_s[ies]) - int(exon_e[ies-1]))
            
            com = "Unannotated SD" if strand == "+" else "Unannotated SA"
            
            output['comments'].append (com)
            
        elif end in exon_e:  # Downstream SS found
        
            iee = exon_e.index (end)
            
            if strand == "+" and iee < len(exon_e)-1:
                output['intron_up'].append (int(exon_s[iee+1]) - int(exon_e[iee]))
            elif iee < len(exon_e)-1:
                output['intron_down'].append (int(exon_s[iee+1]) - int(exon_e[iee]))
                            
            com = "Unannotated SA" if strand == "+" else "Unannotated SD"
            
            output['comments'].append (com)
        

        else:            
            output['comments'].append ("Both splice sites unannotated")
        



    ucols = ['name', 'gids', 'old_annot']

    for u in ucols:
        output[u] = list(set(output[u]))

    new_annot = "N/A"
    if len (output['old_annot']) == 1:
        if output['old_annot'][0] == "5utr-cds":
            new_annot = "AUG circRNA"
        else:
            new_annot = "Other circRNA"
        
    elif len (output['old_annot']) > 1:

        if "5utr-cds" in output['old_annot']:
            new_annot = "Ambiguous AUG circRNA"

        else:
            new_annot = "Ambiguous circRNA"

    output['circ_type'].append (new_annot)
    
    # process outputs - only unique entries:
    
    if len (output['intron_up']) > 0 and len (output['intron_down']) > 0: 
        output['flank_intron'].append (statistics.mean (output['intron_up']) + statistics.mean (output['intron_down']))
        
    for o in add_cols:
        cols += [",".join ([str(s) for s in output[o]])]

    print ("\t".join (cols))
