#Get transcripts with exons matching start / end
#then annotate based on 

import argparse
import gtf2
import time
import sys
import statistics
from collections import defaultdict

parser = argparse.ArgumentParser(description='Get GTF info')

parser.add_argument('-a', "--annot", dest='annot', type=str, help="gtf file")

args = parser.parse_args()

start = time.time()

gtf = gtf2.gtf2 (args.annot)

iline = 0

add_cols = ["name", "gids", "ntids", "tids", 
            "exons", "total_exon_length", "total_intron_length", "longest_exon", "longest_intron", "ORF_length", "ORF_exon_start"]


for line in sys.stdin:

    if len (line) > 0 and line[0] == "#":
        print (line.strip() + " with annotation by annotate_fc.py using gtf:", args.annot)
        continue
    
    cols = line.strip().split("\t")
    
    if cols[0] == "Geneid": # featureCounts

        print ("\t".join (cols + add_cols))
        continue

    output = defaultdict (list)


    iline += 1

    gid = cols[0]
    
    itids = gtf.get_tids_from_gid (gid)

    if len (itids) == 0:
        print (f"gene_id not found: {gid}", file=sys.stderr)
        continue


    output['tids'] = itids
    output['ntids'].append (len (itids))

    name = []

    gene_start, gene_end = [], []
    chrom, strand = cols[1], cols[4]


    for ts in itids:

        output['gids'].append (gtf.get_gene_id (ts))
        output['name'].append (gtf.get_name (ts))
        output['type'].append (gtf.get_type (ts))

        # get gene coordinates

        chrom,  strand = gtf.get_chrom(ts), gtf.get_strand(ts) 
        exon_s, exon_e = gtf.get_exons (ts)

        gene_start.append (min (exon_s))
        gene_end.append (max (exon_e))

        #exoninfo = gtf.fetch (ts, "exon")
        #coord = exoninfo.get_coordinates ()
        #r_chrom, r_start, r_end, r_strand = coord['chrom'], coord['start'], coord['end'], coord['strand'] #info.get_start (), info.get_end ()
        
        #min_start = r_start if min_start == 0 else min (min_start, r_start)
        #max_end = r_end if max_end == 0 else max (max_end, r_end)
        #strand.append (r_strand)
        #chrom.append (r_chrom)

        # fc overwrite

        #get start and stop codons

        cds_s, cds_e = gtf.get_startcodon(ts), gtf.get_stopcodon(ts)

        #exon_s, exon_e = exoninfo.get_exons ()

        nexons = len (exon_s)

        output['exons'].append (nexons)

        aug_exon = -1
        sum_orf, sum_intron, sum_exon, max_intron, max_exon = 0,0,0,0,0

        utr3_start, utr3_end = [],[]

        for i in range (nexons):

            exon_length = abs (int (exon_e[i]) - int (exon_s[i]))
            sum_exon += exon_length
            max_exon = max (max_exon, exon_length)

            if (i >= 1):
                intron_length = abs (int (exon_s[i]) - int (exon_e[i-1]))
                sum_intron += intron_length
                max_intron = max (max_intron, intron_length)


            if cds_s != cds_e and int (exon_e[i]) > cds_s and int(exon_s[i]) < cds_e:
    
                sum_orf += min (int (exon_e[i]), cds_e) - max (int (exon_s[i]), cds_s)

            aug = cds_s if strand == "+" else cds_e 

            if aug >= exon_s[i] and aug <= exon_e[i] and cds_s != cds_e:
                aug_exon = i+1


            utr_up_s, utr_up_e, utr_down_s, utr_down_e = -1,-1,-1,-1
            # utr_up 
            if int (exon_s[i]) < cds_s:
                utr_up_s = int(exon_s[i])
                utr_up_e = min (int(exon_e[i]), cds_s)

            # utr_up 
            if int (exon_e[i]) > cds_e:
                utr_down_s = max (int(exon_s[i]), cds_e)
                utr_down_e = int(exon_e[i])

            utr3_start.append (utr_down_s if strand == "+" else utr_up_s)
            utr3_end.append (utr_down_e if strand == "+" else utr_down_e)
            



        if strand == "-":
            aug_exon = nexons - aug_exon + 1
                 

        output['total_exon_length'].append (sum_exon)
        output['longest_exon'].append (max_exon)
        
        output['total_intron_length'].append (sum_intron)
        output['longest_intron'].append (max_intron)
                
        output['ORF_length'].append (sum_orf)
        output['ORF_exon_start'].append (aug_exon)


    ucols = ['name', 'gids']

    for u in ucols:
        output[u] = set(list(output[u]))

    for o in add_cols:
        cols += [",".join ([str(s) for s in output[o]])]

    # overwrite fc columns 1,2,3,4

    cols[1] = chrom #",".join (list(set(chrom)))
    cols[2] = str(min(gene_start)) #str(min_start)
    cols[3] = str(max(gene_end)) #str(max_end)
    cols[4] = strand #",".join (list(set(strand)))
    

    print ("\t".join (cols))



