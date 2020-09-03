"""
Annotate inverted Alu elements (or other repeat elements) based on repeatmask (sqlite db)

Convert mySQL to sqlite: https://gist.github.com/esperlu/943776. For UCSC: 

usage:
python annotate_repeat.py -rdb <path/to/repeatmask.sqlite> < circRNA.bed > circRNA.IAE.bed


"""

__author__ = "Thomas Hansen (tbh@mbg.au.dk)"
__lab__    = "ncRNAlab"
__version__ = "1.0.0"


import sys
import argparse
#import utils
import sqlutil


parser = argparse.ArgumentParser(description='Annotate circRNAs')

parser.add_argument('-rdb', "--rdb", dest='rdb', help="sqlite database with repeat annotation")
parser.add_argument('-rf', "--repeat_family", dest='repeat_family', nargs="+", default=["Alu"], help="Repeat families to analyse")  # ["MIR", "L1", "L2"]

args = parser.parse_args()

sqlite_repeat_db = args.rdb
repeat_family = args.repeat_family 

sql_rep = sqlutil.sqlutil (type="sqlite3", url=sqlite_repeat_db)
repeat_db = sql_rep.get_table ()

iline = 0

for line in sys.stdin:

    iline += 1
            
    cols = line.rstrip('\r\n').split ('\t')

    if len (cols) < 3:
        continue
    
    # header
    if cols[0][0] == "#" or not cols[1].isdigit () or not cols[2].isdigit:
        
        for family in repeat_family:
            
            f = family.lower()
            cols = cols + [f+"dist", f+"dist_up_down", f+"name_up_down", f+"score_up_down"]
    
        print ("\t".join (cols))
        continue

    chrom, start, end, strand = cols[0], int(cols[1]), int(cols[2]), cols[5]
    
    for family in repeat_family:
    
        #downstream
        
        query = "SELECT *, start-{} AS rel FROM {} WHERE family='{}' AND chr='{}' AND start > '{}' ORDER BY start ASC LIMIT 20".format (
            end, repeat_db, family, chrom, end)
        
        sql_rep.query(query)  
                            
        inf = 99999999
        
        dist_down_p, dist_down_m, dist_up_p, dist_up_m = inf, inf, inf, inf
        name_down_p, name_down_m, name_up_p, name_up_m = "", "", "", ""
        score_down_p, score_down_m, score_up_p, score_up_m = 0,0,0,0
        
        for rr in sql_rep.fetchall():
        
            if rr["strand"] == "+" and int(rr["rel"]) < dist_down_p:
                dist_down_p = int(rr["rel"])
                name_down_p = rr["name"]
                score_down_p = rr["score"]
                
            elif rr["strand"] == "-" and int(rr["rel"]) < dist_down_m:
                dist_down_m = int(rr["rel"])
                name_down_m = rr["name"]
                score_down_m = rr["score"]
                
            if dist_down_p < inf and dist_down_m < inf: 
                break

        query = "SELECT *, {}-end AS rel FROM {} WHERE family='{}' AND chr='{}' AND end < '{}' ORDER BY end DESC LIMIT 20".format (
            start, repeat_db, family, chrom, start)
        
        sql_rep.execute(query)  
                    
        for rr in sql_rep.fetchall():

            if rr["strand"] == "+" and int(rr["rel"]) < dist_up_p:                                  
                dist_up_p = int(rr["rel"])
                name_up_p = rr["name"]
                score_up_p = int(rr["score"])
                
            elif rr["strand"] == "-" and int(rr["rel"]) < dist_up_m:
            
                dist_up_m = int(rr["rel"])
                name_up_m = rr["name"]
                score_up_m = (rr["score"])
                
            
            if dist_up_p < inf and dist_up_m < inf: 
                break
        
        prox_up = min (dist_up_p, dist_up_m)             
        prox_down = min (dist_down_p, dist_down_m)
        
        dist = min (dist_up_p+dist_down_m, dist_up_m+dist_down_p)
        dist_min = min (dist_up_p, dist_up_m) + min (dist_down_p, dist_down_m)
                    
        if dist > inf:
            dist = ""
        
        (dist_up, dist_down) = (dist_up_p,dist_down_m) if (dist_up_p+dist_down_m < dist_up_m+dist_down_p) else (dist_up_m,dist_down_p)
        (name_up, name_down) = (name_up_p,name_down_m) if (dist_up_p+dist_down_m < dist_up_m+dist_down_p) else (name_up_m,name_down_p)               
        (score_up, score_down) = (score_up_p,score_down_m) if (dist_up_p+dist_down_m < dist_up_m+dist_down_p) else (score_up_m,score_down_p)
        
        if strand == "-":
            dist_up, dist_down = dist_down, dist_up
            name_up, name_down = name_down, name_up
            score_up, score_down = score_down, score_up
            
        cols.append (str(dist))
        cols.append (",".join ([str(s) for s in [dist_up, dist_down]]))
        cols.append (",".join ([str(s) for s in [name_up, name_down]]))
        cols.append (",".join ([str(s) for s in [score_up, score_down]]))
        
    print ("\t".join (cols))

sql_rep.close()

      
