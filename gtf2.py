"""
Module for handling gtf files 

"""

__author__ = "Thomas Hansen (tbh@mbg.au.dk)"
__lab__    = "ncRNAlab"
__version__ = "1.0.0"



import sys
from collections import defaultdict
import time
import mmap


rids = ["gene_id", "transcript_id", "gene_name"]


class gtf2_info (object):
    
    def __init__ (self, gtf2, pos, feature):
        self.positions = pos
        self.gtf2 = gtf2
        self.feature = feature
        
        #self.fetch = fetch


    def get_info (self):

        nlines = 0

        with open(self.gtf2.gtf_file) as fgtf:

            for id in self.positions:

                for pos in self.positions[id]:

                    fgtf.seek (int(pos))

                    columns = fgtf.readline().strip().split ("\t")

                    
                    chrom, feature, start, end, strand, info = columns[0], columns[2], columns[3], columns[4], columns[6], columns[8]

                    
                    if feature == self.feature:
                        
                        nlines += 1
                        yield {"id": id, "chrom":chrom, "start":int(start), "end":int(end), "strand":strand, "info":info}
                        

        #if nlines == 0:
        #    yield {"id": "N/A", "chrom": "N/A", "start":-1, "end":-1, "strand":"N/A", "info":"N/A"}
                        

    def get_start (self):
        
        return (min ([info["start"] for info in self.get_info ()]))


    def get_end (self):

        return (max ([info["end"] for info in self.get_info ()]))
        
    def get_coordinates (self):

        start = self.get_start ()
        end   = self.get_end () 
        chrom = ",".join (list (set ([info["chrom"] for info in self.get_info ()])))
        strand = ",".join (list (set ([info["strand"] for info in self.get_info ()])))

        return {"chrom":chrom, "start":start, "end":end, "strand":strand}       

    def get_exons (self):

        exons_s = [info["start"] for info in self.get_info ()]
        exons_e = [info["end"] for info in self.get_info ()]
            
        exons_s.sort()
        exons_e.sort()

        return (exons_s, exons_e)

    def get_introns (self):
        
        exons_s, exons_e = self.get_exons ()

        introns_s = [exons_e[i] for i in range (len (exons_s)-1)]
        introns_e = [exons_s[i+1] for i in range (len (exons_s)-1)]

        return (introns_s, introns_e)

    def get_length (self):
        return len([info for info in self.get_info ()])

    def print (self, feature):

        for info in self.get_info ():
            self.print_info (info, feature)

    
    def print_info (self, info, feature):
            print ("\t".join ([info['chrom'], "gtf2", feature, str (info['start']), str(info['end']), ".", info['strand'], ".", info['info']]))



class gtf2 (object):

    def __init__ (self, gtf_file, head = -1, verbose = 0):

        self.verbose = verbose
        self.nexon_lines = 0
        self.dtid2pos = defaultdict (list)
        self.dgid2tid = defaultdict (list)  
        self.dgid2name = defaultdict (str)  
        self.dtid2gid = defaultdict (str)  
        self.dtid2type = defaultdict (str)  
        self.dchromstart2tid = defaultdict (str)  
        self.dchromend2tid = defaultdict (str)  
        self.dtid2startcodon = defaultdict (int)  
        self.dtid2startcodon_pos = defaultdict (list)  
        
        self.dtid2stopcodon = defaultdict (int)  
        self.dtid2coord = defaultdict (list)  
        
        self.gtf_file = gtf_file

        if self.verbose > 0:
            print (f"Reading gtf {gtf_file}...", file=sys.stderr)
       
        with open(gtf_file, "r+b") as fgtf:

            self.mm = mmap.mmap(fgtf.fileno(), 0, prot=mmap.PROT_READ)
            
            iline = 0

            no_transcript = True

            for line in iter(self.mm.readline, b""):
                
                iline += 1

                # if iline % 10000 == 0:
                #     print (f"...line {iline}", file=sys.stderr)
         
                if iline > head and head != -1:
                    break

                if line == '':
                    break
            
                columns = line.decode().strip().split ('\t')
                

                if len (columns) < 8:
                     continue

                
                chrom, feature, start, end, strand, info =  columns[0], columns[2], int(columns[3]), int(columns[4]), columns[6], columns[8]

                info_proc = self.get_exon_id (info)

                pos = fgtf.tell() - len(line)
                pos = self.mm.tell() - len(line)

                if info_proc["transcript_id"] == "":
                    self.dtid2pos[info_proc["gene_id"]].append (pos) 
                else:
                    self.dtid2pos[info_proc["transcript_id"]].append (pos) 
                
                if feature == "transcript":

                    no_transcript = False
                    
                    self.dgid2tid[info_proc["gene_id"]].append (info_proc["transcript_id"])
                    self.dgid2tid[info_proc["gene_name"]].append (info_proc["transcript_id"])

                    self.dtid2gid[info_proc["transcript_id"]] = info_proc["gene_id"]
                    self.dgid2name[info_proc["gene_id"]] = info_proc["gene_name"]

                    self.dtid2type[info_proc["transcript_id"]] = info_proc["transcript_type"]
                    

                elif feature == "exon":

                    # if transcript feature not in gtf...
                    self.dtid2gid[info_proc["transcript_id"]] = info_proc["gene_id"]
                    self.dgid2name[info_proc["gene_id"]] = info_proc["gene_name"]
                    
                    if no_transcript:

                        if not info_proc["transcript_id"] in self.dgid2tid[info_proc["gene_id"]]:
                            self.dgid2tid[info_proc["gene_id"]].append (info_proc["transcript_id"])
                    
                    self.dchromstart2tid[(chrom,start)] += info_proc["transcript_id"] + ","
                    self.dchromend2tid[(chrom,end)] += info_proc["transcript_id"] + ","

                    self.nexon_lines += 1

                    self.dtid2coord[info_proc["transcript_id"]].append ((chrom, start, end, strand))

                elif feature == "start_codon":
                    self.dtid2startcodon[info_proc["transcript_id"]] = start if strand == "+" else end
                    self.dtid2startcodon_pos[info_proc["transcript_id"]].append ((start, end))
                
                elif feature == "stop_codon":
                    self.dtid2stopcodon[info_proc["transcript_id"]] = start if strand == "+" else end


            if self.verbose > 0:  
                print (f"...found {self.nexon_lines} exons", file=sys.stderr)
            


    def fetch (self, query, feature = "exon"):
        
        q = query
        if isinstance (query, str):
            q = [query]

        positions = {}

        for qi in q:
            positions[qi] = self.dtid2pos[qi]

        # for qi in q: 
        #     if qi in self.dgid2tid:
        #         for tid in self.dgid2tid[qi]:
        #             positions[tid] = self.dtid2pos[tid]
        #     else:        
        #         positions[qi] = self.dtid2pos[qi]

        return gtf2_info (self, positions, feature)
    
   
    def get_exon_id (self, info):
    
        r = {"gene_id": "", \
            "transcript_id":  "", \
            "gene_name":  "",
            "transcript_type": ""}

        for i in info.split ("; "):
            j = i.strip (" ").split (" ")
            if len(j) != 2:
                continue            
            for n in r:
                if n == j[0]:
                    r[n] = j[1].strip("\"")
                    break

        return (r)


    #def get_coordinates (self, tid):    
    #    return (self.fetch (tid, "exon").get_coordinates ())
    
    
    def get_name (self, tid):
        
        gid = self.dtid2gid[tid]
        return (self.dgid2name[gid])

    def get_type (self, tid):
        return (self.dtid2type[tid])

    def get_gene_id (self, tid):        
        return (self.dtid2gid[tid])

    
    def get_startcodon (self, tid):
        return self.dtid2startcodon[tid]


    def get_startcodon_pos (self, tid):
        return (sorted(self.dtid2startcodon_pos[tid],key=lambda x: x[1], reverse=False))
            

    def get_stopcodon (self, tid):        
        return self.dtid2stopcodon[tid]
        
    
    def get_tids_from_gid (self, gid):
        
        if gid in self.dgid2tid:
            return (self.dgid2tid[gid])
        if gid in self.dtid2gid:
            return [gid]
        return []

    def get_tids_from_start (self, chrom, start):
        start = int(start)
        if (chrom, start) in self.dchromstart2tid: 
            return [s for s in self.dchromstart2tid[(chrom, start)].split(",") if s != ""]
        else:
            return []

    def get_tids_from_end (self, chrom, end):
        end = int(end)
        if (chrom, end) in self.dchromend2tid: 
            return [s for s in self.dchromend2tid[(chrom, end)].split(",") if s != ""]
        else:
            return []

    def get_all_tids (self):

        for tid in self.dtid2pos:
            if tid != None:
                yield (tid)

    def get_all_gids (self):

        return [gid for gid in self.dgid2name if gid != None]
        # for gid in self.dgid2name:
        #     if gid != None:
        #         yield (gid)


    def get_exon_coords (self, tid):
        return (sorted(self.dtid2coord[tid],key=lambda x: x[1], reverse=False))


    def get_chrom (self, tid):
        return ("|".join (list(set([c for (c,s,e,t) in self.dtid2coord[tid]]))))


    def get_strand (self, tid):
        return ("|".join (list(set([t for (c,s,e,t) in self.dtid2coord[tid]]))))


    def get_nexons (self, tid):
        return (len (self.get_exon_coords (tid)))


    def get_exons (self, tid):

        sorted_exons = self.get_exon_coords (tid)
        return ([s for (c,s,e,t) in sorted_exons], [e for (c,s,e,t) in sorted_exons])


    def get_introns (self, tid):

        exons_s, exons_e = self.get_exons (tid)
        return (exons_e[:-1], exons_s[1:])

    def get_start (self, tid):
        sorted_exons = self.get_exon_coords (tid)
        return (min ([s for (c,s,e,t) in sorted_exons]))        

    def get_end (self, tid):
        sorted_exons = self.get_exon_coords (tid)
        return (max ([e for (c,s,e,t) in sorted_exons]))        

    def print (self, tid, feature):

        self.fetch (tid, feature).print (feature)

