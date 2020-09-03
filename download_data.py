"""
Download RNAseq data from GEO or ENA - NB: still under development!

usage:
python download_data.py -q <accession_no> -d

The resulting fastq.gz files will be downloaded to location ./<accession_no>/fastq 

If found on ENA, md5sum data will be downloaded, checked, and outputted in ./<accession_no>/fastq/md5sum.txt.report.


"""

__author__ = "Thomas Hansen (tbh@mbg.au.dk)"
__lab__    = "ncRNAlab"
__version__ = "1.0.0"


import argparse
import sys
import os
import ssl
from bs4 import BeautifulSoup as BS  # Also needs lxml: pip install --cert /com/etc/ssl-proxy-cert.pem lxml
import hashlib
import requests


parser = argparse.ArgumentParser(description='Download data from GEO/ENA')

parser.add_argument('-q', "--query", dest='query', type=str, default="",help='Query - typically GSE accession')
parser.add_argument('-g', "--grep", dest='grep', type=str, default="", help='Only download matching dataset')
parser.add_argument('-vg', "--vgrep", dest='vgrep', type=str, default="", help='select non-matching lines')
parser.add_argument('-o', "--output", dest='output', type=str, default="", help='output for download-commands')
parser.add_argument('-i', "--info", dest='info', type=str, default="", help='output for sample info')
parser.add_argument("--verbose", dest='verbose', type=int, default=0, help='verbose')
parser.add_argument('-db', "--database", dest='db', type=str, default="", help='database')
#parser.add_argument('-md5', dest='md5', type=str, default="md5sum.txt", help='md5 output file')

parser.add_argument("-d", "--download", dest="download", help="Download data from GEO. If this option is not specified, only show metadata.", action="store_true", default=False)  

args = parser.parse_args()



# INCLUDE ENA search: https://www.ebi.ac.uk/ena/browse/search-rest
# INCLUDE ARRAY EXPRESS: https://www.ebi.ac.uk/arrayexpress/help/programmatic_access.html
class gid:

    def __init__(self, db, id, alias = ""):
        self.db = db
        self.gid = id
        self.alias = alias
        

class sid ():

    def __init__(self, gid, acc, alias, desc=""):
        self.gid = gid
        self.db = gid.db
        self.acc = acc
        self.alias = alias
        self.desc = desc
        self.exp = gid.gid

        
class uid:

    def __init__(self, sid, uri, filename = "", md5 = "", dtype="wget"):
        
        self.sid = sid
        self.db = sid.db

        self.uri = uri        
        self.desc = sid.desc
        self.md5 = md5
        self.dtype = dtype
        self.filename = filename



class dgeo:

    def __init__(self, query = "", grep = "", vgrep = "", verbose=0):
        self.dbs = ["arrayexpress", "ena", "ena_study", "gds", "sra"]
        self.found_id_db = ""
        self.grep = grep
        self.vgrep = vgrep
        self.query = query
        self.verbose = verbose
        self.md5sum = []
        #self.finfo = open (query + ".txt", "w")
        


    def get_xml (self, url, raw = False):
        #requests
        #print >>sys.stderr, "requesting", url
        xml = requests.get (url)
        return xml if raw else xml.text

    
    def query2expId (self, query, dbs = ""):

        if dbs != "":
            self.dbs = [dbs]

        self.found_in_db = ""

        gid_array = []

        for db in self.dbs:

            print ("...in", db, file=sys.stderr)
           
            if db == "gds":
                
                esearch = self.get_xml ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=(gse%5BEntry%20Type%5D)%20AND%20({})&retmode=xml".format (query))

            elif db == "sra": 

                esearch = self.get_xml ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={}&retmode=xml".format (query))

            elif db == "ena":
                            
                esearch = self.get_xml (f'https://www.ebi.ac.uk/ena/portal/api/search?result=read_study&query=accession="{query}" OR sample_alias="{query}" OR study_alias="{query}"')
                                 
       
            elif db == "ena_study":
                esearch = self.get_xml ("http://www.ebi.ac.uk/ena/data/search?query={}&result=study&display=xml".format (query))
                # PROJECT STUDY: PR..... TODO

       

            elif db == "arrayexpress":
                # TODO
                continue

            
            
            if db == "ena":

                print (esearch)
                for r in esearch.split ("\n"):

                    cols = r.strip().split("\t")
                    if cols[0] == "study_accession":
                        continue

                    if cols[0] != "" and len (cols) > 1:

                        print ("...found in", db, ": ", cols[1], file=sys.stderr)
                        
                        gid_array.append (gid (db, cols[0], cols[1])) #r.attrs['accession'], r.attrs['alias']))

                        self.found_in_db = db

                              
                # for r in search_xml.findAll("STUDY"):

                #     print ("...found in", db, "!", r.attrs, file=sys.stderr)
                    
                #     gid_array.append (gid (db, r.attrs['accession'], r.attrs['alias']))

                #     self.found_in_db = db

            elif db == "ena_study":

                search_xml = BS(esearch, "xml")

                for r in search_xml.findAll("PROJECT_LINK"):

                    print ("...found in", db, "!", r.attrs, file=sys.stderr)
                    
                    rDB = r.find("DB")

                    if rDB:
                        
                        #print (rDB)
                        accs = rDB.text
                        if accs == "ENA-RUN":
                            self.found_in_db = db
                            gid_array.append (gid (db, query))



            

            elif db == "gds" or db == "sra":

                search_xml = BS(esearch, "xml")

                for r in search_xml.findAll("Id"):
                    
                    
                    print ("...found in", db, "!", r.string, file=sys.stderr)
                    
                    gid_array.append (gid (db, r.string))
                                        
                    self.found_in_db = db
                    
                    
            if self.found_in_db != "":
                return (gid_array)

        return (gid_array)



    def expId2SampleId (self, gid):

        sid_array = []

        if gid.db == "ena_study":

            sid_array.append (sid (gid, gid.gid, alias=""))


        elif gid.db == "ena":



            esearch = self.get_xml (f'https://www.ebi.ac.uk/ena/portal/api/search?query=study_accession="{gid.gid}"&result=read_run')

            for r in esearch.split ("\n"):

                cols = r.strip().split("\t")
                if cols[0] == "run_accession":
                   continue

                if cols[0] != "" and len (cols) > 1:
                   sid_array.append (sid(gid, cols[0], alias="", desc=cols[1]))  

            

            # esearch = self.get_xml ("http://www.ebi.ac.uk/ena/data/search?query={}&result=read_run&display=xml".format (gid.gid))
                          
            # search_xml = BS(esearch, "xml")

            # for r in search_xml.findAll("RUN"):

            #     #print >>sys.stderr, "...ENA RUN!", r.attrs

            #     acc = r.attrs['accession']
                
            #     sid_array.append (sid(gid, acc, alias=r.attrs['alias'], desc=r.TITLE.text))  



        elif gid.db == "gds" or gid.db == "sra":

            esummary_gds = self.get_xml ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db={}&id={}&retmode=xml".format (gid.db, gid.gid))   
            esummary_gds_xml = BS(esummary_gds, "xml")
                            
            summary_title = "N/A"
            
            if gid.db == "gds":   
            
                summary_title = esummary_gds_xml.find("Item", {"Name" : "title"}).string
                
                gse = esummary_gds_xml.find("Item", {"Name" : "Accession"}).string
                                    
                if gse[:3] == "GSE": # Parent id
                    
                    current_gse = gse
                    print ("#ACCESSION: " + gse, gid.alias, summary_title.encode('utf-8'), file=sys.stderr)
                    
                    gid.alias = gse  #working, multiple gse's?
                    samples = esummary_gds_xml.findAll ("Item", {"Name" : "Sample"})
                    
                    for sample in samples:
                            
                        gsm_acc = sample.find ("Item", {"Name" : "Accession"}).string
                        gsm_title = sample.find ("Item", {"Name" : "Title"}).string
                        

                        sid_array.append (sid(gid, gsm_acc, alias=gse, desc=gsm_title))  
                                                

               
        return sid_array


    def get_ena_uri (self, sid):

        url_ena = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes".format (sid.acc)
            

        uri_array = []
       
        ena = self.get_xml (url_ena)
            
        for line in ena.split ("\n"):
                                            
            cols = line.split ("\t")
            
            if cols[0] == "run_accession" or len (cols) < 3:
                continue
            
            fastq = cols[1]
            fastq_md5 = cols[2]
                            
            for fq, fqm in zip(fastq.split (";"), fastq_md5.split (";")):
                
                uri_array.append (uid(sid=sid, uri=fq, md5=fqm, dtype="wget"))

        return uri_array

    def sampleId2Uri (self, sid):

        uri_array = []

        # --------------------------------------
        # ENA
        # --------------------------------------

        if sid.db == "ena" or sid.db == "ena_study":
            
            uri_array = self.get_ena_uri (sid)

            if len(uri_array) > 0: 
                return uri_array

        # --------------------------------------
        # GDS / SRA
        # --------------------------------------

        if sid.db == "gds" or sid.db == "sra":

                        
            esearch_sra = self.get_xml ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={}&term={}&retmode=xml".format (sid.db, sid.acc))
            search_sra_xml=BS(esearch_sra, "xml")

            #print "#", search_sra_xml.findAll("Id")
                
            for r in search_sra_xml.findAll("Id"):

                gid_sra = r.string

                print (gid_sra)
                
                efetch_sra =  self.get_xml ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={}&id={}&retmode=xml".format (sid.db, gid_sra))
                
                efetch_sra_xml=BS(efetch_sra, "xml")
                
                sref = efetch_sra_xml.find ("STUDY_REF")

                ref = ""
                
                if "refname" in sref:
                    ref = sref["refname"]
                

                #print "#{} of {}...".format (igsm, len(gsm_array))
                #print "#TITLE: " + efetch_sra_xml.find ("TITLE").string
                    
                for srr in efetch_sra_xml.findAll("RUN"):
                
                    gse = sid.exp
                    output = gse
                    
                    if ref != gse and ref != "":
                        output = gse +"/" + ref
                    
                    sra_xml_title = efetch_sra_xml.find ("TITLE").string
                    
                    #print sra_xml_title.encode('utf-8')                    
                    #print >>sys.stderr, "#" + "\t".join ([ref, gid_sra, sid.exp, sid.acc, srr["accession"]])
                    #print >>self.finfo, "\t".join ([ref, gid_sra, sid.exp, sid.acc, srr["accession"]])


                    sraid = srr["accession"]

                    sid.acc = sraid

                    #from_ena = self.get_ena_uri (sid)

                    sra = []
                    fastq = []

                    for sraf in srr.findAll("SRAFile"):
                        
                        if sraf['semantic_name'] == "fastq":
                            fastq.append ({'filename':sraf['filename'], 'url':sraf['url']})
                        else:
                            sra.append  ({'filename':sraf['filename'], 'url':sraf['url']})


                    # if len (from_ena) > 0:
                    #     uri_array.extend (from_ena)
                    #     print >>sys.stderr, "From ena"

                    # else:

                    use = fastq

                    if len (fastq) == 0:
                        use = sra
                    	
                    for u in use:
                        #uri_array.append (uid (sid = sid, uri = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/{}/{}/{}.sra".format (output, sraid[:6], sraid)))
                        uri_array.append (uid (sid = sid, uri = u['url'], filename=u['filename']))

                    #command.append (["N/A", "N/A", "fastq-dump --split-3 --gzip -O {} -A {}".format (output, srr["accession"])  +  " #" + efetch_sra_xml.find ("TITLE").string.encode('utf-8')])

         



        return uri_array
        #for id in sampleIds:

    def download (self, uri, download = False, foutput=None):
        
                
        folder = uri.sid.exp

        if uri.sid.gid.alias[:3] == "GSE":
            folder = uri.sid.gid.alias

        folder += "/fastq"
                
        # grep here !?

        url = uri.uri

        if url.find ("//") == -1:
            url = "ftp://" + url


        if uri.filename == "":
            cmd = "wget -P {} {} #{}".format (folder, url, uri.desc)
        else:
            cmd = "wget -P {} -O {} {} #{}".format (folder, uri.filename, url, uri.desc)
        
        print ("...", uri.desc, file=sys.stderr)

        if self.verbose > 0:
            print (cmd, file=sys.stderr)

        if foutput != None:
            print (cmd, file=sys.stderr)

        if (self.grep == "" or cmd.find (self.grep) != -1) and (self.vgrep == "" or cmd.find (self.vgrep) == -1): 




            
            if download:

                os.system (cmd)

                if uri.md5 != "":

                    md5 = folder + "/md5sum.txt"             

                    f = open (md5, "a")

                    self.md5sum.append (md5)

                    print ("\t".join ([uri.md5, os.path.basename (uri.uri)]), file=f)
                
                    f.close()

            if download or args.info != "":

            #if (self.grep == "" or uri.desc.find (self.grep) != -1) and (self.vgrep == "" or uri.desc.find (self.vgrep) == -1): 
                
                
                finfo = open (folder + "/info.txt" if args.info == "" else args.info, "a")
        
                #write info
                print ("\t".join ([uri.sid.gid.gid, uri.sid.gid.alias, uri.sid.acc, uri.sid.alias,  uri.uri, uri.desc]), file=finfo)
                finfo.close()

        else:
            print ("... SKIPPED", file=sys.stderr) 
        
        print ("\t".join ([uri.sid.gid.gid, uri.sid.gid.alias, uri.sid.acc, uri.sid.alias,  uri.uri, uri.desc]))
        
        #print "\t".join ([os.path.basename (uri.uri), uri.sid.acc, uri.sid.alias, uri.desc])
              




cwd = os.getcwd()



print ("#Searching for {}".format (args.query), file=sys.stderr)

foutput = None
if args.output != "":
    foutput = open (args.output, "w")
    print ("#Commands for downloading", args.query, file=sys.stderr)



dg = dgeo (query = args.query, grep=args.grep, vgrep=args.vgrep, verbose=args.verbose)

for eid_entry in dg.query2expId (args.query, dbs=args.db):

    #print >>sys.stderr, "eid", id

    sids = dg.expId2SampleId (eid_entry)

    
    for isid, sid_entry in enumerate (sids):

        print ("Retrieving", isid+1, "of", len(sids), "...", file=sys.stderr)
    
        #print >>sys.stderr, "sid", sid.acc


        for uid_entry in dg.sampleId2Uri (sid_entry):

            #print >>sys.stderr, "uid", uid.uri

            dg.download (uid_entry, download=args.download, foutput=foutput)



# finally check md5sum 

for md5 in list(set(dg.md5sum)):


    folder, md5_file = os.path.split(md5)
    
    os.chdir (os.getcwd() +"/" + folder)

    print (f"Checking {md5}", file=sys.stderr)
    os.system (f"md5sum -c {md5_file} > {md5_file}.report")


