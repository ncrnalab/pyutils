
"""
Module for handling sqlite file (or mySQL databases...)

"""


import sys
import os

class sqlutil (object):

    def __init__ (self, type, url, user = "", pw = "",  db = ""):

        self.type = type
        self.url = url
        self.user = user
        self.pw = pw
        self.db = db
        self.cursor = None

        print ("Connecting to ", self.url, file=sys.stderr)

        try:

            if type == "mysql":

                try:
                    import MySQLdb as mdb
                    import MySQLdb.cursors
                except ImportError:
                    print ("No MySQLdb module", file=sys.stderr)
                    sys.exit()

                self.connection = mdb.connect(self.url, self.user, self.pw, self.db) #, cursorclass=cursors.SSDictCursor)
            
                with self.connection:
                    self.cursor = self.connection.cursor(MySQLdb.cursors.DictCursor)
                    #self.cursor = self.connection.cursor()
                    print ("Connected to", self.type, "database", file=sys.stderr)
            
            elif type == "sqlite" or type == "sqlite3":
                try:
                    import sqlite3
                except ImportError:
                    print ("No sqlite3 module", file=sys.stderr)
                    sys.exit()

                try:
                    self.connection = sqlite3.connect(self.url)
                    self.connection.row_factory = sqlite3.Row        ## Here's the magic!

                    self.cursor = self.connection.cursor()
                    print ("Connected to", self.type, "database", file=sys.stderr)
                
                except sqlite3.Error as e:
                    print ("Database error: {}".format (e), file=sys.stderr)
                
                            
        except:

            print ("Connection error", url, file=sys.stderr)
      

    def get_columns (self, table):
 
        
        self.query ('select * from ' + table + ' LIMIT 1')

        return [description[0] for description in self.cursor.description]





    def get_column2index (self):

        columns = self.cursor.description 
        result = [column[0] for column in columns] 
        index = 0
        self.col = {}
        for r in result:
            self.col[r] = index
            index+=1

    def get_headers (self):

        return [description[0] for description in self.cursor.description]

    def query (self, query):
        
        self.cursor.execute (query)
        self.get_column2index ()
        self.results = self.cursor.fetchall()

    def execute (self, query):
        
        self.query (query)

    def fetchall (self):

        for r in self.results:
            yield r

    def fetchone (self):

        for r in self.results:
            return r
            #return
        #return self.results[0]

    def get_db (self):  # get db name from sqlite path

        base_split = os.path.basename(self.url).split (".")

        return "_".join (base_split[:-1])
       
            

    def get (self, row, col):
        
        try:
            return row[col]
        except:
            return row[self.col[col]]

    
    def rowcount (self):

        return len(self.results)


    def get_table (self):

        print ("getting table", file=sys.stderr)
        
        file_title = os.path.splitext (os.path.basename (self.url))[0]
        
        self.execute ("SELECT name FROM sqlite_master WHERE type='table'")

        for r in self.fetchall ():
            #print (r, file=sys.stderr)

            if r['name'][-4:] != "info" and r['name'] != 'sqlite_sequence':
                return r['name']
        
                
        return file_title.replace (".", "_")


    def get_all_tables (self):

        print ("getting all tables", file=sys.stderr)
        
        self.execute ("SELECT name FROM sqlite_master WHERE type='table'")

        for r in self.fetchall ():
            yield r['name']
        
                

    def close (self):

        self.cursor.close()