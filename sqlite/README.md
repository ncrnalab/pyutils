# UCSC mySQL to sqlite3

To convert UCSC SQL tables to local sqlite files, first dump SQL file using mysqldump (here using hg38 repeatmask table as an example):

*Download sql*
 
``` bash

mysqldump --skip-lock-tables --skip-extended-insert --compact -u genome -h genome-mysql.soe.ucsc.edu hg38 rmsk > hg38.rmsk.sql

```

Then, use the available bash-script to convert to SQL to sqlite 

```
./sql2sqlite.sh hg38.rmsk.sql | sqlite3 hg38.rmsk.sqlite

```

make sure sql2sqlite is executable (chmod +x sql2sqlite)