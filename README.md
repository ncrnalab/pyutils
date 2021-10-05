# pyutils

### Introduction

Various python scrips for RNAseq analysis used in [Stagsted *et al*, "The RNA-binding protein SFPQ preserves long-intron splicing and regulates circRNA biogenesis in mammals", eLIFE, 2021](https://elifesciences.org/articles/63088)

### Merge and annotate circRNAs

First, to combine output from circRNA prediction algorithm, use circM:

**circM.py**
 
``` python
Usage:
  circM.py [ARGUMENTS] > [OUTPUT]
  
Options:
  -f <files>        Input files
  -a <string>       Algorithm 
                    [acsf, circexplorer, circexplorer2, 
                    ciri, ciri2, circrna_finder, dcc 
                    find_circ, knife, mapsplice, uroborus] 
  -cf <int>         Expression threshold for each sample, back-splice spanning reads (default=2)
  -cs <int>         Merged expression threshold, backsplice-spanning reads (default=2)
```


Then, to annotate the merged list of circRNAs:

**annotate_circ.py**
 
``` python
Usage:
  annotate_circ.py [ARGUMENTS] < [INPUT BED] > [OUTPUT]
  
Options:
  -a <file>       Path to GTF annotation file (preferably from GENCODE)
```

Finally, to retrieve flanking alu distances (or other repetitive element): 


**annotate_repeat.py**
 
``` python
Usage:
  annotate_repeat.py [ARGUMENTS] < [INPUT BED] > [OUTPUT]
  
Options:
  -rdb <file>       Path to sqlite database containing the UCSC repeatmask track
  -rdf <string>     Comma-separated string of repeat families to analyze (default=Alu)
  
```

Help to retrieve and convert UCSC repeatmask tables to sqlite can be found in the *sqlite* sub-folder.


#### Example

To merge and annotate output from find_circ:

```bash
python circM.py -f *.circ_candidates.bed -a find_circ | python annotate_circ.py -a /path/to/gtf | python annotate_rep.py -rdb /path/to/sqlite > find_circ.bed
```


### Get spliced reads from bam

**get_spliced_reads.py**
 
Outputs bed-like file with all detectable introns and their corresponding read count
 
``` python
Usage:
  get_spliced_reads.py [ARGUMENTS] > [OUTPUT]
  
Mandatory arguments:
  -f <file>       Bam file
  -g <file>       GTF file (preferably from GENCODE)
  -fa <file>      Genome fasta-file (should be indexed with samtools faidx)
Optional arguments:  
  -s <string>     Library strandedness (as in htseq-count). Default='reverse'
  -t <int>        Overhang. Minimun number of aligned reads in exon. Default=8
  -mi <int>       Maximum intron length. Default=500000
```

#### Example

To retreive spliced reads from all bam:

```bash
for b in *.bam;
do
   python get_spliced_reads.py -f $b -g /path/to/gtf -fa /path/to/genome.fa > $b.spliced_reads.bed;
done;
```

### Using spliced reads for flanking linear spliced reads for circRNAs

This also provides linear reads for zero-count circRNAs, which is not typically provided by the circRNA prediction algorithms  

**get_flanking_spliced_reads.py**

``` python
Usage:
  get_flanking_spliced_reads.py [ARGUMENTS] < [INPUT BED] > [OUTPUT]
  
Mandatory arguments:
  -f <files>      All the Bam-derived spliced reads (as produced by get_spliced_reads.py)
  -s <bool>       If set, separate SA and SD spanning linear spliced reads. If not set, only one total linear count per circRNA per sample
```
#### Example

```bash
python get_flanking_spliced_reads.py -f *.spliced_reads.bed < find_circ.bed > find_circ_with_flanking_spliced_reads.bed
```


### Get alternative splicing

**get_flanking_spliced_reads.py**

``` python
Usage:
  get_alternative_splicing.py [ARGUMENTS] > [OUTPUT]
  
Mandatory arguments:
  -f <files>      All the Bam-derived spliced reads (as produced by get_spliced_reads.py)
```
#### Example

```bash
python get_alternative_splicing.py -f *.spliced_reads.bed  > alternative_splicing.bed
```




### License

Copyright (C) 2021 ncRNALab.  
