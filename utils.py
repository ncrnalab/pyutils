"""
Few common python functionalities

"""

import sys


COMPLEMENT = {
    'a' : 't',
    't' : 'a',
    'c' : 'g',
    'g' : 'c',
    'k' : 'm',
    'm' : 'k',
    'r' : 'y',
    'y' : 'r',
    's' : 's',
    'w' : 'w',
    'b' : 'v',
    'v' : 'b',
    'h' : 'd',
    'd' : 'h',
    'n' : 'n',
    'A' : 'T',
    'T' : 'A',
    'C' : 'G',
    'G' : 'C',
    'U' : 'A',    
    'K' : 'M',
    'M' : 'K',
    'R' : 'Y',
    'Y' : 'R',
    'S' : 'S',
    'W' : 'W',
    'B' : 'V',
    'V' : 'B',
    'H' : 'D',
    'D' : 'H',
    'N' : 'N',
}


def complement(s):
    return "".join([COMPLEMENT[x] for x in s])

def rev_comp(seq):
    return complement(seq)[::-1]

	
	
codontable = {
       'ATA':'I',
       'ATC':'I',
       'ATT':'I',
       'ATG':'M',
       'ACA':'T',
       'ACC':'T',
       'ACG':'T',
       'ACT':'T',
       'AAC':'N',
       'AAT':'N',
       'AAA':'K',
       'AAG':'K',
       'AGC':'S',
       'AGT':'S',
       'AGA':'R',
       'AGG':'R',
       'CTA':'L',
       'CTC':'L',
       'CTG':'L',
       'CTT':'L',
       'CCA':'P',
       'CCC':'P',
       'CCG':'P',
       'CCT':'P',
       'CAC':'H',
       'CAT':'H',
       'CAA':'Q',
       'CAG':'Q',
       'CGA':'R',
       'CGC':'R',
       'CGG':'R',
       'CGT':'R',
       'GTA':'V',
       'GTC':'V',
       'GTG':'V',
       'GTT':'V',
       'GCA':'A',
       'GCC':'A',
       'GCG':'A',
       'GCT':'A',
       'GAC':'D',
       'GAT':'D',
       'GAA':'E',
       'GAG':'E',
       'GGA':'G',
       'GGC':'G',
       'GGG':'G',
       'GGT':'G',
       'TCA':'S',
       'TCC':'S',
       'TCG':'S',
       'TCT':'S',
       'TTC':'F',
       'TTT':'F',
       'TTA':'L',
       'TTG':'L',
       'TAC':'Y',
       'TAT':'Y',
       'TAA':'*',
       'TAG':'*',
       'TGC':'C',
       'TGT':'C',
       'TGA':'*',
       'TGG':'W',}


def translate(sequence, orf=True):  # if orf = False, just translate entire seq
    
    proteinsequence = ''
    
    if orf:
       start = sequence.find('ATG')
       if start == -1:
         return ""
    else:
      start=0
       
    sequencestart = sequence[int(start):].upper()
    
    #print start, sequencestart
    
    for n in range(0,len(sequencestart),3):
        #print n, sequencestart[n:n+3] in codontable        
        if sequencestart[n:n+3] in codontable:
            codon = codontable[sequencestart[n:n+3]]
            proteinsequence += codon
            if codon == "*" and orf:
               return proteinsequence
            

    return proteinsequence


    