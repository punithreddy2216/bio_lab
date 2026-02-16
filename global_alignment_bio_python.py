from Bio import pairwise2
from Bio.pairwise2 import format_alignment

seq1="GATTACA"
seq2="GCATGCU"

alignments=pairwise2.align.globalms(seq1,seq2,1,-1,-2,-2)
print(len(alignments))
for an in alignments[:1]:
    print(format_alignment(*an))
print("helo")