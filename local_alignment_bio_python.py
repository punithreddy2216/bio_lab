from Bio import pairwise2
from Bio.pairwise2 import format_alignment

seq1="ACACACTA"
seq2="AGCACACA"

alignments=pairwise2.align.localms(seq1,seq2,2,-1,-1,-1)

for alignment in alignments:
    print(format_alignment(*alignment))
    break