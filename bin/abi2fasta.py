from Bio import SeqIO

import sys


in_f = sys.argv[1] 

out_f = sys.argv[2]

#filename = rec.split('/')[-1].replace('abi','fasta').replace('ab1','fasta')
    


record = SeqIO.read(in_f, "abi")
    
    
with open(out_f,'w')as o:
    SeqIO.write(record,o,'fasta')

