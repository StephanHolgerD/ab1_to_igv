
import pysam
from glob import glob
import sys
from Bio import SeqIO


bam_in1 = sys.argv[1]
bam_in2 = sys.argv[2]


bam_out = sys.argv[3]
bam_primer_out = sys.argv[4]





fastas = sys.argv[5:]


fastas = [x for xx in fastas for x in glob(f'{xx}/*fasta')]

print(fastas)
fasta_ids = []


for f in fastas:
    x= SeqIO.parse(f,'fasta')
    for rec in x:
        fasta_ids.append(rec.id)


done_ids = set()
with pysam.AlignmentFile(bam_in1)as bam1,pysam.AlignmentFile(bam_in2)as bam2:
    with pysam.AlignmentFile(bam_out,'wb',template=bam1) as out_bam:
        for r in bam2.fetch():
            if r.qname in fasta_ids and r.qname not in done_ids:
                out_bam.write(r)
                done_ids.add(r.qname)
    
    with pysam.AlignmentFile(bam_primer_out,'wb',template=bam1) as out_bam:

        for r in bam1.fetch():
            if r.qname in fasta_ids and r.qname not in done_ids:
                done_ids.add(r.qname)
                out_bam.write(r)        
    