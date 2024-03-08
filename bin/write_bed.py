import gffutils
import pysam
import sys



bam_in = sys.argv[3:]
#bam_in = sys.argv[1]


bed_out = sys.argv[1]

gffutils_db = sys.argv[2]



#gffutils_db = '/mnt/d/2023/ITD/ITD_2317_codeclub_gffutils/Code_Club/26_gffutils/GCF_000001405.40_GRCh38.p14_genomic.gff.db'



db = gffutils.FeatureDB(gffutils_db)

seq_id_dict = {}
for r in db.seqids():
    if r=='NC_000023.11':
        seq_id_dict[f'chrX'] = r
        continue
    if r=='NC_000024.10':
        seq_id_dict[f'chrY'] = r
        continue
    if 'NC' in r:
        i = r.split('.')[0].lstrip('NC_').lstrip('0')
        seq_id_dict[f'chr{i}'] = r



#seq_id_dict = {k.split('.'[0]):v for k,v in seq_id_dict.items}

        


ROI = set()

for bam_file in bam_in:
    with pysam.AlignmentFile(bam_file) as bam:
        for r in bam:
            ROI.add((r.reference_name,r.reference_start,r.reference_end))

padding = 2000

bed = {}


for region in ROI:
    c=region[0]
    s=region[1]
    e=region[2]
    if c not in bed:
        bed[c]=[set(range(s-padding,e+padding))]
    else:
        cov_reg = set(range(s-padding,e+padding))
        overlapp = [len(x.intersection(cov_reg)) for x in bed[c]]
        if sum(overlapp)>0:
            continue
        else:
            bed[c].append(cov_reg)

out = open(bed_out,'w')           
for c,rr in bed.items():
    
    if c not in seq_id_dict:
        continue
    new_c = seq_id_dict[c]
    
    for r in rr:
        #c = c.lstrip('chr')
        s = min(r)
        e = max(r)
        genes = set()
        for feat in db.all_features((new_c,s,e)):
            if feat.featuretype=='gene':
                genes.add(feat.id)
        out.write(f'{c}\t{min(r)}\t{max(r)}\t{",".join(genes)}\n')        

        
        
#gene_symbol_wanted = 'GALC'


#gene = db[f'gene-{gene_symbol_wanted}']


#with open('out.bed','w') as o:
#    o.write(f'{seq_id_dict[gene.seqid]}\t{gene.start}\t{gene.end}')