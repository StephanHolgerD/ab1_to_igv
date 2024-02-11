import os
import pandas as pd
from glob import glob
from collections import Counter
import shutil
DIRECTION=["1","2"]
GROUPS=set()
sample_df=pd.read_csv("samples.tsv", sep="\t",dtype=object)

SAMPLES =sample_df['ID']
genome = '/home/stephano/Documents/02_gitHubProjects/00_testData/GRCh38_full_analysis_set_plus_decoy_hla.mmi'

bwt2_index = '/home/stephano/Documents/02_gitHubProjects/00_testData/GRCh38_full_analysis_set_plus_decoy_hla_btw2/GRCh38_full_analysis_set_plus_decoy_hla'

def check_symlink(file1, file2):
    try:
        shutil.copytree(file1,file2)
        #os.symlink(file1, file2)
    except FileExistsError:
        print("Link for " + file1 + " is already present in 01_raw")





for s,d in zip(sample_df['folder'],sample_df['ID']):
    os.makedirs("../01_data/",exist_ok=True)
    check_symlink(s,f'../01_data/{d}')
 #   os.makedirs("../02_trimmed/" + b +"/fastqc", exist_ok=True)
 #   os.makedirs("../03_mapped/"+b, exist_ok=True)







ab1_folder, ab1_file = glob_wildcards("../01_data/{tumor}/{sample}.ab1")



fasta_folder, fasta_file = glob_wildcards("../01_data/{tumor}/{sample}.fasta")




ab1_p = [f'{x}/{xx}' for x,xx in zip(ab1_folder,ab1_file)]


fasta_p = [f'{x}/{xx}' for x,xx in zip(fasta_folder,fasta_file)]



rule all:
    input:
        expand('../02_fasta/{sample}.fasta',sample=ab1_p),
        '../03_merge_fasta/primer.fasta',
        '../03_merge_fasta/ab1.fasta',
        '../04_mapped/primer.bam',
        '../04_mapped/ab1.bam',
        expand('../05_result/{sample}/{sample}.bam',sample=SAMPLES),
        expand('../05_result/{sample}/{sample}.bam',sample=SAMPLES)
        
        



rule ab1_to_fasta:
    input:
        ab1_f = '../01_data/{sample}.ab1'
    threads: 2
    priority: 50
    output:
        fasta_f = '../02_fasta/{sample}.fasta'
        
    conda:
        "envs/bio.yaml"
    shell:
        'python bin/abi2fasta.py {input.ab1_f} {output.fasta_f}'

rule create_fasta1:
    input:
        fasta_f = expand('../02_fasta/{sample}.fasta',sample=ab1_p)


    output:
        allfasta = '../03_merge_fasta/ab1.fasta'

    shell:
        'cat {input.fasta_f}  >>{output.allfasta}'


rule create_fasta2:
    input:
        fasta_primer = expand('../01_data/{sample}.fasta',sample=fasta_p)


    output:
        allfasta = '../03_merge_fasta/primer.fasta'

    shell:
        'cat {input.fasta_primer} >>{output.allfasta}'



rule bwt2_map:
    input:
        allfasta = '../03_merge_fasta/primer.fasta'
    conda:
        "envs/bwt2_samtools.yaml"
    output:
        sam = '../04_mapped/primer.sam',
        bam = '../04_mapped/primer.bam',

    shell:
        'bowtie2 -L 4 --gbar 1 -R 10 -D 50 -f -x {bwt2_index} -U  {input.allfasta} -S {output.sam}; \
        samtools view -b  {output.sam} | samtools sort /dev/stdin >  {output.bam};samtools index {output.bam}' 

rule minimap2_map:
    input:
        allfasta = '../03_merge_fasta/ab1.fasta'
    conda:
        "envs/minimap2_samtools.yaml"
    output:
        bam = '../04_mapped/ab1.bam',
        sam = '../04_mapped/ab1.sam',
    shell:
        'minimap2 -ax splice -a {genome} {input.allfasta} > {output.sam} ;\
            samtools sort -O BAM {output.sam}  >  {output.bam} ; samtools index {output.bam}' 
        
        
        

rule split_bam:
    input:
        primer_bam = '../04_mapped/primer.bam',
        ab1_bam = '../04_mapped/ab1.bam',
    output:
        result_bam = '../05_result/{sample}/{sample}.bam'
    conda:
        "envs/bio_pysam.yaml"
    shell:
        'python bin/split_bam.py {input.primer_bam} {input.ab1_bam} {output.result_bam} ../01_data/{wildcards.sample}/ \
            ../02_fasta/{wildcards.sample}/'