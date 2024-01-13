import pandas as pd 
# import psycopg
import re 
import os
import argparse
# import sys
# def get_all_genome_rlkdb_genes_genome(gtf_file):
parser = argparse.ArgumentParser(description='get location of genes on chromosomes')
parser.add_argument('gtf_file',help='Genome gtf annotaion file')
args = parser.parse_args()
gene_id_list=[]
genome_list=[]
with open(args.gtf_file,'r') as fr:
    for line in fr:
        line=line.strip().split('\t')
        if len(line)>1:
            if len(re.findall('transcript_id \"(.*?)\";',line[-1]))>0:
                gene_id_list.append(re.findall('transcript_id \"(.*?)\";',line[-1])[0])
                genome_list.append(line[0])
            if len(re.findall('transcript_id \"transcript:(.*?)\";',line[-1]))>0:
                gene_id_list.append(re.findall('transcript_id \"transcript:(.*?)\";',line[-1])[0])
                genome_list.append(line[0])
            if len(re.findall('protein_id \"(.*?)\";',line[-1]))>0:
                gene_id_list.append(re.findall('protein_id \"(.*?)\";',line[-1])[0])
                genome_list.append(line[0])
            if len(re.findall('transcript_id \"mRNA:(.*?)\";',line[-1]))>0:
                gene_id_list.append(re.findall('transcript_id \"mRNA:(.*?)\";',line[-1])[0])
                genome_list.append(line[0])
            if len(re.findall('gene_id \"(.*?)\";',line[-1]))>0:
                gene_id_list.append(re.findall('gene_id \"(.*?)\";',line[-1])[0])
                genome_list.append(line[0])
            if len(re.findall('gene_name \"(.*?)\";',line[-1]))>0:
                gene_id_list.append(re.findall('gene_name \"(.*?)\";',line[-1])[0])
                genome_list.append(line[0])
result=pd.DataFrame({'gene_id':gene_id_list,'genome':genome_list})
result=result.drop_duplicates()
result.to_csv('genes_genome.csv',index=None,header=None)
# if __name__=='__main__':
#     get_all_genome_rlkdb_genes_genome(sys.argv[1])
