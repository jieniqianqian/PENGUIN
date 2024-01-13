import Bio.SeqIO
# import psycopg
import pandas as pd
import os
# import sys
import argparse
import itertools
letters=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
letters_tuple=list(itertools.permutations(letters,2))
letters_list=[]
for i in letters_tuple:
    letters_list.append(''.join(i))
def rename_new_id(df):
    global letters_list
    if len(df.index)==1:
        return df
    else:
        df.index=range(len(df.index))
        result_list=[]
        if len(df.index)<=26:
            for i in range(len(df.index)):
                result_list.append([df.iloc[i,0],df.iloc[i,1],df.iloc[i,2]+letters[i]])
        else:
            for i in range(len(df.index)):
                result_list.append([df.iloc[i,0],df.iloc[i,1],df.iloc[i,2]+letters_list[i]])
        result=pd.DataFrame(result_list,columns=['Query_accession','Target_accession','new_id'])
        return result
parser = argparse.ArgumentParser(description='Rename the gene ID of chromosome-level plant based on the results of Comparison with model plants')
parser.add_argument('seqfile',help='Genome protein annotaion')
parser.add_argument('all_genome_protein_matches_fmt6',help='The fmt6 results of Comparison with model plants')
parser.add_argument('genes_genome',help='The chromosome location of genes')
parser.add_argument('uniprot_abbreviation_upper',help='Species id abbreviation')
parser.add_argument('mode',help='Comparison model')
args = parser.parse_args()
# def get_all_genome_new_id_split_diamond_result(seqfile,all_genome_protein_matches_fmt6,genes_genome,uniprot_abbreviation_upper,mode):
gene_id_length={}
genome=Bio.SeqIO.parse(args.seqfile,'fasta')
for i in genome:
    gene_id_length[i.id]=len(str(i.seq))
Query_len=pd.DataFrame(gene_id_length,index=[0])
Query_len=Query_len.T 
Query_len.columns=['Query_len']
Query_len['Query_accession']=Query_len.index.tolist()
genes_cds=pd.read_csv(args.genes_genome,names=['Query_accession','genome'])
def get_similiar_id(df):
    df=df.loc[df['Sequence_identity']==df['Sequence_identity'].max(),:]
    df=df.iloc[0,:]
    return df
diamond_result=pd.read_table(args.all_genome_protein_matches_fmt6,names=['Query_accession','Target_accession','Sequence_identity','Length','Mismatches','Gap_openings',\
'Query_start','Query_end','Target_start','Target_end','E_value','Bit_score'])
if args.mode=='O':
    Target_len=pd.read_table('Orysa_seq_len.tsv')
else:
    Target_len=pd.read_table('Arath_seq_len.tsv')
Target_len.columns=['Target_accession','Target_len']
diamond_result=pd.merge(diamond_result,Target_len,on='Target_accession')
diamond_result=pd.merge(diamond_result,Query_len,on='Query_accession')
id_conversion_result=diamond_result.loc[(diamond_result['Sequence_identity']==100.0)&(diamond_result['Target_len']==diamond_result['Query_len'])&(diamond_result['Mismatches']==0)&(diamond_result['Gap_openings']==0)&(diamond_result['E_value']==0),:]
id_conversion_others=diamond_result.loc[~diamond_result['Query_accession'].isin(id_conversion_result['Query_accession']),:].sort_values(by=['Query_accession','Sequence_identity'],axis=0,ascending=False).groupby('Query_accession').apply(get_similiar_id)
diamond_result=pd.concat([id_conversion_result,id_conversion_others])
diamond_result.index=range(len(diamond_result.index))
identity=[]
for j in range(len(diamond_result.index)):
    identity_list=diamond_result.loc[j,['Sequence_identity','Target_len','Query_len']].values
    if identity_list[1]>identity_list[2]:
        identity.append(identity_list[0]/(identity_list[1]/identity_list[2]))
    else:
        identity.append(identity_list[0]/(identity_list[2]/identity_list[1]))
diamond_result=pd.DataFrame({'Query_accession':diamond_result['Query_accession'].values,'Id_identity':identity}).join(diamond_result[['Target_accession','Length','Mismatches','Gap_openings',\
'Query_start','Query_end','Target_start','Target_end','E_value','Bit_score']])
diamond_result.index=range(len(diamond_result.index))
id_conversion_result=pd.merge(diamond_result,genes_cds)
id_conversion_result.index=range(len(id_conversion_result.index))
id_conversion_result['Target_accession']=id_conversion_result['Target_accession'].apply(lambda x:x.split('.')[0])
new_id=[]
if args.mode=='O':
    for j in range(len(id_conversion_result.index)):
        new_id.append(args.uniprot_abbreviation_upper+str(id_conversion_result.loc[j,'genome'])+'G'+id_conversion_result.loc[j,'Target_accession'][6:8].lstrip('0')+'O'+id_conversion_result.loc[j,'Target_accession'][9:])
else:
    for j in range(len(id_conversion_result.index)):
        new_id.append(args.uniprot_abbreviation_upper+str(id_conversion_result.loc[j,'genome'])+'G'+id_conversion_result.loc[j,'Target_accession'][2]+'A'+id_conversion_result.loc[j,'Target_accession'][4:])
id_conversion_result=id_conversion_result.loc[:,['Query_accession','Target_accession']]
id_conversion_result['new_id']=new_id
id_conversion_result=id_conversion_result.groupby('new_id').apply(rename_new_id)
id_conversion_result.to_csv('change_rlkdb_id.csv',index=None)
# if __name__=='__main__':
#     get_all_genome_new_id_split_diamond_result(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])