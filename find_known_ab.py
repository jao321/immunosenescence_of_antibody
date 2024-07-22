import pandas as pd
from anarci import number
import glob
from Bio.Seq import Seq
from Bio import SeqIO
import os

def renumber(seqAA):
    try:
        numbered = number(seqAA)[0]
        cdr3 = [numbered[x] for x in range(len(numbered)) if (numbered[x][0][0]>=104) and (numbered[x][0][0]<=117)]
        cdr3 = ''.join([cdr3[x][1] for x in range(len(cdr3))]).replace('-','')
        cdr2 = [numbered[x] for x in range(len(numbered)) if (numbered[x][0][0]>=54) and (numbered[x][0][0]<=66)]
        cdr2 = ''.join([cdr2[x][1] for x in range(len(cdr2))]).replace('-','')
        cdr1 = [numbered[x] for x in range(len(numbered)) if (numbered[x][0][0]>=23) and (numbered[x][0][0]<=40)]
        cdr1 = ''.join([cdr1[x][1] for x in range(len(cdr1))]).replace('-','')
        # rep.loc[rep["sequence_id"]==ID]["sequence"].item()
        return cdr3+","+cdr2+","+cdr1
    except:
        return ",,"


repertoires = glob.glob("path/to/folder/with/tsv_YClon_clonotyped_files/**/*YClon_clonotyped.tsv")

general_df = pd.DataFrame()
for i in repertoires:
    df = pd.read_csv(i,sep="\t")
    df = df[['sequence_id', 'sequence','clone_id']]
    df['sample_ID'] = i.split('/')[-1].replace('.tsv','')
    try:
        df['AA'] = df['sequence'].apply(lambda x: str(Seq(x).translate()))
        df['cdr_renumbered'] = df['AA'].apply(lambda x: renumber(x))
    except:
        print(i)
        continue
    df['cdr1_renumb'] = df['cdr_renumbered'].str.split(',').str[2]
    df['cdr2_renumb'] = df['cdr_renumbered'].str.split(',').str[1]
    df['cdr3_renumb'] = df['cdr_renumbered'].str.split(',').str[0]
    general_df = pd.concat([general_df,df])

general_df = general_df.dropna(subset='cdr3_renumb')
general_df = general_df[general_df.cdr3_renumb!='']


path_knowm_ab = 'path/to/fasta/with/sequences_of_known_ab.fasta'
sequences = []

for record in SeqIO.parse(path_knowm_ab,"fasta"):
    sequences.append([record.id,str(record.seq)])
known_ab_df = pd.DataFrame(sequences, columns=["id","sequence"])
known_ab_df['AA'] = known_ab_df['sequence'].apply(lambda x: str(Seq(x).translate()))
known_ab_df['cdr_renumbered'] = known_ab_df['sequence'].apply(lambda x: renumber(str(Seq(x).translate())))
known_ab_df['cdr1_renumb'] = known_ab_df['cdr_renumbered'].str.split(',').str[2]
known_ab_df['cdr2_renumb'] = known_ab_df['cdr_renumbered'].str.split(',').str[1]
known_ab_df['cdr3_renumb'] = known_ab_df['cdr_renumbered'].str.split(',').str[0]


plabdab = pd.read_csv('plabdab_result.csv')
plabdab['cdr_renumbered'] = plabdab['heavy_sequence'].apply(lambda x: renumber(x))
plabdab['cdr1_renumb'] = plabdab['cdr_renumbered'].str.split(',').str[2]
plabdab['cdr2_renumb'] = plabdab['cdr_renumbered'].str.split(',').str[1]
plabdab['cdr3_renumb'] = plabdab['cdr_renumbered'].str.split(',').str[0]
plabdab = plabdab[plabdab['targets_mentioned']!='Unidentified']
plabdab = plabdab.dropna(subset='cdr3_renumb')

cdr1_cdr2_cdr3_match_plabdab = pd.merge(general_df,plabdab,on=['cdr1_renumb','cdr2_renumb','cdr3_renumb'])
cdr3_match_plabdab = pd.merge(general_df,plabdab,on=['cdr3_renumb'])

cdr1_cdr2_cdr3_match_known = pd.merge(general_df,known_ab_df,on=['cdr1_renumb','cdr2_renumb','cdr3_renumb'])
cdr3_match_known = pd.merge(general_df,known_ab_df,on=['cdr3_renumb'])







