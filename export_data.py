import pandas as pd
import numpy as np
import glob
import os
import sys

animal_name1 = sys.argv[1]
animal_name2 = sys.argv[2]
queryid = sys.argv[3]
outFile = animal_name1 + "_" + animal_name2 + "_" + queryid + ".txt"

#Writing pipeline output to one file
path = r'output'
all_files = sorted(glob.glob(path + "/*.kaks"))
all_files = sorted(all_files, key=lambda name: int(name[13:-5]))

#Concatenates all the Ka/Ks data for each gene pair into one dataframe and saves it to file
li = []
for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=0, sep='\t')
    li.append(df)
data = pd.concat(li, axis=0, ignore_index=True)
data.to_csv('kaksData/' + outFile, sep='\t')

sig_data = data[(data['Ka/Ks'] > 1) & (data['P-Value(Fisher)'] < 0.1)]

#Outputs genes with evidence of positive selection to output file and formats parameters
gene_df = pd.DataFrame(columns=['Protein', 'Symbol', 'Pair Number', animal_name1 + 'Protein', animal_name2 + 'Protein','Ka/Ks', 'P-Value'])
for index, row in sig_data.iterrows():
    
    file = open("OMA/Output/OrthologousGroupsFasta/OG{0}.fa".format(index), mode='r')
    text = file.read()
    
    seqs = text.split('>')
    seq1 = seqs[1]
    seq2 = seqs[2]

    desc1 = seq1.split('|')
    gene_name = desc1[4]
    prot_id1 = desc1[0]

    desc2 = seq2.split('|')
    prot_id2 = desc2[0]

    stop = 0
    prot_desc = ''
    if len(desc1[5]) > 0:
        stop = desc1[5].find(' [')
        prot_desc = desc1[5][:stop]
    elif len(desc2[5]) > 0:
        stop = desc2[5].find(' [')
        prot_desc = desc2[5][:stop]
    
    gene_df.loc[len(gene_df)] = [prot_desc, gene_name, index, prot_id1, prot_id2, row['Ka/Ks'], row['P-Value(Fisher)']]

gene_df['Pair Number'] = gene_df['Pair Number'].astype('int64')
gene_df.set_index('Pair Number', inplace=True)

gene_df = gene_df.dropna()
gene_df.to_csv('filteredGenes/' + outFile, sep='\t')