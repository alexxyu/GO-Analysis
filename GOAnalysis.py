#!/usr/bin/env python

import pandas as pd
import numpy as np
import os, glob
from scipy.stats import linregress, norm, spearmanr
from goatools import obo_parser
import sys

pd.options.mode.chained_assignment = None

#Requires 1 input: output filename
#Call: python GOAnalysis.py [Lifespan measure] [output filename]
if len(sys.argv) < 3:
    print('Error: Have to specify which lifespan measurement to use and output filename.')
    print('')
    exit()

l_symbol = sys.argv[1]
life_measure = ''
if l_symbol == 'ML':
    life_measure = 'Maximum longevity (yrs)'
elif ((l_symbol != 'NLa') and (l_symbol != 'NLn1')):
    print('Error: Not a valid measurement specified.')
    print('Maximum longevity: ML')
    print('Normalized lifespan by alpha: NLa')
    print('Normalized lifespan by -1: NLn1')
    print('')
    exit()

out_filename = sys.argv[2]

#Reads different data files relevant to analysis
raw_data = pd.read_csv('data/anage_data.txt', sep="\t")
genome = pd.read_csv('data/genomes.txt', sep='\t', encoding='latin-1')
GO_data = pd.read_csv('data/animalsGO.txt', sep='\t')

#Filters for animals with finite data values and merges dataframes together
data = raw_data[np.isfinite(raw_data['Body mass (g)'])]
data['Binomial'] = data['Genus'] + ' ' + data['Species']
data = pd.merge(data, genome, how="outer", on='Binomial')
data['Genome'] = data['Genome'].fillna(value=0)
data = pd.merge(data, GO_data, how="outer", on='Common name')

data = data[np.isfinite(data['Maximum longevity (yrs)'])] 
data = data[np.isfinite(data['Metabolic rate (W)'])]

#Calculates normalized lifespan measurement
mammals_and_birds = data[(data['Class']!='Reptilia') & (data['Class']!='Amphibia')]
slope, intercept, r_value, p_value, std_err = linregress(np.log(mammals_and_birds['Metabolic rate (W)']/mammals_and_birds['Body mass (g)']), np.log(mammals_and_birds['Maximum longevity (yrs)']))
data['NLa'] = data['Maximum longevity (yrs)'] / (data['Metabolic rate (W)']/data['Body mass (g)']) ** (-slope)
data['NLn1'] = data['Maximum longevity (yrs)'] / (data['Metabolic rate (W)']/data['Body mass (g)']) ** (-1)
mammals_and_birds = data[(data['Class']!='Reptilia') & (data['Class']!='Amphibia')]

#Filters for animals with genomes sequenced and cleans up data
genned_data = data[data['Genome']==1]
genned_data = genned_data.drop_duplicates(subset=['Common name'], keep="first")
genned_data = genned_data.reset_index(drop=True)
genned_data = genned_data.reindex([5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,4,1,2,0,3,24])

oboDat = obo_parser.GODag('data/go-basic.obo')

#Creates new dataframe with GO terms with significant Spearman correlation
path = r'countData'
path_len = len(path) + len('count_') + 1
all_files = glob.glob(path + "/*.txt")
go_df = pd.DataFrame(columns=['GO ID', 'GO Definition', 'Median Gene #', 'Spearman Coef', 'P-Value'])

if not os.path.exists('output'):
    os.makedirs('output')

for filename in all_files:

    #Reads and formats data for data export
    #Drops reptiles from dataset
    try:
        genned_data[filename[path_len:-4]] = pd.read_csv(filename, sep='\t', header = None).values
        coef, p = spearmanr(genned_data[l_symbol].drop([44,45]), genned_data[filename[path_len:-4]].drop([44,45])/genned_data['Protein_Coding'].drop([44,45]))

        go_df.loc[len(go_df)] = ["GO:"+filename[path_len:-4], oboDat["GO:"+filename[path_len:-4]].name, np.median(genned_data[filename[path_len:-4]]), coef, p]

    except:
        print('Error: %s contains a wrong number of terms. Term skipped.' % filename)

print('Found %i term(s) to output.' % len(go_df))
go_df.to_csv('output/' + out_filename, sep='\t')
