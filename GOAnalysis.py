#!/usr/bin/env python

import pandas as pd
import numpy as np
import os, glob
from scipy.stats import linregress, norm, spearmanr
from goatools import obo_parser
import sys

pd.options.mode.chained_assignment = None

#Requires 2 inputs: alpha value and output filename
#Call: python GOAnalysis.py [alpha value] [output filename]
if len(sys.argv) < 3:
    print('Error: Have to specify alpha value and output filename.')
    exit()

alpha = float(sys.argv[1])
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
data['NL'] = data['Maximum longevity (yrs)'] / (data['Metabolic rate (W)']/data['Body mass (g)']) ** (-alpha)
mammals_and_birds = data[(data['Class']!='Reptilia') & (data['Class']!='Amphibia')]

#Filters for animals with genomes sequenced and cleans up data
genned_data = mammals_and_birds[mammals_and_birds['Genome']==1]
genned_data = genned_data[np.isfinite(genned_data['Total'])]
genned_data = genned_data.drop(198)

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
    try:
        genned_data[filename[path_len:-4]] = pd.read_csv(filename, sep='\t', header = None).values[:-2]
        coef, p = spearmanr(genned_data['NL'], genned_data[filename[path_len:-4]]/genned_data['Total'])

        #Adds terms only if it is significant
        if p < alpha and coef > 0:
            go_df.loc[len(go_df)] = ["GO:"+filename[path_len:-4], oboDat["GO:"+filename[path_len:-4]].name, np.median(genned_data[filename[path_len:-4]]), coef, p]
    except:
        print('Error: %s contains a wrong number of terms. Term skipped.' % filename)

print('Found %i term(s) to be significant.' % len(go_df))
go_df.to_csv('output/' + out_filename, sep='\t')
