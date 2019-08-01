import pandas as pd
import numpy as np
import os, glob
from scipy.stats import linregress, norm, spearmanr

raw_data = pd.read_csv('data/anage_data.txt', sep="\t")
genome = pd.read_csv('data/genomes.txt', sep='\t', encoding='latin-1')
GO_data = pd.read_csv('data/animalsGO.txt', sep='\t')

data = raw_data[np.isfinite(raw_data['Body mass (g)'])]
data['Binomial'] = data['Genus'] + ' ' + data['Species']
data = pd.merge(data, genome, how="outer", on='Binomial')
data['Genome'] = data['Genome'].fillna(value=0)
data = pd.merge(data, GO_data, how="outer", on='Common name')

data = data[np.isfinite(data['Maximum longevity (yrs)'])]
data = data[np.isfinite(data['Metabolic rate (W)'])]

mammals_and_birds = data[(data['Class']!='Reptilia') & (data['Class']!='Amphibia')]
alpha, intercept, r_value, p_value, std_err = linregress(np.log(mammals_and_birds['Metabolic rate (W)']/mammals_and_birds['Body mass (g)']), np.log(mammals_and_birds['Maximum longevity (yrs)']))
data['NL'] = data['Maximum longevity (yrs)'] / (data['Metabolic rate (W)']/data['Body mass (g)']) ** (-alpha)
mammals_and_birds = data[(data['Class']!='Reptilia') & (data['Class']!='Amphibia')]

genned_data = mammals_and_birds[mammals_and_birds['Genome']==1]
genned_data = genned_data[np.isfinite(genned_data['Total'])]
genned_data = genned_data.drop(198)

#filename = "count_0000002.txt"

path = r'/Users/alexyu/Downloads/kaks/countData'
all_files = glob.glob(path + "/*.txt")
go_df = pd.DataFrame(columns=['GO ID', 'Median Gene #', 'Spearman Coef', 'P-Value'])
for filename in all_files:
    
    term = pd.read_csv(filename, sep='\t', header=None).values

    genned_data[filename[45:-4]] = pd.read_csv(filename, sep='\t', header=None).values[:-2]
    coef, p = spearmanr(genned_data['NL'], genned_data[filename[45:-4]]/genned_data['Total'])

    if p < 0.1:
        go_df.loc[len(go_df)] = ["GO:"+filename[45:-4], np.median(term), coef, p]

go_df.to_csv('sigGOTerms.txt', sep='\t')