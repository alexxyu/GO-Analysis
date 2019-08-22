#!/usr/bin/env python

import pandas as pd
import glob
import os

GO_data = pd.read_csv('../data/dataSpecies.txt', sep='\t', names=['Animal'])

path_to_count_data = r'countData/'
count_files = glob.glob(os.path.join(path_to_count_data, '*.txt'))
for file in count_files:
    
    term_data = pd.read_csv(file, sep='\t', names=['GO:'+file[len(path_to_count_data)+6:-4]])
    GO_data = GO_data.join(term_data)

GO_data_t = GO_data.T
GO_data_t.columns = GO_data_t.iloc[0]
GO_data_t = GO_data_t.iloc[1:]
GO_data_t = GO_data_t.sort_index()
GO_data_t = GO_data_t.rename_axis(None, axis=1).rename_axis('GO Term', axis=0)
GO_data_t = GO_data_t.reindex(sorted(GO_data_t.columns), axis=1)
GO_data_t.to_csv('output/GO_count_data.txt', sep='\t')
