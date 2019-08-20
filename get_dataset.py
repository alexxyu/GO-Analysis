#!/usr/bin/env python

import numpy as np
import pandas as pd
import sys

#Parses animal species name to ensembl dataset name
datasets = pd.read_csv('data/ensembl_datasets.txt', sep="\t")
animal_name = sys.argv[1]
print(datasets['Dataset'][datasets['Name'] == animal_name].iloc[0])