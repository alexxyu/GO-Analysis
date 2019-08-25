#! /usr/bin/env python

from ete3 import Tree, TreeStyle
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, hex2color, rgb2hex
import sys
import pandas as pd
import numpy as np

def layout(node):
    
    global use_detail
    if node.is_leaf():
        
        name = node.name
        vals = name.split("--")

        if use_detail == "false":
            node.name = vals[0]

        else:

            lsm = vals[1]
            gene_count = vals[2]

            rgb_tuple = cmap(norm(float(lsm)))
            rgb_hex = rgb2hex(rgb_tuple)
        
            node.name = vals[0] + " (L: %s) (G: %s)" % (lsm, gene_count)
            node.img_style["bgcolor"] = rgb_hex

use_detail = sys.argv[1]

if use_detail != "false":
    l_symbol = sys.argv[2]
    spec_ls = pd.read_csv('data/species_lifespans.txt', sep='\t')
    l_min, l_max = np.min(spec_ls[l_symbol].apply(np.min)), np.max(spec_ls[l_symbol].apply(np.max))
    norm = Normalize(float(l_min), float(l_max))

    cmap = plt.get_cmap('Blues')

tree = Tree("tree.nwk")
ts = TreeStyle()
ts.layout_fn = layout
tree.render("tree.png", w = 800, tree_style = ts)
