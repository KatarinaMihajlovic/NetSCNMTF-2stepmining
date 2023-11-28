# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 12:38:35 2021

@author: kmihajlo
"""

import pandas as pd
import os, pickle
from sklearn.metrics.pairwise import euclidean_distances
import numpy as np


def write_txt(file_path, name, list_l):
    if not os.path.exists(file_path):
        os.makedirs(file_path) 
    with open(f'{file_path}/{name}', 'w') as f: 
        for x in list_l:
            f.write(f'{x}\n')
in_dir = 'input'   


base_dir = 'input/NMTF_G1s'
cell_conds = os.listdir(base_dir)

for C1_i in range(len(cell_conds)): 
    cell_cond_1 = cell_conds[C1_i]
    print(cell_cond_1)
    cc_1d = cell_cond_1.split('_')[1]
    G1 = pd.read_csv(f'{base_dir}/{cell_cond_1}/ALL_G1_with_headers.csv', header=0, index_col=0, delimiter='\t')
    S5 = pd.read_csv(f'{base_dir}/{cell_cond_1}/ALL_S_EXP.csv', header=0, index_col=0, delimiter='\t')

    G1_vals = G1.to_numpy()
    S5_vals = S5.to_numpy()
    G1_S5 = G1_vals.dot(S5_vals)

    # Normalized Euclidean distance
    EuclidDist_all = euclidean_distances(G1_S5)
    feature_norms = np.linalg.norm(EuclidDist_all) #, axis = 0)
    print(feature_norms)
    EuclidDist_all_normalized = EuclidDist_all/feature_norms
    print(EuclidDist_all_normalized)

    file_path = 'output/EuclideanDistance'
    if not os.path.exists(file_path):
        os.makedirs(file_path) 
    save_file = f'{file_path}/{cell_cond_1}_PWEucDist'
    np.save(save_file, EuclidDist_all_normalized)  
    
    with open(f'output/{cell_cond_1}_Genes.pkl', 'wb') as handel:
        pickle.dump(G1.index, handel)    
