# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 15:16:29 2024

@author: Katarina
"""

from itertools import combinations
import pandas as pd
import os, random, pickle

def divide_chunks(l, n):      
    # looping till length l 
    for i in range(0, len(l), n):  
        yield l[i:i + n]

legend = {'Control_D06':'C6','Control_D15':'C15', 'Control_D21':'C21', 'Control_IPSCs':'C0',
          'PINK1_D06':'PD6', 'PINK1_D15':'PD15', 'PINK1_D21':'PD21', 'PINK1_IPSCs':'PD0'}

folds = 5
in_dir= 'input/Geneslist'
for file in os.listdir(in_dir):
    cc = file.split('_')[1] + '_' + file.split('_')[2].split('.')[0]
    cc = legend[cc]
    genelist = pd.read_csv(f'{in_dir}/{file}', header=0, index_col=0)
    genelist = list(genelist.index)
    random.shuffle(genelist)
    n = int(len(genelist)/folds)
    
    fold5_sets = list(divide_chunks(genelist, n+1)) 
    Training_sets = []
    Test_sets = []
    
    Training = combinations(fold5_sets, 4)
    for combo in Training:
        Training_sets.append(list(set().union(*combo)))  
    
    for TS in Training_sets:
        Test_sets.append(list(set(genelist)-set(TS)))
    
    Training_Test = [[Training_sets[i],Test_sets[i]] for i in range(len(Test_sets))]
    
    with open(f'output/{cc}_Train_Test5FOLD.pkl', 'wb') as handle:
        pickle.dump(Training_Test, handle)
        
