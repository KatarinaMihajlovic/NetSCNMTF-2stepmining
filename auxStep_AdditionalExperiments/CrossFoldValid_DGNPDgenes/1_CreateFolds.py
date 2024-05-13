# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 15:29:29 2024

@author: Katarina
"""

import pickle, random
from itertools import combinations

def divide_chunks(l, n):      
    # looping till length l 
    for i in range(0, len(l), n):  
        yield l[i:i + n]
        
with open('input/PD_genes_DGN.pkl', 'rb') as handle:
    PD_genes = pickle.load(handle) 


random.shuffle(PD_genes)

   
n = 276
fold5_DGN = list(divide_chunks(PD_genes, n)) 

Training_sets = []
Test_sets = []

Training = combinations(fold5_DGN, 4)
for combo in Training:
    Training_sets.append(list(set().union(*combo)))  

for TS in Training_sets:
    Test_sets.append(list(set(PD_genes)-set(TS)))

Training_Test = [[Training_sets[i],Test_sets[i]] for i in range(len(Test_sets))]

with open('output/Training_Test5FOLD_DGN.pkl', 'wb') as handle:
    pickle.dump(Training_Test, handle)

