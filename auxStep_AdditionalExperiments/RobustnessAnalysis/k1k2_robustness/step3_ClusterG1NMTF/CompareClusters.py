# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 17:10:24 2024

@author: Katarina
"""
import pickle
from sklearn.metrics import rand_score
import os, sys
import pandas as pd
import numpy as np
import itertools

def find_list_containing_elements(list1, list2):
    """
    Finds in which list each element of list1 is present in list2.

    Args:
    list1: The list of elements to search for.
    list2: The list of lists where the elements are searched.

    Returns:
    A list where each element represents the index of the list in list2 
    where the corresponding element from list1 is found.
    """
    indices_list = []
    for elem in list1:
        for idx, sublist in enumerate(list2):
            if elem in sublist:
                indices_list.append(idx)
                break  # Break once the element is found in one list
    return indices_list

def list_files_in_subdirectories(directory):
    file_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_list.append(os.path.join(root, file))
    return file_list

labels_key2 = {'C0':'Control_IPSCs', 'C6':'Control_D06', 'C15':'Control_D15', 'C21':'Control_D21', 
              'PD0':'PINK1_IPSCs', 'PD6':'PINK1_D06', 'PD15':'PINK1_D15', 'PD21':'PINK1_D21'}
#'C0':'Control_IPSCs',
cc = str(sys.argv[1])

print(cc)
genes = pd.read_csv(f'input/Geneslist/Geneslist_{cc}.csv', header=0, index_col=0, sep='\t')    
genes = list(genes.index)
runs = 10
dims = len(os.listdir(f'output/{cc}'))

cc_RS = np.empty((runs*dims, runs*dims))     

files = list_files_in_subdirectories(f'output/{cc}')

for i in range(len(files)):
    print(files[i])
    with open(files[i], 'rb') as handle:
        G1_clust1 = pickle.load(handle)  
    for j in range(i, len(files)):
        print(files[j])
        with open(files[j], 'rb') as handle:
            G1_clust2 = pickle.load(handle)  
        belonging1 = find_list_containing_elements(genes, G1_clust1)
        belonging2 = find_list_containing_elements(genes, G1_clust2)
        
        RS = rand_score(belonging1,belonging2)
        cc_RS[i][j] = RS
        cc_RS[j][i] = RS


clsts_dim = os.listdir(f'output/{cc}')
labels =  list(itertools.chain.from_iterable(itertools.repeat(x, runs) for x in clsts_dim))
cc_RS_df = pd.DataFrame(cc_RS)#,columns=labels,index=labels)


sd = 'output/RandScore'
if not os.path.exists(sd):
    os.makedirs(sd)
cc_RS_df = pd.DataFrame(cc_RS,columns=labels,index=labels)
cc_RS_df.to_csv(f'{sd}/{cc}_RandScore.csv')


ccs_RI = []           
for cc in labels_key2.keys():
    print(cc)
    cc_RS_df = pd.read_csv(f'{sd}/{labels_key2[cc]}_RandScore.csv',index_col=0, header=0)
    cc_RS_df = cc_RS_df.rename(columns=dict(zip(cc_RS_df.columns, labels)))
      
    
    # median, stdv of RS
    upper_triangle_values = np.triu(cc_RS_df.values, k=1)
    # Convert upper triangle values to a 1D array
    upper_triangle_values = upper_triangle_values[np.triu_indices(cc_RS_df.values.shape[0])]
    upper_triangle_values = upper_triangle_values[upper_triangle_values != 0]
    ccs_RI.append(upper_triangle_values)
    
    median = np.median(upper_triangle_values)
    stdv = np.std(upper_triangle_values)
    print(median, stdv)

print('stats for RI across all CCs')        
ccs_RI =  [item for row in ccs_RI for item in row]

median = np.median(ccs_RI)
stdv = np.std(ccs_RI)
print(median, stdv)

















