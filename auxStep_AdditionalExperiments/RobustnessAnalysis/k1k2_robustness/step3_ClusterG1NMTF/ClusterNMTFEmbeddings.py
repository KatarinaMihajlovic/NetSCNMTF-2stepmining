# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 11:02:34 2021

@author: kmihajlo
"""


import pickle, os
from sklearn.cluster import KMeans
import pandas as pd
import numpy as np
import sys

def kMeans(M, nodes, n_clusters):
    M = M.tolist()
    nodes2coordinates = dict(zip(nodes, M))  
    Cluster_belonging = {k: [] for k in range(n_clusters)}
    
    kMeans = KMeans(n_clusters=n_clusters, init = 'random').fit(M)
    KMeans_labels = list(kMeans.labels_)
    
    for cluster_index in range(len(KMeans_labels)):
        cluster = KMeans_labels[cluster_index]
        node_coords = M[cluster_index]
        for node, coordinates in nodes2coordinates.items():
            if node_coords == coordinates:
                Cluster_belonging[cluster].append(node)
    
    #print(Cluster_belonging)
    Cluster_belonging_list = []
    for _, values in Cluster_belonging.items():
        Cluster_belonging_list.append(values)
    #print(Cluster_belonging_list)
    return Cluster_belonging_list


in_dir = 'input/NMTF_G1' 
out_dir = 'output'
km_runs = 10

cc = str(sys.argv[1])
k1k2 = str(sys.argv[2])

G1 = np.load(f'{in_dir}/{cc}/{k1k2}_ALL_G1.npy')
numClst = G1.shape[1]   
genes = pd.read_csv(f'input/Geneslist/Geneslist_{cc}.csv', header=0, index_col=0, sep='\t')
genes = list(genes.index)
  
save_dir = f'{out_dir}/{cc}/{k1k2}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir) 
for kmrun in range(km_runs):
    print(kmrun)
    clustersG1_kmean = kMeans(G1, genes, numClst)
    #print(clustersG1_kmean)
    with open(f'{save_dir}/kMeans_G1_{numClst}cls_{kmrun}.pkl', 'wb') as handle:
        pickle.dump(clustersG1_kmean, handle)




    
    