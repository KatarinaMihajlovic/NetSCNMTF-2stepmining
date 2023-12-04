# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 11:02:34 2021

@author: kmihajlo
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 11:28:29 2021

@author: bscuser
"""
import pickle, os
from sklearn.cluster import KMeans
import pandas as pd


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


in_dir = 'input' 
out_dir = 'output'
km_runs = 10
   
for root, dirs, files in os.walk(in_dir):
    for file in files:
        print(file)
        cell_cond = root.split(os.sep)[1]
        nets = file.split('_')[0]
        G1_df = pd.read_csv(f'{root}/{file}', header=0, index_col=0, sep='\t')
        G1_NumClusters = len(list(G1_df.columns))
        # G1_NumClusters = 50
        genes = list(G1_df.index)
        G1 = G1_df.values
        # if 'ALL' in file:
        print(cell_cond, nets)
        save_dir = f'{out_dir}/{cell_cond}/{nets}'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)             
        for kmrun in range(km_runs):
            print(kmrun)
            clustersG1_kmean = kMeans(G1, genes, G1_NumClusters)
            #print(clustersG1_kmean)
            with open(f'{save_dir}/kMeans_G1_{G1_NumClusters}cls_{kmrun}.pkl', 'wb') as handle:
                pickle.dump(clustersG1_kmean, handle)

    
    
