# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 13:09:01 2021

@author: kmihajlo
"""

import os, pickle
from scipy.stats import hypergeom
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from shutil import copyfile

legend = {'Control_D06':'C6','Control_D15':'C15', 'Control_D21':'C21', 'Control_IPSCs':'C0',
          'PINK1_D06':'PD6', 'PINK1_D15':'PD15', 'PINK1_D21':'PD21', 'PINK1_IPSCs':'PD0'}
nets_1 = ['PPI','GI','MI','COEX']
nets_2 = ['PPI+GI','PPI+MI','PPI+COEX','COEX+GI','COEX+MI','GI+MI']
nets_3 = ['GI+COEX+MI','PPI+GI+MI','PPI+COEX+MI','PPI+GI+COEX']
nets_everything = ['EXP'] + nets_1 + nets_2 + nets_3 + ['ALL']

def Entrez2Symbols(gene_info): #dict with entrez as keys and symbol as values (key:value)
    Entrez2Sym = {}
    gene_info.readline()
    for line in gene_info.readlines():
        lspt=line.strip().split('\t')
        Symbol = lspt[2]

        Entrez2Sym[lspt[1]] = Symbol
    return(Entrez2Sym)

Entrez2Symbol_file = open('input/Homo_sapiens.gene_info', 'r')
Entrez2Sym = Entrez2Symbols(Entrez2Symbol_file)
Entrez2Symbol_file.close()

# Benjamini-Hochberg p-value correction for multiple hypothesis testing
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
    #return p


def plotPD_EnrichClst_kmeans(all_nets, cell_cond, clustmeth, perc_enr_cluster_all, var_clusts_all, percPDgenes_all, var_genes_all, legend = legend):
   
    N = len(all_nets)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.2       # the width of the bars
    save_dir = f'output/{cell_cond}'

    if not os.path.exists(save_dir):
        os.makedirs(save_dir) 

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.grid(color='black', axis='y', linestyle='-', linewidth=0.4, alpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    for i in range(N):
        ae1 = np.array(list(zip([var_clusts_all[i-1][0]], [var_clusts_all[i-1][1]]))).T
        rects1 = ax.bar(ind[i-1], perc_enr_cluster_all[i-1], yerr=ae1, width=width, color='limegreen')
        ae2 = np.array(list(zip([var_genes_all[i-1][0]], [var_genes_all[i-1][1]]))).T
        rects2 = ax.bar(ind[i-1]+width, percPDgenes_all[i-1], yerr=ae2, width=width, color='r')

    ax.set_ylabel('%', fontsize = 20, fontweight = 'bold')
    ax.set_xticks(ind)
    ax.set_xticklabels(all_nets,  fontsize = 18, rotation=90) 
    ax.tick_params(axis='y', which='major', labelsize=16)
    ax.legend( (rects1[0], rects2[0]), ('PD enrich clusts', 'PD genes'),  fontsize = 18)

    plt.ylim(0,60)
    plt.title(f'{legend[cell_cond]} - {clustmeth}',  fontsize = 30, pad = 24)



    plt.tight_layout()    
    plt.savefig(f'{save_dir}/{clustmeth}.jpg', dpi = 350, format='jpg') 
    plt.show()  
    plt.close()    

def plotPD_EnrichClst(cell_conds, clust_meth, perc_enr_cluster_all, percPDgenes_all, legend = legend, rotation = 0):
    N = len(cell_conds)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.2       # the width of the bars
    save_dir = 'output'

    if not os.path.exists(save_dir):
        os.makedirs(save_dir) 

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.grid(color='black', axis='y', linestyle='-', linewidth=0.4, alpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    for i in range(N):
        print(i)
        rects1 = ax.bar(ind[i], perc_enr_cluster_all[i], width=width, color='limegreen')
        rects2 = ax.bar(ind[i]+width, percPDgenes_all[i], width=width, color='r')

    ax.set_ylabel('%', fontsize = 26, fontweight = 'bold')
    ax.set_xticks(ind+width/2)
    ax.set_xticklabels(cell_conds,  fontsize = 22, rotation = rotation) #nets_everythinga should be all_nets

    ax.tick_params(axis='y', which='major', labelsize=24)
    ax.legend((rects1[0], rects2[0]), ('PD enriched clusters', 'PD genes'),  fontsize = 24, loc = 1)

    plt.ylim(0,40)
    plt.title('Enrichment in PD genes',  fontsize = 30, pad = 24)

    plt.tight_layout()    
    plt.savefig(f'{save_dir}/{clust_meth}.jpg', dpi = 350, format='jpg')    
    plt.show()  
    plt.close()    
    


def PDenrichedClusts(G1_clust, PD_genes, Entrez2Sym = Entrez2Sym):
    genes = [str(item) for sublist in G1_clust for item in sublist]
    # genes = [Entrez2Sym[gene] for gene in genes]                    
    Possible_PDgenes = list(set(genes) & set(PD_genes))
                          
    Enriched_clusts = []
    PDpreds_clusts = []
    pvals = []
    Enriched_clusts_i = []
    PDgenes_cc = []
    fold_clust = []
    for i in range(len(G1_clust)): 
        genes_clust = G1_clust[i]
        # genes_clust = [Entrez2Sym[str(x)] for x in genes_clust]
        PDgenes_clust = [x for x in Possible_PDgenes if x in genes_clust]                            
        
        M = len(genes)
        K = len(Possible_PDgenes)
        N = len(genes_clust)
        X = len(PDgenes_clust)
        try:
            fold = (X/N)/(K/M)
        except ZeroDivisionError:
            fold = 0
        if fold >= 1:
            # print(fold)
            pval = hypergeom.sf(X-1, M, K, N)
            if pval <= 0.05: 
                Enriched_clusts_i.append(i)
                Enriched_clusts.append(genes_clust)
                fold_clust.append(fold)
                PDpreds = [x for x in genes_clust if x not in Possible_PDgenes]
                Clst_PDgenes = [x for x in genes_clust if x in Possible_PDgenes]
                PDgenes_cc.append(Clst_PDgenes)
                
                PDpreds_clusts.append(PDpreds)
                pvals.append(pval)

    # Benjamini-Hochberg p-value correction   
    pvals_adj =  p_adjust_bh(pvals)
    indexes = []
    for i in range(len(pvals_adj)):
        if pvals_adj[i] > 0.05:
            indexes.append(i)
    
    if len(indexes) >= 1:
        for index in sorted(indexes, reverse=True):
            del PDpreds_clusts[index]
            del Enriched_clusts[index]
            del Enriched_clusts_i[index]
            del PDgenes_cc[index]
            del fold_clust[index]
                
    # perc of clusters enriched in PD genes and percentage of total PD genes these clusters capture and percentage of preds, compared to all other genes in the condition
    enr_cluster = len(Enriched_clusts)
    total_cluster = len(G1_clust)
    perc_enr_cluster = 100.*enr_cluster/total_cluster    
           
    PDgenes_cc = [item for sublist in PDgenes_cc for item in sublist]
    percPDgenes = len(PDgenes_cc)/len(Possible_PDgenes)*100
    
    fold_clust_avg = np.mean(fold_clust)
    return perc_enr_cluster, percPDgenes, Enriched_clusts, Enriched_clusts_i, fold_clust_avg

def variation_whiskers(variable, mean_val):
    q3, q1 = np.percentile(variable, [84.13 ,15.87])
    min_error = mean_val - q1
    max_error = q3 - mean_val   
    variations = [min_error, max_error]
    return variations



in_dir = 'input/Clusters' 
out_dir = 'output'
with open('input/PD_genes_DGN.pkl', 'rb') as handle:
    PD_genes = pickle.load(handle)

for root, dirs, files in os.walk(in_dir):
    for file in files:
        print(root, file)
        cell_cond = root.split('/')[2]
        nets = root.split('/')[3]
    
        save_dir = f'{out_dir}/{cell_cond}/{nets}'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir) 

        clust_n = file.split('_')[3].split('.')[0]   
        with open(f'{root}/{file}', 'rb') as handle:
            G1_clust = pickle.load(handle)
   
          
        perc_enr_cluster, percPDgenes, Enriched_clusts, Enriched_clusts_i, _ = PDenrichedClusts(G1_clust, PD_genes)

        save_file1 = f'{save_dir}/kMeans_{clust_n}_PDEnrichClusts.pkl'
        save_file2 = f'{save_dir}/kMeans_{clust_n}_PDEnrichClusts_indices.pkl'
                               
        with open(save_file1, 'wb') as handle:
            pickle.dump(Enriched_clusts, handle)   
        with open(save_file2, 'wb') as handle:
            pickle.dump(Enriched_clusts_i, handle)                  
                    

#PART2 - Best run
   
cell_conds_ord = ['Control_IPSCs',
 'Control_D06',
 'Control_D15',
 'Control_D21',
 'PINK1_IPSCs',
 'PINK1_D06',
 'PINK1_D15',
 'PINK1_D21']

runs = 10

nclusts = {}
for root, dirs, files in os.walk('input/Clusters'):
    for file in files:
        if 'kMeans' in file:
            cell_cond = root.split('/')[2]
            cls_n = file.split('_')[2].split('c')[0]
            nclusts[cell_cond] = cls_n
print(nclusts)            
if not os.path.exists('output/best_runs'):
    os.makedirs('output/best_runs')

for nets in nets_everything:
    for cell_cond in cell_conds_ord:  
        print(cell_cond)
        percPDgenes_cc = []
        perc_enr_cluster_cc = []
        
        for run in range(runs):
            print(run)            
            with open(f'input/Clusters/{cell_cond}/{nets}/kMeans_G1_{nclusts[cell_cond]}cls_{run}.pkl', 'rb') as handle:
                G1_clust = pickle.load(handle)                
            perc_enr_cluster, percPDgenes,_ ,_,_ = PDenrichedClusts(G1_clust, PD_genes)
                           
            percPDgenes_cc.append(percPDgenes)
            perc_enr_cluster_cc.append(perc_enr_cluster)
               
        tmp = max(percPDgenes_cc)
        index_max = percPDgenes_cc.index(tmp)
        print(tmp, index_max, perc_enr_cluster_cc[index_max])
    
        save_dir = f'{out_dir}/{cell_cond}/{nets}'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        clst_file = f'kMeans_G1_{nclusts[cell_cond]}cls_{index_max}.pkl'
        
        with open(f'input/Clusters/{cell_cond}/{nets}/{clst_file}', 'rb') as handle:
            G1_clust = pickle.load(handle)
        perc_enr_cluster, percPDgenes, Enriched_clusts, Enriched_clusts_i,_ = PDenrichedClusts(G1_clust, PD_genes)
    
    
        save_file1 = f'{save_dir}/kMeans_{index_max}_PDEnrichClusts.pkl'
        save_file2 = f'{save_dir}/kMeans_{index_max}_PDEnrichClusts_indices.pkl'
    
        with open(save_file1, 'wb') as handle:
            pickle.dump(Enriched_clusts, handle)   
        with open(save_file2, 'wb') as handle:
            pickle.dump(Enriched_clusts_i, handle)  
    
        save_dir = f'output/best_runs/{cell_cond}/{nets}'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        copyfile(save_file1, f'{save_dir}/kMeans_{index_max}_PDEnrichClusts.pkl')
        copyfile(f'input/Clusters/{cell_cond}/{nets}/{clst_file}', f'{save_dir}/{clst_file}')
            
            
            
### plotting
percPDgenes_1meth = []
perc_enr_cluster_1meth = []
fold_clust_avg = []
for cell_cond in cell_conds_ord: 
    direct = f'output/best_runs/{cell_cond}/ALL/'   
    for file in os.listdir(direct):
        if 'kMeans_G1_' in file:
            file_n = os.path.join(direct, file)
                    
    with open(file_n, 'rb') as handle:
        G1_clust = pickle.load(handle)
    perc_enr_cluster, percPDgenes,_ ,_,fold_clust = PDenrichedClusts(G1_clust, PD_genes)
    fold_clust_avg.append(fold_clust)
    percPDgenes_1meth.append(percPDgenes)
    perc_enr_cluster_1meth.append(perc_enr_cluster)
cell_conds_ord_s = [legend[x] for x in cell_conds_ord]
plotPD_EnrichClst(cell_conds_ord_s, 'kMeans', perc_enr_cluster_1meth, percPDgenes_1meth)            


