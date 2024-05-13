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
import seaborn as sns
from statannotations.Annotator import Annotator

legend = {'Control_D06':'C6','Control_D15':'C15', 'Control_D21':'C21', 'Control_IPSCs':'C0',
          'PINK1_D06':'PD6', 'PINK1_D15':'PD15', 'PINK1_D21':'PD21', 'PINK1_IPSCs':'PD0'}
nets_1 = ['PPI','GI','MI','COEX']
nets_2 = ['PPI+GI','PPI+MI','PPI+COEX','COEX+GI','COEX+MI','GI+MI']
nets_3 = ['GI+COEX+MI','PPI+GI+MI','PPI+COEX+MI','PPI+GI+COEX']
nets_everything = ['EXP'] + nets_1 + nets_2 + nets_3 + ['ALL']


# Benjamini-Hochberg p-value correction for multiple hypothesis testing
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
    #return p


def PDenrichedClusts(G1_clust, PD_genes):
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



in_dir = 'input/Clusters' 
out_dir = 'output'
with open('output/Training_Test5FOLD_DGN.pkl', 'rb') as handle:
    Training_Test5FOLD_DGN = pickle.load(handle)

ccs = os.listdir(in_dir)
nets = os.listdir(f'{in_dir}/{ccs[0]}')
nets = ['ALL','PPI+GI+COEX']

for net in nets:
    print(net)
    net_testenr, net_testenr_avg = [],[]
    net_test_notenr, net_test_notenr_avg = [],[]
    
    net_foldTest, netFold_BG = [],[]
    for cc in ccs:
        for file in os.listdir(f'{in_dir}/{cc}/{net}'):
            
            with open(f'{in_dir}/{cc}/{net}/{file}', 'rb') as handle:
                G1_clust = pickle.load(handle)    
        
        
            save_dir = f'{out_dir}/{net}/{cc}'
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)         
            clust_n = file.split('_')[3].split('.')[0]   
            
            Validation_run = []
            for i in range(len(Training_Test5FOLD_DGN)):
                Train = Training_Test5FOLD_DGN[i][0]
                Test = Training_Test5FOLD_DGN[i][1]
                Total_genes = [item for sublist in G1_clust for item in sublist]
                BG = list(set(Total_genes)-set(Train)-set(Test))
                              
                perc_enr_cluster, percPDgenes, Enriched_clusts, Enriched_clusts_i, _ = PDenrichedClusts(G1_clust, Train)
                Enriched_clusts_size = len([item for sublist in Enriched_clusts for item in sublist])
                
                Test_enr = []
                BG_enr = []
                for eclst in Enriched_clusts:
                    test_genes = list(set(eclst)&set(Test))
                    bg_genes = list(set(eclst)-set(Test))
                    Test_enr.append(test_genes)
                    BG_enr.append(bg_genes)
                Test_enr = list(set().union(*Test_enr))
                # Test_fold = len(Test_enr)/len(Test)
                BG_enr = list(set().union(*BG_enr))
                # BG_fold = len(BG_enr)/len(BG)

                Test_fold = (len(Test_enr)/Enriched_clusts_size)/(len(Test)/len(Total_genes))
                BG_fold = (len(BG_enr)/Enriched_clusts_size)/(len(BG)/len(Total_genes))

                net_foldTest.append(Test_fold)                
                netFold_BG.append(BG_fold)                
     
        

    net_foldTest_avg = sum(net_foldTest) / len(net_foldTest) 
    net_foldBG_all = sum(netFold_BG) / len(netFold_BG) 
    
    print(net_foldTest_avg/net_foldBG_all)
        
    foldTestBGgenes = net_foldTest + netFold_BG
    case = ['FoldTest']*len(net_foldTest) + ['FoldBG']*len(netFold_BG)
    avg_enrichments_df = pd.DataFrame({'value': foldTestBGgenes, 'case': case})
    
    
    fig, ax = plt.subplots(figsize=(12, 8))
    # plt.grid(axis='y', alpha = 0.5)
    ax = sns.boxplot(data=avg_enrichments_df, x='case', y='value', palette=sns.color_palette('pastel'))
    annotator = Annotator(ax, [['FoldTest','FoldBG']], data=avg_enrichments_df, x='case', y='value', order=['FoldTest','FoldBG'])
    annotator.configure(test='Mann-Whitney-gt', text_format='full',fontsize=22, loc='outside')
    annotator.apply_and_annotate()
    ax.set_ylabel('Fold', fontsize = 30)#, pad = 24)
    ax.set_xlabel('')
    
    ax.set_xticklabels(['test DisGeNet PD genes','Background'])
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=20)
    plt.savefig(f'{out_dir}/{net}/FoldDifferences.jpg',  dpi = 350, format='jpg')  

    plt.show()
    plt.close()

# do pd genes from test populate enriched clusters more than the background genes populate it         
