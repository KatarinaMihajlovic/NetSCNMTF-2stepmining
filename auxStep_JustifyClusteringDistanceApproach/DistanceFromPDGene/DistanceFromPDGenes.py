# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 11:58:30 2022

@author: kmihajlo
"""

import pandas as pd
import os, pickle, math, sys
import numpy as np
from math import sqrt
from scipy.stats import hypergeom, mannwhitneyu
from scipy.spatial import distance
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def euclidean_distance(a, b):
    return sqrt(sum((e1-e2)**2 for e1, e2 in zip(a,b)))

def Entrez2Symbols(gene_info): #dict with entrez as keys and symbol as values (key:value)
    Entrez2Sym = {}
    gene_info.readline()
    for line in gene_info.readlines():
        lspt=line.strip().split('\t')
        Symbol = lspt[2]

        Entrez2Sym[int(lspt[1])] = Symbol
    return(Entrez2Sym)


def LitEnrich(All_litValid, gene_list, outf):  
    M = 0
    K = 0
    N = 0
    X = 0
    for gene in All_litValid:
        M+=1
        if gene in gene_list:
            N+=1
        if All_litValid[gene][0] > 0 or All_litValid[gene][1] != '':
            K+=1
            if gene in gene_list:
                X+=1
    
    print(M,K,N,X)
    outf.write(f'{M}\t{K}\t{N}\t{X}\n')
    
    perc = X/N*100
    print(perc)
    outf.write(f'{perc}%\n')
    
    try:
        fold = (X/N)/(K/M)
        print(fold)
        outf.write(f'{fold}\n')
    
    except ZeroDivisionError:
        fold = 0
    
    if fold >= 1:
        pval = hypergeom.sf(X-1, M, K, N)
        print(pval)
        outf.write(f'{pval}\n')
    
        if pval <= 0.05: 
            #print(f'Enriched clust {i},  fold = {fold},  pval = {pval}')
            #print(N, X)
            print('Enriched in LitValid - PubMed,Gene4PD')    
            outf.write('Enriched in LitValid - PubMed,Gene4PD\n')    
    outf.write('\n')
    print('\n')   

def logPubmedn(cc1cc2_PDpreds_LitValid):
    cc1cc2_pubmed = []
    for gene in cc1cc2_PDpreds_LitValid.keys():
        logn = math.log(cc1cc2_PDpreds_LitValid[gene][0] + 1, 10)
        cc1cc2_pubmed.append(logn)  
    return cc1cc2_pubmed

def PlotHist_MWU(cc1cc2_pubmed, cc1cc2_AllGenes_pubmed, cc_pair, save_dir):
    ##### MWU test - shows that clustered genes have higher number of citations 
    x = cc1cc2_pubmed
    y = cc1cc2_AllGenes_pubmed
    statmwu,pvalmwu = mannwhitneyu(x,y, alternative='greater')
    # print(pvalmwu)
    
    ##computing the histograms
    num_bin = 50  
    lst = cc1cc2_pubmed + cc1cc2_AllGenes_pubmed
    minv = min(lst)
    maxv = max(lst)
        
    bin_lims = np.linspace(minv,maxv,num_bin+1)
    bin_centers = 0.5*(bin_lims[:-1]+bin_lims[1:])
    bin_widths = bin_lims[1:]-bin_lims[:-1]
    
    hist1, _ = np.histogram(cc1cc2_AllGenes_pubmed, bins=bin_lims)
    hist2, _ = np.histogram(cc1cc2_pubmed, bins=bin_lims)
    
    ##normalizing
    hist1b = hist1/np.max(hist1)
    hist2b = hist2/np.max(hist2)
    
    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_title(f'{cc_pair}: PD preds vs Other genes', fontsize=32, pad=20)    
    #ax.set_title('Core PD preds vs Other Shared genes', fontsize=22, pad=20)    
    
    ax.bar(bin_centers, hist2b, width = bin_widths, align = 'center', alpha = 0.5, color='red')
    ax.bar(bin_centers, hist1b, width = bin_widths, align = 'center', alpha = 0.5, color='cornflowerblue')
    
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    ax.set_xlabel('log(num_PubMed + 1)', fontsize=28)
    ax.legend(handles= [Patch(facecolor='red', edgecolor='red',label=f'(1) PD preds ({len(x)})'),
                        Patch(facecolor='cornflowerblue', edgecolor='cornflowerblue',label=f'(2) Other Genes ({len(y)})')],
              loc='best', labelspacing=1, prop={'size': 24})
    pvalmwu_s = "{:.3e}".format(pvalmwu)
    if pvalmwu <= 0.05:
        fig.text(0.56, 0.63, f'MWU (1) > (2)\np-value = {pvalmwu_s}*', fontsize=24)
    else:
        fig.text(0.56, 0.63, f'MWU (1) > (2)\np-value = {pvalmwu_s}', fontsize=24)
    
    plt.savefig(f'{save_dir}/{cc_pair}_PDpredsvsOGs', dpi = 600)  
    plt.show()
    plt.close()    
 
    
### MAIN CODE
    
Entrez2Symbol_file = open('input/Homo_sapiens.gene_info', 'r')
Entrez2Sym = Entrez2Symbols(Entrez2Symbol_file)

in_dir = 'input'   
day_dict = {'0':'IPSCs', '6':'D06', '15':'D15', '21':'D21'}


base_dir = 'input/NMTF_G1s'
cell_conds = os.listdir(base_dir)

with open(f'{in_dir}/PD_genes_DGN.pkl', 'rb') as handle:
    PD_genes = pickle.load(handle)
    
Flag = bool(sys.argv[0])
if Flag == True:    
    DistanceToPDGenes_dict = {key:None for key  in cell_conds if 'PINK1' in key}    
    for CC in range(len(cell_conds)): 
        DistanceToPDGenes = [] 
        
        cell_cond = cell_conds[CC]
        if 'PINK1' in cell_cond:
            print(cell_cond)
            stage = cell_cond.split('_')[1]
            G1 = pd.read_csv(f'{base_dir}/{cell_cond}/ALL_G1_with_headers.csv', header=0, index_col=0, delimiter='\t')
            G1_Ctrl = pd.read_csv(f'{base_dir}/Control_{stage}/ALL_G1_with_headers.csv', header=0, index_col=0, delimiter='\t')
            
            G1 = G1.rename(index=Entrez2Sym)
            G1_Ctrl = G1_Ctrl.rename(index=Entrez2Sym)
            G1_genes = list(G1.index)
            G1_CtrlGenes = list(G1_Ctrl.index)
            
            PD_genes_cc = list(set(PD_genes) & set(G1_genes))
            G1_PDgenes = G1.loc[PD_genes_cc]
            G1_OtherGenes = G1.drop(PD_genes_cc)
            G1_genes = list(G1_OtherGenes.index)
            G1_PDCtrlGenes = [x for x in G1_CtrlGenes if x in G1_genes]
            G1_OtherGenes = G1_OtherGenes.loc[G1_PDCtrlGenes]
            
            cnt = 0
            for gene in G1_OtherGenes.index.to_list():
                cnt+=1
                if cnt % 1000 == 0:
                    print(cnt)
                min_dist = 1000
                PD_gene = None
                gene_embed = G1_OtherGenes.loc[gene].to_numpy()        
                for genePD in PD_genes_cc:
                    genePD_embed = G1_PDgenes.loc[genePD].to_numpy()
                    # cos_dist = distance.cosine(gene_embed, genePD_embed)
                    Euc_dist = distance.euclidean(gene_embed, genePD_embed)

                    if Euc_dist < min_dist:
                        min_dist = Euc_dist
                        PD_gene = genePD
                DistanceToPDGenes.append([gene, min_dist, PD_gene])
            DistanceToPDGenes.sort(key=lambda x: x[1])
        
            DistanceToPDGenes_dict[cell_cond] = DistanceToPDGenes

    with open('output/DistanceToPDGenes.pkl', 'wb') as fp:   
        pickle.dump(DistanceToPDGenes_dict, fp)
else:
    with open('output/DistanceToPDGenes.pkl', 'rb') as handle:
        DistanceToPDGenes_dict = pickle.load(handle)    



### Stage Specific PD predictions comparisson

stage_dict = {'0':'IPSCs', '6':'D06', '15':'D15', '21':'D21'}


top_n_genes_cases = ['StageSpec', '5_percentile', '10_percentile']
# DistanceToPDGenes_dict_ccs = {}

for topn_genes in top_n_genes_cases:   
    with open('input/LitValid_AllGenes.pkl', 'rb') as handle:
        LitValid_AllGenes = pickle.load(handle)   
    with open('input/All_CommonGenes_LitValid.pkl', 'rb') as handle:
        All_CommonGenes_LitValid = pickle.load(handle)  
        
    New_SSPDpreds_dict = {key:[] for key in DistanceToPDGenes_dict.keys()}    
    save_dir = f'output/{topn_genes}'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir) 
   
    file_out = 'StageSpecValids.txt'
    outf = open(f'{save_dir}/{file_out}', 'w')
    for stage in stage_dict:
        cc = f'PINK1_{stage_dict[stage]}'
        print(cc)
        DistanceToPDGenes_cc = DistanceToPDGenes_dict[cc]
        
        with open(f'input/StageSpecPreds/C{stage}PD{stage}_PDpreds.pkl', 'rb') as handle:
            StageSpecPreds = pickle.load(handle)         
        distances = [pair[1] for pair in DistanceToPDGenes_cc]
        
        if 'percentile' in topn_genes:
            percentile = int(topn_genes.split('_')[0])
            percentile_x = np.percentile(distances, percentile) 
            genes_topxpercentile = [pair[0] for pair in DistanceToPDGenes_cc if pair[1] <= percentile_x]
            topPDgenes = genes_topxpercentile
        
        else:
            potential_preds = [pair[0] for pair in DistanceToPDGenes_cc]
            topPDgenes = potential_preds[:len(StageSpecPreds)] #for len of predictions the same a StageSpec
            # topPDgenes = [x[0] for x in topPDgenes]  
    
        New_SSPDpreds_dict[cc] = topPDgenes           
        #litvalid
        outf.write(f'{cc}\n')
        OtherGenes = [pair[0] for pair in DistanceToPDGenes_cc if pair[0]] # including my top genes
        cc_dict = {your_key: LitValid_AllGenes[your_key] for your_key in OtherGenes}
        LitEnrich(cc_dict, topPDgenes, outf)
    
        # pubmed
        topPDgenes_dict = {your_key: LitValid_AllGenes[your_key] for your_key in topPDgenes}
        OtherGenes = [pair[0] for pair in DistanceToPDGenes_cc if pair[0] not in topPDgenes]
        cc_dict = {your_key: LitValid_AllGenes[your_key] for your_key in OtherGenes}
        
        topPDgenes_pubmed = logPubmedn(topPDgenes_dict)
        OtherGenes_pubmed = logPubmedn(cc_dict) #all genes a OtherGenes_LitValid
        PlotHist_MWU(topPDgenes_pubmed, OtherGenes_pubmed, f'PD{stage}', save_dir)   

    outf.close()                     
