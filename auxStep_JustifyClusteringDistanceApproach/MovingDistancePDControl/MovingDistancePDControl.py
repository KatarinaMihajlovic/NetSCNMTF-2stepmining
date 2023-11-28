# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 11:19:19 2022

@author: kmihajlo
"""

# Euclidean distance between PD-Control stage-spec pairs, all genes expressed in them
# order according to highest distance
# intersection and rank according to avg highest distance between all

import pandas as pd
import os, pickle, math
import numpy as np
from scipy.stats import hypergeom, mannwhitneyu
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.spatial import distance


def Entrez2Symbols(gene_info): #dict with entrez as keys and symbol as values (key:value)
    Entrez2Sym = {}
    gene_info.readline()
    for line in gene_info.readlines():
        lspt=line.strip().split('\t')
        Symbol = lspt[2]

        Entrez2Sym[int(lspt[1])] = Symbol
    return(Entrez2Sym)

def gene_distances(EucDist_PD_vals, EucDist_C_vals, genes, day, case):
    dists = []
    for i in range(len(EucDist_PD_vals)):
        dist_PD = EucDist_PD_vals[i]
        dist_C = EucDist_C_vals[i]
        if case == 'Eucl_dist': 
            dist_gene = distance.euclidean(dist_PD, dist_C)
        elif case == 'Manh_dist': 
            dist_gene = distance.cityblock(dist_PD, dist_C)
        elif case == 'Euc_dist': 
            dist_gene = distance.cosine(dist_PD, dist_C)
        dists.append(dist_gene)
    
    dists_pd = pd.DataFrame(dists, index = genes, columns = [case])
    dists_pd = dists_pd.sort_values(by=case, ascending=False)
    
    return dists_pd

def logPubmedn(cc1cc2_PDpreds_LitValid):
    cc1cc2_pubmed = []
    for gene in cc1cc2_PDpreds_LitValid.keys():
        logn = math.log(cc1cc2_PDpreds_LitValid[gene][0] + 1, 10)
        cc1cc2_pubmed.append(logn)  
    return cc1cc2_pubmed

def PlotHist_MWU(cc1cc2_pubmed, cc1cc2_AllGenes_pubmed, cc_pair, save_dir, label2 = 'Other Genes'):
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
    ax.set_title(f'{cc_pair}: PD preds vs {label2}', fontsize=32, pad=20)    
    
    ax.bar(bin_centers, hist2b, width = bin_widths, align = 'center', alpha = 0.5, color='red')
    ax.bar(bin_centers, hist1b, width = bin_widths, align = 'center', alpha = 0.5, color='cornflowerblue')
    
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    ax.set_xlabel('log(num_PubMed + 1)', fontsize=28)
    ax.legend(handles= [Patch(facecolor='red', edgecolor='red',label=f'(1) PD preds ({len(x)})'),
                        Patch(facecolor='cornflowerblue', edgecolor='cornflowerblue',label=f'(2) {label2} ({len(y)})')],
              loc='best', labelspacing=1, prop={'size': 24})
    pvalmwu_s = "{:.3e}".format(pvalmwu)
    if pvalmwu <= 0.05:
        fig.text(0.56, 0.6, f'MWU (1) > (2)\np-value = {pvalmwu_s}*\n', fontsize=24)
    else:
        fig.text(0.56, 0.6, f'MWU (1) > (2)\np-value = {pvalmwu_s}\n', fontsize=24)
    
    if label2 == 'Other Genes':
        plt.savefig(f'{save_dir}/{cc_pair}_PDpredsvsOGs', dpi = 600)  
    elif label2 == 'Stage Spec PD preds':
        plt.savefig(f'{save_dir}/{cc_pair}_PDpredsvsSSPs', dpi = 600)  
    elif label2 == 'Core PD preds':
        plt.savefig(f'{save_dir}/{cc_pair}_PDpredsvsCPs', dpi = 600)  

    plt.show()
    plt.close()    
 

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
            print('Enriched in LitValid - PubMed,Gene4PD')    
            outf.write('Enriched in LitValid - PubMed,Gene4PD\n')
 
    outf.write('\n')
    print('\n')   

def topn_corepreds(DistanceToPDGenes_dict_ccs, PD_genes, LitValid_AllGenes, All_CommonGenes_LitValid, New_SSPDpreds_dict, save_dir, case):
    for gene in PD_genes:
        All_CommonGenes_LitValid.pop(gene, None)
    AllCommonGenes = [[pair[0] for pair in DistanceToPDGenes_dict_ccs[cc]] for cc in DistanceToPDGenes_dict_ccs]
    AllCommonGenes = set.intersection(*map(set,AllCommonGenes))   
    AllCommonGenes = [gene for gene in AllCommonGenes if gene not in PD_genes]
    All_CommonGenes_LitValid = {your_key: LitValid_AllGenes[your_key] for your_key in AllCommonGenes}

    New_SSPDpreds = New_SSPDpreds_dict.values()
    common_genes = list(set.intersection(*map(set,New_SSPDpreds)))     
    
    OtherGenes_allCCs_dist_pd = pd.DataFrame(index= common_genes, columns = DistanceToPDGenes_dict_ccs.keys())
    for cc in DistanceToPDGenes_dict_ccs:
        for i in range(len(DistanceToPDGenes_dict_ccs[cc])):
            gene = DistanceToPDGenes_dict_ccs[cc][i][0]
            if gene in common_genes:
                distance_cg = DistanceToPDGenes_dict_ccs[cc][i][1]
                OtherGenes_allCCs_dist_pd.loc[gene,cc] = distance_cg
    
    OtherGenes_allCCs_dist_pd['average'] = OtherGenes_allCCs_dist_pd.mean(axis=1)
    OtherGenes_allCCs_dist_pd = OtherGenes_allCCs_dist_pd.sort_values(by='average', ascending=True)
    common_genes_ordered = OtherGenes_allCCs_dist_pd.index.tolist()
    
    file_out = 'StagesIntersection_litValid_enrich.txt'
    outf = open(f'{save_dir}/{file_out}', 'w')
    print('Core Preds\n')
    LitEnrich(All_CommonGenes_LitValid, common_genes_ordered, outf)
    outf.close()   
            
    topnCGs_dict = {your_key: All_CommonGenes_LitValid[your_key] for your_key in common_genes_ordered}
    with open(f'{save_dir}/topnCGs_{case}.pkl', 'wb') as fp:   
        pickle.dump(topnCGs_dict, fp)
        
    with open(f'{save_dir}/StagesIntersection_LitValid.txt', 'w') as f:
        for gene in topnCGs_dict.keys():
            f.write(f'{gene}\t{topnCGs_dict[gene][0]}\t{topnCGs_dict[gene][1]}\n')
    with open(f'{save_dir}/StagesIntersection_MovDist.txt', 'w') as f:
        for gene in topnCGs_dict.keys():
            f.write(f'{gene}\n')
            
    OtherGenes = [gene for gene in All_CommonGenes_LitValid if gene not in common_genes_ordered]
    cc_dict = {your_key: All_CommonGenes_LitValid[your_key] for your_key in OtherGenes}
    
    with open('input/CorePreds_LitValid.pkl', 'rb') as handle:
        CorePreds_LitValid = pickle.load(handle) 
    CorePreds_LitValid_pubmed = logPubmedn(CorePreds_LitValid)
        
    topnCGs_pubmed = logPubmedn(topnCGs_dict)
    OtherGenes_pubmed = logPubmedn(cc_dict) #all genes a OtherGenes_LitValid
    PlotHist_MWU(topnCGs_pubmed, OtherGenes_pubmed, 'Stages Intersection', save_dir)
    PlotHist_MWU(topnCGs_pubmed, CorePreds_LitValid_pubmed, 'Stages Intersection', save_dir, label2 = 'Core PD preds')        


### MAIN CODE

Entrez2Symbol_file = open('input/Homo_sapiens.gene_info', 'r')
Entrez2Sym = Entrez2Symbols(Entrez2Symbol_file)
Entrez2Symbol_file.close()

stage_dict = {'0':'IPSCs', '6':'D06', '15':'D15', '21':'D21'}
day_l = {'IPSCs':'0', 'D06':'6', 'D15':'15', 'D21':'21'}
type_l = {'Control':'C', 'PINK1':'PD'}


in_dir = 'input'   


with open(f'{in_dir}/PD_genes_DGN.pkl', 'rb') as handle:
    PD_genes = pickle.load(handle)
    
with open('input/All_CommonGenes_LitValid.pkl', 'rb') as handle:
    All_CommonGenes_LitValid = pickle.load(handle)    
common_genes = list(All_CommonGenes_LitValid.keys())


base_dir = 'input/NMTF_G1s'
cell_conds = os.listdir(base_dir)

Euc_dist_common_dict = {key:None for key in cell_conds if 'PD' in key}         
Euc_dist_allccgenes_dict = {key:None for key in cell_conds if 'PD' in key}      

days = list(day_l.keys())


for day in days: 
    in_dir = 'input/SimMeasuresG'
    print(day)
    ccs2compare = f'C{day_l[day]}PD{day_l[day]}'

    EucDist_PD = np.load(f'{in_dir}/PINK1_{day}_PWEucDist.npy')
    with open(f'{in_dir}/PINK1_{day}_Genes.pkl', 'rb') as handle:
        Genes_PD = pickle.load(handle)  
    Genes_PD = list(Genes_PD)  
    EucDist_PD_df = pd.DataFrame(EucDist_PD, index = Genes_PD, columns = Genes_PD)
    
    EucDist_C= np.load(f'{in_dir}/Control_{day}_PWEucDist.npy')
    with open(f'{in_dir}/Control_{day}_Genes.pkl', 'rb') as handle:
        Genes_C = pickle.load(handle)  
    Genes_C = list(Genes_C)   
    EucDist_C_df = pd.DataFrame(EucDist_C, index = Genes_C, columns = Genes_C)

    Genes_cc1cc2 = list(set(Genes_C) & set(Genes_PD))      
    
    EucDist_PD_df = EucDist_PD_df[Genes_cc1cc2]
    EucDist_PD_df = EucDist_PD_df.loc[Genes_cc1cc2]
    EucDist_C_df = EucDist_C_df[Genes_cc1cc2]
    EucDist_C_df = EucDist_C_df.loc[Genes_cc1cc2]   
    
    EucDist_PD_vals = EucDist_PD_df.values
    EucDist_C_vals = EucDist_C_df.values

    Euc_dists_pd = gene_distances(EucDist_PD_vals, EucDist_C_vals, Genes_cc1cc2, day, case = 'Eucl_dist')
    
    Euc_dist_pairs = []
    Euc_dist_pairs_common = []
    for gene in Euc_dists_pd.index:
        if gene not in PD_genes:
            Euc_dist_pairs.append([gene, Euc_dists_pd.loc[gene,'Eucl_dist']])
            if gene in common_genes:
                Euc_dist_pairs_common.append([gene, Euc_dists_pd.loc[gene,'Eucl_dist']])

    Euc_dist_pairs.sort(key=lambda x: x[1] , reverse=True)  
    Euc_dist_pairs_common.sort(key=lambda x: x[1] , reverse=True)  

    Euc_dist_allccgenes_dict[f'PD{day_l[day]}'] = Euc_dist_pairs
    Euc_dist_common_dict[f'PD{day_l[day]}'] = Euc_dist_pairs_common

with open('input/LitValid_AllGenes.pkl', 'rb') as handle:
    LitValid_AllGenes = pickle.load(handle)   
                
########################## Stage Specific predictions

###### movement threshold accorind to the movement of each indiviudal stage

top_n_genes_cases = ['StageSpec', '95_percentile', '90_percentile']
file_out = 'StageSpecValids.txt'

for topn_genes in top_n_genes_cases:   
    print(topn_genes)
    New_SSPDpreds_dict = {key:[] for key in Euc_dist_allccgenes_dict.keys()}
    save_dir = f'output/{topn_genes}'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir) 

    outf = open(f'{save_dir}/{file_out}', 'w')
    for stage in stage_dict:
        # cc = f'PINK1_{stage_dict[stage]}'
        cc = f'PD{stage}'
        print(cc)
        Euc_dist_cc = Euc_dist_allccgenes_dict[cc]
        
        with open(f'input/StageSpecPreds/C{stage}PD{stage}_PDpreds.pkl', 'rb') as handle:
            StageSpecPreds = pickle.load(handle) 
        StageSpecPreds_valid = logPubmedn(StageSpecPreds)
        distances = [pair[1] for pair in Euc_dist_cc]
        
        if 'percentile' in topn_genes:
            percentile = int(topn_genes.split('_')[0])
            percentile_x = np.percentile(distances, percentile) 
            genes_topxpercentile = [pair[0] for pair in Euc_dist_cc if pair[1] >= percentile_x]
            topPDgenes = genes_topxpercentile
        
        else:
            potential_preds = [pair[0] for pair in Euc_dist_cc]
            topPDgenes = potential_preds[:len(StageSpecPreds)] #for len of predictions the same a StageSpec
        
        New_SSPDpreds_dict[cc] = topPDgenes           

        #litvalid
        outf.write(f'{cc}\n')
        OtherGenes = [pair[0] for pair in Euc_dist_cc if pair[0]] # including my top genes
        cc_dict = {your_key: LitValid_AllGenes[your_key] for your_key in OtherGenes}
        LitEnrich(cc_dict, topPDgenes, outf)
    
        # pubmed
        topPDgenes_dict = {your_key: LitValid_AllGenes[your_key] for your_key in topPDgenes}
        OtherGenes = [pair[0] for pair in Euc_dist_cc if pair[0] not in topPDgenes]
        cc_dict = {your_key: LitValid_AllGenes[your_key] for your_key in OtherGenes}
        
        topPDgenes_pubmed = logPubmedn(topPDgenes_dict)
        OtherGenes_pubmed = logPubmedn(cc_dict) #all genes a OtherGenes_LitValid
        PlotHist_MWU(topPDgenes_pubmed, OtherGenes_pubmed, f'PD{stage}', save_dir)   
        PlotHist_MWU(topPDgenes_pubmed, StageSpecPreds_valid, f'PD{stage}', save_dir, label2='Stage Spec PD preds')   

    outf.close()                     
    topn_corepreds(Euc_dist_allccgenes_dict, PD_genes, LitValid_AllGenes, All_CommonGenes_LitValid, New_SSPDpreds_dict, save_dir, topn_genes)

