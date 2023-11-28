# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 10:12:57 2021

@author: kmihajlo
"""
import os, math
import pickle
from scipy.stats import hypergeom, mannwhitneyu
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import pandas as pd
import venn
import csv


with open('input/PD_genes_DGN.pkl', 'rb') as handle:
    PD_genes = pickle.load(handle) 

def logPubmedn(cc1cc2_PDpreds_LitValid):
    cc1cc2_pubmed = []
    cc1cc2_genes_pubmed = []
    for gene in cc1cc2_PDpreds_LitValid.keys():
        logn = math.log(cc1cc2_PDpreds_LitValid[gene][0] + 1, 10)
        cc1cc2_pubmed.append(logn)  
        cc1cc2_genes_pubmed.append([gene,logn])
    return cc1cc2_pubmed, cc1cc2_genes_pubmed


def PlotHist_MWU(cc1cc2_pubmed, cc1cc2_AllGenes_pubmed, cc_pair, save_dir):
    ##### MWU test - shows that clustered genes have higher number of citations 
    xvals = cc1cc2_pubmed
    yvals = cc1cc2_AllGenes_pubmed
    statmwu,pvalmwu = mannwhitneyu(xvals,yvals, alternative='greater')
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
    ax.set_ylabel('', fontsize=28)

    ax.legend(handles= [Patch(facecolor='red', edgecolor='red',label=f'(1) PD preds ({len(xvals)})'),
                        Patch(facecolor='cornflowerblue', edgecolor='cornflowerblue',label=f'(2) Other Genes ({len(yvals)})')],
              loc='best', labelspacing=1, prop={'size': 24})
    pvalmwu_s = "{:.3e}".format(pvalmwu)
    if pvalmwu < 0.05:
        fig.text(0.57, 0.63, f'MWU (1) > (2)\np-value = {pvalmwu_s}*', fontsize=24)
    else:
        fig.text(0.57, 0.63, f'MWU (1) > (2)\np-value = {pvalmwu_s}', fontsize=24)

    
    plt.savefig(f'{save_dir}/{cc_pair}_PDpredsvsOGs.jpg',  dpi = 350, format='jpg')  
    plt.show()
    plt.close()    

    
def BarPlotValidGenes(PDpreds_enrich, file_out):
    ### plot DEGs enrichment of pairwise PD preds
    cc_pairs = [x[0] for x in PDpreds_enrich]
    PD_preds = [x[1] for x in PDpreds_enrich]
    ValidGenes = [x[2] for x in PDpreds_enrich]
    p_vals =  [x[3] for x in PDpreds_enrich]
    cc_pairs.insert(1, cc_pairs[-1])
    cc_pairs.pop(-1)
    PD_preds.insert(1, PD_preds[-1])
    PD_preds.pop(-1)
    ValidGenes.insert(1, ValidGenes[-1])
    ValidGenes.pop(-1)
    p_vals.insert(1, p_vals[-1])
    p_vals.pop(-1)
    
    cc_pairs = ['PD_0','PD_6','PD_15','PD_21']
    
    width = 0.3    
    opacity = 0.7
    
    fig, ax = plt.subplots(figsize=(16, 9))
    # Remove top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    # Remove y-axis tick marks
    ax.yaxis.set_ticks_position('none')
    # Set plot title
    
    ax.set_title('Enrichment of PD predictions in PD-related genes', fontsize = 28)
    # Add major gridlines in the y-axis
    ax.grid(color='black', axis='y', linestyle='-', linewidth=0.4, alpha=0.9)
    ax.set_facecolor('xkcd:white')
    
    bar1 = ax.bar(np.arange(len(PD_preds)), PD_preds, width, align='center', alpha=opacity, color='r', label='PD predictions')
    bar2 = ax.bar(np.arange(len(ValidGenes)) + width, ValidGenes, width, align='center', alpha=opacity,  color='limegreen', label='PD-related genes')
     
    ax.legend(fontsize=22, loc = 'upper left')
    ax.set_ylabel('#genes', fontsize = 24, fontweight = 'bold')
    ax.set_xticks(np.arange(len(PD_preds))+width/2)
    ax.set_xticklabels(cc_pairs,  fontsize = 24) 
    ax.tick_params(axis='y', which='major', labelsize=24)
    plt.ylim(0,1800)
    i = 0
    for rect in bar1:
        height = rect.get_height()
        if float(p_vals[i]) <= 0.05:
            pval = "{:.3e}".format(p_vals[i])
            plt.text(rect.get_x() + rect.get_width(), height+0.1, f'{pval}*', ha='center', va='bottom', fontsize = 24)
        else:
            plt.text(rect.get_x() + rect.get_width(), height+0.1, f'{p_vals[i]}', ha='center', va='bottom', fontsize = 24)
        i+=1
    
    plt.savefig('output/StageSpecPreds/LitValid_Enrich.jpg',  dpi = 350, format='jpg', bbox_inches='tight') 
    plt.show()
    plt.close()


def Entrez2Symbols(gene_info): #dict with entrez as keys and symbol as values (key:value)
    Entrez2Sym = {}
    gene_info.readline()
    for line in gene_info.readlines():
        lspt=line.strip().split('\t')
        Symbol = lspt[2]

        Entrez2Sym[lspt[1]] = Symbol
    return(Entrez2Sym)

def LitEnrich(All_litValid, gene_list, outf, case):  
    M = 0
    K = 0
    N = 0
    X = 0
    LitEnr_genes = []
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
            print(f'percentage of enriched genes: {X/N*100}')
            print(f'percentage of enriched background: {K/M*100}')
            outf.write(f'percentage of enriched genes: {X/N*100}\n')
            outf.write(f'percentage of enriched background: {K/M*100}\n')
            
    print('\n')  
    LitEnr_genes = [case, N, X, pval]
    return LitEnr_genes

colorp = {"Other_gene":'cornflowerblue', "PDpred_valid":'red'}
alphap = {"Other_gene":0.2, "PDpred_valid":0.8}
markerp = {"Other_gene":'o', "PDpred_valid":'x'}


def PDGs2Background_dist(Predictions_pubmed, BG_pubmed, day, save_dir):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)  
        
    statmwu,pvalmwu = mannwhitneyu(Predictions_pubmed, BG_pubmed, alternative='greater')

    ##computing the histograms
    num_bin = 50
    lst = list(Predictions_pubmed)
    minv = min(lst)
    maxv = max(lst)
        
    bin_lims = np.linspace(minv,maxv,num_bin+1)
    bin_centers = 0.5*(bin_lims[:-1]+bin_lims[1:])
    bin_widths = bin_lims[1:]-bin_lims[:-1]

    hist1, _ = np.histogram(Predictions_pubmed, bins=bin_lims)
    hist2, _ = np.histogram(BG_pubmed, bins=bin_lims)
    print(np.sum(hist1), np.sum(hist2))
    ##normalizing
    hist1b = hist1/np.sum(hist1)*100 # divide with sum, the distribution surface sums to 1
    hist2b = hist2/np.sum(hist2)*100

    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_facecolor('xkcd:white')
    if day != 'Core PD predictions':
        title = f'PD$_{{{day}}}$'
    else:
        title = day
        
    ax.set_title(title, fontsize=40, pad=20) #': {PD_Gs} vs Other genes'   

    ax.bar(bin_centers, hist1b, width = bin_widths, align = 'center', alpha = 0.5, color='red')
    ax.bar(bin_centers, hist2b, width = bin_widths, align = 'center', alpha = 0.5, color='cornflowerblue')
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    ax.xaxis.offsetText.set_fontsize(30)
    
    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)
    ax.set_xlabel('log(PubMed co-occurences + 1)', fontsize=30, fontweight = 'bold') #dist_measure_l[dist_meas]
    ax.set_ylabel('relative count (%)', fontsize=30, fontweight = 'bold')
    for label in ax.get_xaxis().get_ticklabels()[::2]:
        label.set_visible(False)
    
    pvalmwu_s = "{:.3e}".format(pvalmwu)
    if pvalmwu <= 0.05:
        fig.text(0.5, 0.52, f'MWU (1) > (2)\np-value = {pvalmwu_s}*\n', fontsize=28)
    else:
        fig.text(0.5, 0.52, f'MWU (1) > (2)\np-value = {pvalmwu_s}\n', fontsize=28)
    
    ax.legend(handles= [Patch(facecolor='red', alpha = 0.5, edgecolor='red',label=f'(1) PD predictions ({len(Predictions_pubmed)})'),
                Patch(facecolor='cornflowerblue', alpha = 0.5, edgecolor='cornflowerblue',label=f'(2) background ({len(BG_pubmed)})')],
              loc=1, labelspacing=1, prop={'size': 28}, facecolor ='white')
    
    if day != 'Core PD predictions':
        plt.savefig(f'{save_dir}/C{day}PD{day}_PDpredsvsOG_MWU.jpg', dpi = 350, format='jpg')        
    else:
        plt.savefig(f'{save_dir}/CorePDpreds_PDpredsvsOG_MWU.jpg', dpi = 350, format='jpg')        
                
    plt.show()
    plt.close()



####### MAIN CODE
day_dict = {'0':'IPSCs', '6':'D06', '15':'D15', '21':'D21'}
from datetime import date

today = date.today()
d1 = today.strftime("%d_%m_%Y")

with open('input/LitValid_AllGenes.pkl', 'rb') as handle:
    LitValid_AllGenes = pickle.load(handle) 
# remove PD genes from background to make it consistent with the background from where you choose the predictions

All_genes = []
for root, dirs, files in os.walk('input/All_genes'):
    for file in files:
        if not file.endswith('.csv'):
            with open(f'{root}/{file}', 'rb') as handle:
                cc1cc2_allgenes = pickle.load(handle) 
        All_genes.append(cc1cc2_allgenes)
All_genes = set().union(*All_genes)
print(len(All_genes))
with open('output/All_genes_Union.txt', 'w') as f: 
    for gene in All_genes:
        f.write(f'{gene}\n')   
        
LitValid_AllGenes = {key: LitValid_AllGenes[key] for key in All_genes}

# All_litValidGenes = {}


### Stage Specific validation 
file_out = 'LitValid_Enirch.txt'
PDpreds_enrich = []
if not os.path.exists('output/StageSpecPreds'):
    os.makedirs('output/StageSpecPreds')  
outf = open(f'output/StageSpecPreds/{file_out}','w')

for root, dirs, files in os.walk('input/StageSpecPreds'):
    for file in files:
        if file.endswith('.pkl'):
            print(root, file)
            sd = root.split('/')[1]
            filen = file.split('.')[0]
            save_dir = f'output/{sd}'
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)  
            day = filen.split('C')[1].split('P')[0]
            cc_pair = file.split('_')[0]
            
            
            with open(f'{root}/{file}', 'rb') as handel:
                cc1cc2_PD_preds = pickle.load(handel) #StageSPec Preds
            cc1cc2_PDpreds_ED = pd.read_csv(f'{root}/EuclideanDistance_kMeans_Eucl_dist/{filen}.csv', header = 0, index_col = 0)
 
            cc1cc2_PDpreds_LitValid = dict((k, LitValid_AllGenes[k]) for k in cc1cc2_PD_preds if k in LitValid_AllGenes)                              
            with open(f'{save_dir}/{filen}.pkl', 'wb') as fp:   
                pickle.dump(cc1cc2_PDpreds_LitValid, fp)
            with open(f'{save_dir}/{filen}.txt', 'w') as f:
                for gene in cc1cc2_PDpreds_ED.index:
                    f.write(f'{gene}\t{cc1cc2_PDpreds_LitValid[gene][0]}\t{cc1cc2_PDpreds_LitValid[gene][1]}\n')
    

            #check with background             
            save_dir2 = 'output/StageSpecPreds/Other_genes'
            if not os.path.exists(save_dir2):
                os.makedirs(save_dir2)  
 
            with open(f'input/All_genes/{cc_pair}_All_genes.pkl', 'rb') as handel:
                cc1cc2_AllGenes = pickle.load(handel)                
            cc1cc2_OtherGenes = [gene for gene in cc1cc2_AllGenes if gene not in cc1cc2_PD_preds]

            cc1cc2_OtherGenes_LitValid = dict((k, LitValid_AllGenes[k]) for k in cc1cc2_OtherGenes if k in LitValid_AllGenes)                    
            with open(f'{save_dir2}/{cc_pair}_OGs_LitValid.pkl', 'wb') as fp:   
                pickle.dump(cc1cc2_OtherGenes_LitValid, fp)
                
                
            # MWU test - shows if clustered genes have higher number of citations 
            for gene in PD_genes:
                cc1cc2_OtherGenes_LitValid.pop(gene, None)
            cc1cc2_LitValid = {**cc1cc2_PDpreds_LitValid, **cc1cc2_OtherGenes_LitValid}
            cc1cc2_pubmed, cc1cc2_genes_pubmed = logPubmedn(cc1cc2_PDpreds_LitValid)
            cc1cc2_OtherGenes_pubmed, _ = logPubmedn(cc1cc2_OtherGenes_LitValid)

            PDGs2Background_dist(cc1cc2_pubmed, cc1cc2_OtherGenes_pubmed, day, save_dir)
            
            ### check enrichment of PD preds with PD proof from Pubmed and Gene4PD, 
            case = filen.split('_')[0]
            outf.write(f'{filen}\n')

            LitEnr_genes = LitEnrich(cc1cc2_LitValid, cc1cc2_PDpreds_LitValid, outf, case)
            PDpreds_enrich.append(LitEnr_genes) 

outf.close()

### plot Validated Genes and p values of enrichments of Stage Specific PD preds
BarPlotValidGenes(PDpreds_enrich, file_out)



######### Core PD predictions Validation
save_dir = 'output/CorePreds'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)  

# Create Background set of Other genes; Genes expressed across all conditinos - includes PD genes from DisGeNet

Entrez2Symbol_file = open('input/Homo_sapiens.gene_info', 'r')
Entrez2Sym = Entrez2Symbols(Entrez2Symbol_file)

all_genes = []
for root, dirs, files in os.walk('input/All_genes'):
    for file in files:
        if file.endswith('.csv') and 'Geneslist' in file:
            cc_genes = []
            with open(f'{root}/{file}') as csv_file:
                csv_reader = csv.reader(csv_file, delimiter='\n')
                next(csv_reader)                
                for row in csv_reader:
                    cc_genes.append(row[0])
            all_genes.append(cc_genes)    
common_genes = set.intersection(*map(set,all_genes)) #genes expressed across all cell conditions
with open('output/CommonlyExpressedGenes.txt', 'w') as f: 
    for gene in common_genes:
        f.write(f'{gene}\n')   
        
with open('input/CorePreds.pkl', 'rb') as handel:
    CGs_PD_preds = pickle.load(handel)



# VENN between DEGs from the original scRNA-seq study and Core PD preds
with open('input/DEGs_Skupin.txt') as f:
    DEGs = f.readlines()
DEGs = [x[:-1] for x in DEGs]
DEGs_proteincoding = [x for x in DEGs if x in LitValid_AllGenes]

genes = [set(DEGs_proteincoding), set(CGs_PD_preds)]
gene_sets = ['DEGs', 'Core PD predictions']
gene_sets_dict = dict(zip(gene_sets, genes)) 

fig, ax = plt.subplots(figsize=(8, 5))
venn.venn(gene_sets_dict, cmap=['gold', 'red'], fontsize=24)
plt.savefig(f'{save_dir}/DEGvsCorePreds.jpg',  dpi = 350, format='jpg', bbox_inches='tight')    
fig.show()


CGs_PD_preds_LitValid = {your_key: LitValid_AllGenes[your_key] for your_key in CGs_PD_preds}
All_CommonGenes_LitValid = {your_key: LitValid_AllGenes[your_key] for your_key in common_genes}
Other_commongenes = [gene for gene in common_genes if gene not in CGs_PD_preds]
Other_CommonGenes_LitValid = {your_key: LitValid_AllGenes[your_key] for your_key in Other_commongenes}

with open(f'{save_dir}/CorePreds_LitValid.txt', 'w') as f:
    for gene in CGs_PD_preds:
        f.write(f'{gene}\t{CGs_PD_preds_LitValid[gene][0]}\t{CGs_PD_preds_LitValid[gene][1]}\n')
with open(f'{save_dir}/CorePreds_LitValid.pkl', 'wb') as fp:   
    pickle.dump(CGs_PD_preds_LitValid, fp)               
with open(f'{save_dir}/All_CommonGenes_LitValid.pkl', 'wb') as fp:   
    pickle.dump(All_CommonGenes_LitValid, fp)

     
for gene in PD_genes:
    All_CommonGenes_LitValid.pop(gene, None)
for gene in PD_genes:
    Other_CommonGenes_LitValid.pop(gene, None)

CGs_PD_preds_pubmed,_ = logPubmedn(CGs_PD_preds_LitValid)
OtherGenes_pubmed,_ = logPubmedn(Other_CommonGenes_LitValid) #all genes a OtherGenes_LitValid
PDGs2Background_dist(CGs_PD_preds_pubmed, OtherGenes_pubmed, 'Core PD predictions', save_dir)


### check enrichment of Cire PD preds with PD proof from Pubmed and Gene4PD, 
print('Core PD Preds')
file_out = 'CorePreds_LitValid_Enirch.txt'
outf = open(f'{save_dir}/{file_out}','w')
LitEnr_genes = LitEnrich(All_CommonGenes_LitValid, CGs_PD_preds_LitValid, outf, 'Core PD preds')
outf.close()

# plot Lit enrichment of Core pd preds 
width = 0.1    
opacity = 0.7
fig, ax = plt.subplots(figsize=(6, 9))
# Remove top and right border
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
# Remove y-axis tick marks
ax.yaxis.set_ticks_position('none')
# Set plot title    
ax.set_title('Enrichment of PD predictions in PD-related genes', fontsize = 28, pad=20)

# Add major gridlines in the y-axis
ax.grid(color='black', axis='y', linestyle='-', linewidth=0.4, alpha=0.9)
ax.set_facecolor('xkcd:white')

bar1 = ax.bar(np.arange(1), LitEnr_genes[1], width, align='center', alpha=opacity, color='r', label='Core PD predictions')
bar2 = ax.bar(np.arange(1) + width, LitEnr_genes[2], width, align='center', alpha=opacity,  color='limegreen', label='PD-related genes')
 
ax.legend(fontsize=22, loc = 'upper left')
ax.set_ylabel('#genes', fontsize = 24, fontweight = 'bold')
ax.set_xticks(np.arange(1)+width/2)
ax.set_xticklabels([LitEnr_genes[0]],  fontsize = 24) 
ax.tick_params(axis='y', which='major', labelsize=24)
plt.ylim(0,320)

pval = LitEnr_genes[3]
# pval = "{:.2f}".format(pval)
for rect in bar1:
    height = rect.get_height()
    if float(pval) <= 0.05:
        plt.text(rect.get_x() + rect.get_width(), height+0.1, "{:.3e}".format(pval) + '*', ha='center', va='bottom', fontsize = 24)
    else:
        plt.text(rect.get_x() + rect.get_width(), height+0.1, "{:.3e}".format(pval), ha='center', va='bottom', fontsize = 24)

plt.savefig(f'{save_dir}/LitValid_Enrich.jpg',  dpi = 350, format='jpg', bbox_inches='tight')   
plt.show()
plt.close()

  