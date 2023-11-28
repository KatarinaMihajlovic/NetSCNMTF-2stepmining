# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 15:44:45 2023

@author: kmihajlo
"""
import networkx as nx
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
from random import sample

def Load_Network(fname):
    net = nx.read_edgelist(fname, data = False)
    nodes = [n for n in net.nodes()]
#   nodeset = set(nodes)
    nb_nodes = len(nodes)
    nodes2ind = {}
    for n in range(nb_nodes):
        nodes2ind[nodes[n]]=n
    return net, nodes, nodes2ind

def BasicStats(G, nodes):
    G_CPD = G.subgraph(nodes)
    
    nodes = nx.number_of_nodes(G_CPD)
    edges = nx.number_of_edges(G_CPD)
    density = nx.density(G_CPD)
    
    print(nodes, 'nodes')
    print(edges, 'edges')
    print(density, 'density') 
    
    GLCC = G_CPD.subgraph(max(nx.connected_components(G_CPD), key=len)).copy()
    largest_cc = max(nx.connected_components(G_CPD), key=len)
    print(len(largest_cc), 'largest connected component')
    GLCC = G_CPD.subgraph(largest_cc).copy() 
    print('density: ', nx.classes.function.density(G_CPD))
    clustercoef = nx.average_clustering(GLCC)
    print('Clustering Coefficient (AVG): ', clustercoef) 
    return largest_cc

def ShortestPath2Gene(G, gene, TargetGenes):
    G2gene = nx.single_source_dijkstra(G, gene)
    G2gene_lens = G2gene[0]
    G2gene_lens.pop(gene)

    TargetGenes_SP = {key:G2gene_lens[key] for key in G2gene_lens if key in TargetGenes}
    BG_SP = {key:G2gene_lens[key] for key in G2gene_lens if key not in TargetGenes}
    TargetGenes_SP_dist = np.array([TargetGenes_SP[key] for key in TargetGenes_SP])
    BG_SP_dist = np.array([BG_SP[key] for key in BG_SP])
    return TargetGenes_SP_dist, BG_SP_dist, TargetGenes_SP, BG_SP

def format_int_with_commas(x):
    """
    Formats an integer with commas as thousand separators.
    """
    return f"{x:,}"


def plot2bar(CorePreds_SP_dist, BG_SP_dist, network, save_file, case, log = False):  
    statmwu,pvalmwu = mannwhitneyu(CorePreds_SP_dist, BG_SP_dist, alternative='less')  
    if log:
        CorePreds_SP_cnt = np.bincount(CorePreds_SP_dist)
        BG_SP_cnt = np.bincount(BG_SP_dist)
    else:
        CorePreds_SP_cnt = np.bincount(CorePreds_SP_dist)/len(CorePreds_SP_dist)*100
        BG_SP_cnt = np.bincount(BG_SP_dist)/len(BG_SP_dist)*100

    CorePreds_SP_cnt = np.delete(CorePreds_SP_cnt, 0)
    CorePreds_SPs = np.unique(CorePreds_SP_dist)
    BG_SP_cnt = np.delete(BG_SP_cnt, 0)
    BG_SPs = np.unique(BG_SP_dist)
        
    N = np.max(BG_SP_dist)
    ind = np.arange(N)   
    width = 0.3       
    
    fig, ax = plt.subplots(figsize=(11, 9))
    bg = format_int_with_commas(len(BG_SP_dist))
    cpd = format_int_with_commas(len(CorePreds_SP_dist))
    _ = ax.bar(BG_SPs, BG_SP_cnt, label=f'(1) background ({bg})', width=width, color='cornflowerblue')#, alpha = 0.5)#, marker='o', markersize=14)#, label=f'{PD_Gs} ({len(GM_BGgenes)})') 
    _ = ax.bar(CorePreds_SPs+width, CorePreds_SP_cnt, label=f'(2) {case} ({cpd})', width=width, color='red')#, alpha = 0.5)#, marker='o')#, markersize=14)#, label=f'{PD_Gs} ({len(GM_BGgenes)})') 

    
    if log:
        ax.set_yscale('log')
        ax.set_ylabel('Gene count', fontsize=26)
    else:
        ax.set_ylabel('Relative gene count (%)', fontsize=26)
        if network == 'PPI':
            ax.set_ylim([0, 100])
        elif network == 'COEX':
            ax.set_ylim([0, 120])
           
    ax.set_xlim(xmin=0.5, xmax= np.max(BG_SPs) + 0.5)

    plt.grid(linestyle = '--', linewidth = 0.5)    
    title = f'Shortest path to PINK1 in {network} network'
    plt.title(title, fontsize=26) 
    ax.set_xlabel('Shortest path', fontsize=26) #dist_measure_l[dist_meas]
    plt.xticks(BG_SPs + width / 2, BG_SPs, fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=26)
        
    pvalmwu_s = "{:.3e}".format(pvalmwu)
    if pvalmwu <= 0.05:
        fig.text(0.55, 0.58, f'MWU (2) < (1)\np-value = {pvalmwu_s}*\n', fontsize=24)
    else:
        fig.text(0.55, 0.58, f'MWU (2) < (1)\np-value = {pvalmwu_s}\n', fontsize=24)
    
    plt.savefig(save_file, dpi = 350, format='jpg')                        
    plt.show()
    plt.close()

def plot2boxplot(CorePreds_SP_dist, BG_SP_dist, network, file_name, showfliers =True):
    fig, ax = plt.subplots(figsize=(16, 12))
   
    # Remove top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    # Remove y-axis tick marks
    ax.yaxis.set_ticks_position('none')
    
    # Set plot title
    title = f'Shortest path to PINK1 in {network} network'
    ax.set_title(title, fontsize = 22)
    # Add major gridlines in the y-axis
    ax.grid(color='grey', axis='y', linestyle='-', linewidth=0.25, alpha=0.5)
    ax.set_facecolor('xkcd:white')
 
    colors = ['cornflowerblue', 'red']
    colors_OG = dict(color=colors[0])
    colors_PD = dict(color=colors[1])  
    positions = [0,0.4]
    
    OG = ax.boxplot(BG_SP_dist, positions=[positions[0]], widths=0.4,
                    boxprops=colors_OG, medianprops=colors_OG, whiskerprops=colors_OG, 
                    capprops=colors_OG, flierprops=dict(markeredgecolor=colors[0]),showfliers =showfliers )   
    PD = ax.boxplot(CorePreds_SP_dist, positions=[positions[1]], widths=0.4,
                    boxprops=colors_PD, medianprops=colors_PD, whiskerprops=colors_PD, 
                    capprops=colors_PD, flierprops=dict(markeredgecolor=colors[1]),showfliers =showfliers )  
 
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=12)
    tick_names = ['background', 'DEGs']
 
    ax.set_xticks(positions)
    ax.set_xticklabels(tick_names, rotation = 90)
 
    xt_labels =  ax.get_xticklabels()
 
    plt.savefig(file_name, dpi = 350, format='jpg',  bbox_inches="tight")
    plt.show()
    plt.close()     



def write_txt(genes, sd, save_file):
    with open(f'{sd}/{save_file}.txt', 'w') as f: 
        for gene in genes:                        
            f.write(f'{gene}\n')
        f.write('\n')
        
         
'''MAIN CODE'''
     
with open('input/CorePreds.txt') as f:
    CorePreds = f.read().splitlines() 
with open('input/DEGs_Skupin.txt') as f:
    DEGs_Skupin = f.read().splitlines()    
with open('input/All_genes_Union.txt') as f:
    All_genes_Union = f.read().splitlines() 


DEGs_Skupin = list(set(All_genes_Union) & set(DEGs_Skupin))
print(len(DEGs_Skupin))

target_gene = 'PINK1'
Statistics = open(f'output/Statistics.txt','w')
 
network = 'PPI'   
print(network)
G = nx.read_edgelist('input/PPI_General_Biogrid_GeneSym.edgelist', data = False)
G  = G.subgraph(All_genes_Union)

CorePreds_SP_dist, BG_SP_dist, CorePreds_SP, BG_SP = ShortestPath2Gene(G, target_gene, CorePreds)
plot2bar(CorePreds_SP_dist, BG_SP_dist, network, f'output/ShortestPath2PINK1_{network}_CorePDPredsvsBG.jpg', 'Core PD predictions')

# are the first neighbors significant
first_neigh = [key for key in CorePreds_SP if CorePreds_SP[key] == 1]
perc1st_neigh = len(first_neigh)/len(CorePreds)*100
write_txt(first_neigh, sd='output', save_file=f'{target_gene}_1neighCorePDPreds_{network}')
first_neigh_BG = [key for key in BG_SP if BG_SP[key] == 1]

X = len(first_neigh)
N = len(CorePreds)
K = len(first_neigh_BG) + len(first_neigh)
M = len(G.nodes())
fold = (X/N)/(K/M)
pval = hypergeom.sf(X-1, M, K, N)
print(f'1st Neighbours to {target_gene} in PPI significance: fold = {fold}, pval = {pval}')

Statistics.write(f'{network}\n') 
Statistics.write(f'1st Neighbors to {target_gene}, n = {len(first_neigh)} ({perc1st_neigh}\%)\n')
Statistics.write(f'Are 1st Neighbors of {target_gene} from Core PD predictions significant?\n')
Statistics.write('Test: Hypergeometric test\n')
Statistics.write(f'fold = {fold}, pval = {pval}\t\n')
Statistics.close()



### Novak DEGs
DEGs_Skupin_SP_dist, BG_SP_dist, DEGs_Skupin_SP, BG_SP = ShortestPath2Gene(G, target_gene, DEGs_Skupin)
plot2bar(DEGs_Skupin_SP_dist, BG_SP_dist, network, f'output/ShortestPath2PINK1_{network}_DEGsvsBG.jpg', 'DEGs')

# are the first neighbors significant
first_neigh = [key for key in DEGs_Skupin_SP if DEGs_Skupin_SP[key] == 1]
perc1st_neigh = len(first_neigh)/len(DEGs_Skupin)*100
write_txt(first_neigh, sd='output', save_file=f'{target_gene}_1neighDEGs_Skupin_{network}')
first_neigh_BG = [key for key in BG_SP if BG_SP[key] == 1]

X = len(first_neigh)
N = len(DEGs_Skupin)
K = len(first_neigh_BG) + len(first_neigh)
M = len(G.nodes())
fold = (X/N)/(K/M)
pval = hypergeom.sf(X-1, M, K, N)
print(f'1st Neighbours to {target_gene} in PPI significance: fold = {fold}, pval = {pval}')

