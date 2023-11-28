import math
from scipy.stats import hypergeom, shapiro
import numpy as np
import os.path
import pandas as pd
import matplotlib.pyplot as plt
import statistics, sys


# OBJ_fig.savefig(k1k2_folder + '/' + nets + '_OBJ.jpg', dpi = 350, format='jpg')

# Benjamini-Hochberg p-value correction for multiple hypothesis testing
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

def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    stjpg = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(stjpg * p[by_descend]))
    return q[by_orig]
    #return p



# Read GO terms annotations
def load_GO(fname, genes):
    geneset=set(genes)
    GO_counter = {}
    gene2GO = {}
    #loading raw gene & go relationships
    fRead = open(fname, "r")
    for line in fRead.readlines():
        lspt = line.strip().split()
        if len(lspt)>1:
            try:
                geneid = Entrez2Sym[lspt[0]]
                term = lspt[1]
                if geneid in geneset:
                    if term not in GO_counter:
                        GO_counter[term] = 0
                    GO_counter[term] += 1
                    if geneid not in gene2GO:
                        gene2GO[geneid] = set()
                    gene2GO[geneid].add(term)
            except KeyError:
                pass
    fRead.close()
    #print(GO_counter)
    #filter out GO terms that annotates only one gene
    GO_counter2 = {}
    gene2GO2 = {}
    removed = set()
    for term in GO_counter:
        if GO_counter[term] > 1:
            GO_counter2[term] = GO_counter[term]
        else:
            removed.add(term)
    for gene in gene2GO:
        genego = gene2GO[gene].difference(removed)
        if len(genego)>0:
            gene2GO2[gene]=genego
    return [GO_counter2, gene2GO2]

def go_enrichment(clusters, gene2go, go_counts):
    M = len(gene2go) # total number of annotated genes in the dataset
    data = []
    i = -1
    enrichment = [[] for j in range(len(clusters))]
    cts = 0
    
    clts = []
    gos = []
    pvals = []
    
    NE_list = []
    
    for cluster in clusters:
        i+=1
        annotated_genes = []
        for gene in cluster:
            if gene in gene2go:
                annotated_genes.append(gene)
        N = len(annotated_genes) # number of annotated genes in the cluster
        #print N
        annotation_set = set()
        for gene in annotated_genes:
            for term in gene2go[gene]:
                annotation_set.add(term)
        #print len(annotation_set)
        for term in annotation_set:
            K = go_counts[term] # number of genes annotated with the given term in all the data
            X = 0   # number of gene in the clusters that are annotated with the go term
            for gene in annotated_genes:
                if term in gene2go[gene]:
                    X += 1
            pval = hypergeom.sf(X-1, M, K, N) # faster and more accurate than 1-cdf
            clts.append(i)
            gos.append(term)
            pvals.append(pval)
            if pval <= 0.05:
                cts += 1
                #print "Cluster %i, term %s: X=%i, K=%i, pval = %s"%(i,term,X,K,str(pval))
                enrichment[i].append([term,pval])
            #print "%s %s"%(term,prb)

        #print(len(enrichment))    
        #Mean normalized entropy:
        d = float(len(annotation_set))  # number of different annotations in cluster c
        nec = 0.
        if d>1.:
            Cl = float(len(cluster))    # number of gene in cluster c
            sum_pi = 0.
            #counting map
            counting_map = {}
            for gene in cluster:
                if gene in gene2go:
                    for go in gene2go[gene]:
                        if go not in counting_map:
                            counting_map[go] = 0.
                        counting_map[go] += 1.
            for go in counting_map.keys():
                pi = counting_map[go]/Cl    # percentage of genes in c annotated with the considered go term
                sum_pi += pi*math.log(pi)
            nec = (-1./(math.log(d)))*sum_pi
        NE_list.append( nec )
    
    #applying BH correction (Benjamini-Hochner correction)
    BHpvals = p_adjust_bh(pvals)

    BHcts = 0
    BH_enrichment = [[] for j in range(len(clusters))]
    enriched_genes = []
    for i in range(len(clts)):
        #print pvals[i], BHpvals[i]
        if BHpvals[i]<=0.05:
            BHcts += 1
            BH_enrichment[clts[i]].append([gos[i],BHpvals[i]])
    for i in range(len(clusters)):
            cluster_set = set()
            enriched_gos = BH_enrichment[i]
            for gene in clusters[i]:
                for go, pval in enriched_gos:
                    if gene in gene2go:
                        if go in gene2go[gene]:
                            cluster_set.add(gene)
            enriched_genes.append(list(cluster_set)) #genes that are enriched in each cluster
    
    #print cts, BHcts
    MNE = sum(NE_list)/float(len(NE_list))
    #print "Mean Normalized Entropy = %s"%(str(MNE))
    #print(BH_enrichment)
    enr_cluster = 0
    total_cluster = 0
    for i in range(len(clusters)):
        if len(clusters[i])>0:
            total_cluster += 1.
            if len(BH_enrichment[i])>0:
                enr_cluster += 1
    perc_enr_cluster = 100.*enr_cluster/total_cluster
    perc_genes = sum([len(enriched_genes[i]) for i in range(len(clusters))]) #number of genes enrihced in all clusters together (sum number of genes for each individual cluster)
    #print(perc_genes)
    perc_genes = 100.*perc_genes/float(len(gene2go))
    #print(perc_genes, float(len(gene2go)))    
    return [BH_enrichment, MNE, perc_genes, perc_enr_cluster]



def avg_percentile_enrich(ER_ccs):
    for i in range(len(ER_ccs)):
        mean_val = np.percentile(ER_ccs[i], 50)
        ER_ccs[i].append(mean_val)
        
        q3, q1 = np.percentile(ER_ccs[i], [84.13, 15.87])
        min_error = mean_val - q1
        max_error = q3 - mean_val 
        ER_ccs[i].append(min_error)
        ER_ccs[i].append(max_error)

    return ER_ccs



def plotGeneClusters(all_nets, ER_list_BP_ccs, ER_list_KP_ccs, ER_list_RP_ccs, case, cell_cond, clust_meth):
    capsize = 5 
    alpha = 0.7
    N = len(all_nets)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.2       # the width of the bars
    if cell_cond != 'All_Cell_conditions':
        save_dir = f'output/{cell_cond}'
    else:
        save_dir = 'output'
    
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.grid(color='black', axis='y', linestyle='-', linewidth=0.4, alpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    for i in range(N):
        # print(all_nets[i-1], ae1, ER_list_KP_ccs[i-1])
        rects1 = ax.bar(ind[i-1], ER_list_KP_ccs[i-1], width=width, color='r', alpha = alpha, ecolor='black', capsize=capsize)
        rects2 = ax.bar(ind[i-1]+width, ER_list_RP_ccs[i-1], width=width, color='g', alpha = alpha, ecolor='black', capsize=capsize)
        rects3 = ax.bar(ind[i-1]+width*2, ER_list_BP_ccs[i-1], width=width, color='b', alpha = alpha, ecolor='black', capsize=capsize)

        # print(all_nets[i-1])
        # print(ER_clusters_KP[i-1])
        # print(ER_clusters_RP[i-1])
        # print(ER_clusters_BP[i-1])
        # print('\n')    
        
    ax.set_ylabel(f'{case} with enriched annotations (%)', fontsize = 20, fontweight = 'bold')
    ax.set_xticks(ind+width)
    if cell_cond != 'All_Cell_conditions':
        if labels_key[cell_cond] == 'C21' or labels_key[cell_cond] == 'PD21':
            ax.set_xticklabels(all_nets,  fontsize = 18, rotation=90) #nets_everythinga should be all_nets
        else:
            ax.set_xticklabels(len(all_nets)*[''])
        plt.title(f'{labels_key[cell_cond]} - {clust_meth}',  fontsize = 30, pad = 24)

    else:
        ax.set_xticklabels(all_nets,  fontsize = 18, rotation=90) #nets_everythinga should be all_nets
        plt.title(f'Average over all Cell conditions - {clust_meth}',  fontsize = 30, pad = 24)

    ax.tick_params(axis='y', which='major', labelsize=16)
    ax.legend( (rects1[0], rects2[0], rects3[0]), ('KP', 'RP', 'GO'),  fontsize = 18)

    if case == 'Clusters':       
        plt.ylim(0,105)
        plt.tight_layout()    
        plt.savefig(f'{save_dir}/{cell_cond}_ECs_{clust_meth}.jpg', dpi = 350, format='jpg')    
    elif case == 'Genes':   
        plt.ylim(0,60)
        plt.tight_layout()    
        plt.savefig(f'{save_dir}/{cell_cond}_EGs_{clust_meth}.jpg', dpi = 350, format='jpg')    
    plt.show()  
    plt.close()

def plotGeneClusters_avg(all_nets, ER_list_BP_ccs, ER_list_KP_ccs, ER_list_RP_ccs, case, cell_cond, clust_meth):
    capsize = 5 
    alpha_b = 0.8
    N = len(all_nets)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.2       # the width of the bars

    save_dir = 'output'
    
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.grid(color='black', axis='y', linestyle='-', linewidth=0.4, alpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    for i in range(N):
        # print(all_nets[i-1], ae1, ER_list_KP_ccs[i-1])
        ae1 = np.array(list(zip([ER_list_KP_ccs[i-1][-2]], [ER_list_KP_ccs[i-1][-1]]))).T
        rects1 = ax.bar(ind[i-1], ER_list_KP_ccs[i-1][-3], yerr=ae1, width=width, color='r', alpha = alpha_b, ecolor='black', capsize=capsize)
        ae2 = np.array(list(zip([ER_list_RP_ccs[i-1][-2]], [ER_list_RP_ccs[i-1][-1]]))).T
        rects2 = ax.bar(ind[i-1]+width, ER_list_RP_ccs[i-1][-3], yerr=ae2, width=width, color='g', alpha = alpha_b, ecolor='black', capsize=capsize)
        ae3 = np.array(list(zip([ER_list_BP_ccs[i-1][-2]], [ER_list_BP_ccs[i-1][-1]]))).T
        rects3 = ax.bar(ind[i-1]+width*2, ER_list_BP_ccs[i-1][-3], yerr=ae3, width=width, color='b', alpha = alpha_b, ecolor='black', capsize=capsize)

        # print(all_nets[i-1])
        # print(ER_clusters_KP[i-1])
        # print(ER_clusters_RP[i-1])
        # print(ER_clusters_BP[i-1])
        # print('\n')    
        
    ax.set_ylabel(f'{case} with enriched annotations (%)', fontsize = 24, fontweight = 'bold', loc = 'top')
    ax.set_xticks(ind+width)

    if cell_cond == 'All_Cell_conditions ':
        if case == 'Genes':
            ax.set_xticklabels(N*'',  fontsize = 28, rotation=10) #nets_everythinga should be all_nets
        else:
            ax.set_xticklabels(N*'',  fontsize = 28, rotation=10) #nets_everythinga should be all_nets
            plt.title('Average over all Cell conditions',  fontsize = 30, pad = 24)

    elif cell_cond == 'All_Cell_conditions_3scen':
        ax.set_xticklabels(all_nets,  fontsize = 28) #nets_everythinga should be all_nets        
    else:
        if case == 'Genes':
            ax.set_xticklabels(all_nets,  fontsize = 28, rotation=90) #nets_everythinga should be all_nets
        else:
            ax.set_xticklabels(N*'',  fontsize = 28, rotation=90) #nets_everythinga should be all_nets
       
    # plt.title(f'Average over all Cell conditions - {clust_meth}',  fontsize = 30, pad = 24)

    ax.tick_params(axis='y', which='major', labelsize=28)
    ax.legend( (rects1[0], rects2[0], rects3[0]), ('KP', 'RP', 'GO'),  fontsize = 28, loc = 2)

    if case == 'Clusters':       
        plt.ylim(0,140)
        plt.tight_layout()    
        plt.savefig(f'{save_dir}/{cell_cond}_ECs_{clust_meth}.jpg', dpi = 350, format='jpg', bbox_inches='tight')   
    elif case == 'Genes':   
        plt.ylim(0,80)
        plt.tight_layout()    
        plt.savefig(f'{save_dir}/{cell_cond}_EGs_{clust_meth}.jpg', dpi = 350, format='jpg', bbox_inches='tight')   
    plt.show()  
    plt.close()    
            
            
"""
    Main code starts here
"""
import pickle

indir = 'input'

nets_1 = ['PPI','GI','MI','COEX']
nets_2 = ['PPI+GI','PPI+MI','PPI+COEX','COEX+GI','COEX+MI','GI+MI']
nets_3 = ['GI+COEX+MI','PPI+GI+MI','PPI+COEX+MI','PPI+GI+COEX']
nets_everything = ['EXP'] + nets_1 + nets_2 + nets_3 + ['ALL']

labels_key = {'Control_IPSCs':'C0', 'Control_D06':'C6', 'Control_D15':'C15', 'Control_D21':'C21', 
              'PINK1_IPSCs':'PD0', 'PINK1_D06':'PD6', 'PINK1_D15':'PD15', 'PINK1_D21':'PD21'}

nets_leg = {'EXP':'E', 'PPI':'E+P','GI':'E+G','MI':'E+M','COEX':'E+C',
            'PPI+GI':'E+P+G','PPI+MI':'E+P+M','PPI+COEX':'E+P+C','COEX+GI':'E+C+G','COEX+MI':'E+C+M','GI+MI':'E+G+M',
            'GI+COEX+MI':'E+G+C+M','PPI+GI+MI':'E+P+G+M','PPI+COEX+MI':'E+P+C+M','PPI+GI+COEX':'E+P+G+C', 'ALL':'E+P+G+C+M', 'NetSC-NMTF':'E+P+G+C+M'}


cell_conds = labels_key.keys()
Flag = bool(sys.argv[1])
if Flag == True:
    print(Flag)
    # Computing Enrichments
    for cell_cond in cell_conds:  
        print(cell_cond)
        save_dir = f'output/{cell_cond}'
        try:
            os.makedirs(save_dir)
        except FileExistsError:
            # directory already exists
            pass
        
        genes = pd.read_csv(f'{indir}/Geneslist/Geneslist_{cell_cond}.csv', header = 0)
        genes = genes.values
        genes = [str(item) for sublist in genes for item in sublist]
        PPIlist = genes
        
        fgoname_kp = f'{indir}/HSA_Kegg_Pathways.lst'
        go_counts_kp,gene2go_kp = load_GO(fgoname_kp, PPIlist)
        print("%i KP annotated genes"%(len(gene2go_kp)))
        
        fgoname_rp = f'{indir}/HSA_Reactome_Pathways.lst'
        go_counts_rp,gene2go_rp = load_GO(fgoname_rp, PPIlist)
        print("%i RP annotated genes"%(len(gene2go_rp)))
        
        fgoname_bp = f'{indir}/HSA_GO-BP.lst'
        go_counts_bp,gene2go_bp = load_GO(fgoname_bp, PPIlist)
        print("%i BP annotated genes"%(len(gene2go_bp)))
        

        ER_clusters_KP = []
        ER_clusters_RP = []
        ER_clusters_BP = []

        ER_genes_KP = []
        ER_genes_RP = []
        ER_genes_BP = []   
        all_nets = []

        for nets in nets_everything:
            for root, dirs, files in os.walk(f'{indir}/Clusters/{cell_cond}'):
                for file in files:
                    lspt = root.split('/')
                    if len(lspt) > 2 and lspt[3] == nets:
                        print(root, file)
                        nets_f = lspt[3]
                        all_nets.append(nets_leg[nets_f])
                        
                        with open(f'{root}/{file}', 'rb') as handle:
                            G1clusters = pickle.load(handle)
                        Cp = [[str(y) for y in x] for x in G1clusters]           
       
                        bh_kp, mne_kp, eg_kp, perc_cluster_kp = go_enrichment(Cp, gene2go_kp, go_counts_kp)
                        bh_rp, mne_rp, eg_rp, perc_cluster_rp = go_enrichment(Cp, gene2go_rp, go_counts_rp)
                        bh_bp, mne_bp, eg_bp, perc_cluster_bp = go_enrichment(Cp, gene2go_bp, go_counts_bp)
                        #print(f'Function Enirchment : {len(bh_bp)}')
       
                        KP_terms = [[i[0] for i in nested] for nested in bh_kp]
                        RP_terms = [[i[0] for i in nested] for nested in bh_rp]
                        BP_terms = [[i[0] for i in nested] for nested in bh_bp]
       
                        print(perc_cluster_kp, eg_kp)
                        print(perc_cluster_rp, eg_rp)
                        print(perc_cluster_bp, eg_bp)
       
                        ER_clusters_KP.append(perc_cluster_kp)
                        ER_clusters_RP.append(perc_cluster_rp)
                        ER_clusters_BP.append(perc_cluster_bp)
       
                        ER_genes_KP.append(eg_kp)
                        ER_genes_RP.append(eg_rp)
                        ER_genes_BP.append(eg_bp)

        clusters_output_file = open(f'{save_dir}/{cell_cond}_kMeans_ECs.txt','w')
        genes_output_file = open(f'{save_dir}/{cell_cond}_kMeans_EGs.txt','w')
        for i in range(len(all_nets)):
            clusters_output_file.write(f'{all_nets[i]}\n{ER_clusters_KP[i]}\n{ER_clusters_RP[i]}\n{ER_clusters_BP[i]}\n')
            genes_output_file.write(f'{all_nets[i]}\n{ER_genes_KP[i]}\n{ER_genes_RP[i]}\n{ER_genes_BP[i]}\n')
        clusters_output_file.close()
        genes_output_file.close()
          

# Plotting 
ER_clusters_BP_avgccs = [[] for i in range(len(nets_everything))]
ER_clusters_RP_avgccs = [[] for i in range(len(nets_everything))]
ER_clusters_KP_avgccs = [[] for i in range(len(nets_everything))]
 
ER_genes_BP_avgccs = [[] for i in range(len(nets_everything))]
ER_genes_RP_avgccs = [[] for i in range(len(nets_everything))]
ER_genes_KP_avgccs = [[] for i in range(len(nets_everything))]

for cell_cond in cell_conds: #['PINK1_D15', 'PINK1_D21']:
    print(cell_cond)
    save_dir = f'output/{cell_cond}'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir) 
    
    all_nets = []
    ER_clusters_KP, ER_clusters_RP, ER_clusters_BP = [],[],[]
    ER_genes_KP, ER_genes_RP, ER_genes_BP = [],[],[]
    with open(f'output/{cell_cond}/{cell_cond}_kMeans_ECs.txt') as ECf:
        # print(ECf.readlines())
        for i in range(len(nets_everything)):
            x = ECf.readline()[:-1]
            if x == 'NetSC-NMTF':
                x = nets_leg[x]      
            all_nets.append(x)    
            ER_clusters_KP.append(float(ECf.readline()[:-1]))
            ER_clusters_RP.append(float(ECf.readline()[:-1]))
            ER_clusters_BP.append(float(ECf.readline()[:-1]))
    with open(f'output/{cell_cond}/{cell_cond}_kMeans_EGs.txt') as EGf:
        for i in range(len(nets_everything)):
            EGf.readline()
            ER_genes_KP.append(float(EGf.readline()[:-1]))
            ER_genes_RP.append(float(EGf.readline()[:-1]))
            ER_genes_BP.append(float(EGf.readline()[:-1]))

    for i in range(len(ER_clusters_KP)):               
        ER_clusters_BP_avgccs[i].append(ER_clusters_BP[i])
        ER_clusters_RP_avgccs[i].append(ER_clusters_RP[i])
        ER_clusters_KP_avgccs[i].append(ER_clusters_KP[i])
         
        ER_genes_BP_avgccs[i].append(ER_genes_BP[i])
        ER_genes_RP_avgccs[i].append(ER_genes_RP[i])
        ER_genes_KP_avgccs[i].append(ER_genes_KP[i])
        
    plotGeneClusters(all_nets, ER_clusters_BP, ER_clusters_KP, ER_clusters_RP, case = 'Clusters', cell_cond = cell_cond, clust_meth = 'kMeans')
    plotGeneClusters(all_nets, ER_genes_BP, ER_genes_KP, ER_genes_RP, case = 'Genes', cell_cond = cell_cond, clust_meth = 'kMeans')
    


avg_percentile_enrich(ER_clusters_BP_avgccs)
avg_percentile_enrich(ER_clusters_RP_avgccs)
avg_percentile_enrich(ER_clusters_KP_avgccs)
avg_percentile_enrich(ER_genes_BP_avgccs)
avg_percentile_enrich(ER_genes_RP_avgccs)
avg_percentile_enrich(ER_genes_KP_avgccs)        

plotGeneClusters_avg(all_nets, ER_clusters_BP_avgccs, ER_clusters_KP_avgccs, ER_clusters_RP_avgccs, case = 'Clusters', cell_cond = 'All_Cell_conditions', clust_meth = 'kMeans')
plotGeneClusters_avg(all_nets, ER_genes_BP_avgccs, ER_genes_KP_avgccs, ER_genes_RP_avgccs, case = 'Genes', cell_cond = 'All_Cell_conditions', clust_meth = 'kMeans')
