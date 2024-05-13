import math
from scipy.stats import hypergeom
import numpy as np
import os.path
import pandas as pd


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
Entrez2Symbol_file.close()

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
			X = 0	# number of gene in the clusters that are annotated with the go term
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
		d = float(len(annotation_set))	# number of different annotations in cluster c
		nec = 0.
		if d>1.:
			Cl = float(len(cluster))	# number of gene in cluster c
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
				pi = counting_map[go]/Cl	# percentage of genes in c annotated with the considered go term
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




def Test_CorrectClust(G1clusters, gene2go_true, bh, Test_annotated):
    Test_CorrectClust = []
    for gene in Test_annotated:
       for j in range(len(G1clusters)):
           if gene in G1clusters[j]:
                gene_annots = gene2go_true[gene]
                clust_Enrannots = bh[j]
                overlap_annots = list(set(gene_annots)&set([x[0] for x in clust_Enrannots]))
                if len(overlap_annots) > 0:
                    Test_CorrectClust.append(gene)                          
    Test_percCorrectClust = len(Test_CorrectClust)/len(Test_annotated)*100
    return Test_percCorrectClust

	        
"""
	Main code starts here
"""
import pickle

indir = 'input/Clusters'

labels_key = {'Control_IPSCs':'C0', 'Control_D06':'C6', 'Control_D15':'C15', 'Control_D21':'C21', 
              'PINK1_IPSCs':'PD0', 'PINK1_D06':'PD6', 'PINK1_D15':'PD15', 'PINK1_D21':'PD21'}
labels_key2 = {'C0':'Control_IPSCs', 'C6':'Control_D06', 'C15':'Control_D15', 'C21':'Control_D21', 
              'PD0':'PINK1_IPSCs', 'PD6':'PINK1_D06', 'PD15':'PINK1_D15', 'PD21':'PINK1_D21'}




Flag = True
if Flag == True:
    # Computing Enrichments
    KPperc = []
    RPperc = []
    BPperc = []
    
    for cc in os.listdir(indir):  
        print(cc)
        
        cc_KPperc = []
        cc_RPperc = []
        cc_BPperc = []
        
        genelist = pd.read_csv(f'input/Geneslist/Geneslist_{cc}.csv', header=0, index_col=0)
        genelist = list(genelist.index)     
        fgoname_kp = 'input/HSA_Kegg_Pathways.lst'
        go_counts_kp_true,gene2go_kp_true = load_GO(fgoname_kp, genelist)
        print("%i KP annotated genes"%(len(gene2go_kp_true)))
        
        fgoname_rp = 'input//HSA_Reactome_Pathways.lst'
        go_counts_rp_true,gene2go_rp_true = load_GO(fgoname_rp, genelist)
        print("%i RP annotated genes"%(len(gene2go_rp_true)))
        
        fgoname_bp = 'input//HSA_GO-BP.lst'
        go_counts_bp_true,gene2go_bp_true = load_GO(fgoname_bp, genelist)
        print("%i BP annotated genes"%(len(gene2go_bp_true)))
        
            
        with open(f'input/Folds/{labels_key[cc]}_Train_Test5FOLD.pkl', 'rb') as handle:
            Train_Test5FOLD = pickle.load(handle)   
    
        for i in range(len(Train_Test5FOLD)):
            fold_KPperc = []
            fold_RPperc = []
            fold_BPperc = []      
            
            Train = Train_Test5FOLD[i][0]
            Test = Train_Test5FOLD[i][1]
          
            Test_KPannotated = list(set(Test)&set(gene2go_kp_true.keys()))
            Test_RPannotated = list(set(Test)&set(gene2go_rp_true.keys()))
            Test_BPannotated = list(set(Test)&set(gene2go_bp_true.keys()))
    
            go_counts_kp,gene2go_kp = load_GO(fgoname_kp, Train)
            print("%i KP annotated genes"%(len(gene2go_kp)))
            
            go_counts_rp,gene2go_rp = load_GO(fgoname_rp, Train)
            print("%i RP annotated genes"%(len(gene2go_rp)))
            
            go_counts_bp,gene2go_bp = load_GO(fgoname_bp, Train)
            print("%i BP annotated genes"%(len(gene2go_bp)))
            
            for file in os.listdir(f'input/Clusters/{cc}/FOLD{i}'):
                print(file)
                with open(f'input/Clusters/{cc}/FOLD{i}/{file}', 'rb') as handle:
                    G1clusters = pickle.load(handle)

                bh_kp, mne_kp, eg_kp, perc_cluster_kp = go_enrichment(G1clusters, gene2go_kp, go_counts_kp)
                bh_rp, mne_rp, eg_rp, perc_cluster_rp = go_enrichment(G1clusters, gene2go_rp, go_counts_rp)
                bh_bp, mne_bp, eg_bp, perc_cluster_bp = go_enrichment(G1clusters, gene2go_bp, go_counts_bp)
                
                Test_percKPCorrectClust = Test_CorrectClust(G1clusters, gene2go_kp_true, bh_kp, Test_KPannotated)
                Test_percRPCorrectClust = Test_CorrectClust(G1clusters, gene2go_rp_true, bh_rp, Test_RPannotated)
                Test_percBPCorrectClust = Test_CorrectClust(G1clusters, gene2go_bp_true, bh_bp, Test_BPannotated)
                
                fold_KPperc.append(Test_percKPCorrectClust)
                fold_RPperc.append(Test_percRPCorrectClust)
                fold_BPperc.append(Test_percBPCorrectClust)
                print(len(fold_KPperc))       
                
                KPperc.append(Test_percKPCorrectClust)
                RPperc.append(Test_percRPCorrectClust)
                BPperc.append(Test_percBPCorrectClust)
                print(len(KPperc))       

    
            cc_KPperc.append(np.array(fold_KPperc))
            cc_RPperc.append(np.array(fold_RPperc))
            cc_BPperc.append(np.array(fold_BPperc))
            print(len(cc_KPperc))
    
        np.save(f'output/{cc}_PartMImissing_KPperc_TestCorrectClst.npy', np.array(cc_KPperc))
        np.save(f'output/{cc}_PartMImissing_RPperc_TestCorrectClst.npy', np.array(cc_RPperc))
        np.save(f'output/{cc}_PartMImissing_BPperc_TestCorrectClst.npy', np.array(cc_BPperc))
    
        
    
    np.save('output/PartMImissing_KPperc_TestCorrectClst.npy', np.array(KPperc))
    np.save('output/PartMImissing_RPperc_TestCorrectClst.npy', np.array(RPperc))
    np.save('output/PartMImissing_BPperc_TestCorrectClst.npy', np.array(BPperc))
    
print(np.median(KPperc))    
print(np.median(BPperc))
print(np.median(RPperc))


    
    


