#import networkx as nx
import math, json
from scipy.stats import hypergeom
import numpy as np
import matplotlib.pyplot as plt
import os.path
import pandas as pd
from matplotlib_venn import venn2, venn2_circles


def parseKEGGCats(file, KEGG2names):
    with open(file) as json_file:
        data = json.load(json_file)
    
    resAll = []
    for i in range(len(data['children'])):
        name = ' '.join(data['children'][i]['name'].split(' ')[1:])
        
        for subgroup in data['children'][i]['children']:
            subgroupName = ' '.join(subgroup['name'].split(' ')[1:])
            for s in subgroup['children']:
                resTmp = [name, subgroupName]
                
                if 'PATH' in s['name']:
                    resTmp.append(s['name'].split(':')[-1][:-1])
                    resTmp.append(KEGG2names[s['name'].split(':')[-1][:-1]])
                    resAll.append(resTmp)
            
    df = pd.DataFrame(resAll, columns=['Group', 'Subgroup', 'Pathway', 'PathwayName'])
    return df

def Entrez2Symbols(gene_info): #dict with entrez as keys and symbol as values (key:value)
    Entrez2Sym = {}
    gene_info.readline()
    for line in gene_info.readlines():
        lspt=line.strip().split('\t')
        Symbol = lspt[2]
        Entrez2Sym[lspt[1]] = Symbol
    return(Entrez2Sym)


# Benjamini-Hochberg p-value correction for multiple hypothesis testing
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
    #return p

# Read gene ontology
def load_gene_ontology(fname):
	go2name = {}
	ifile = open(fname,'r')
	go = ""
	name = ""
	for line in ifile.readlines():
		try:
			if len(line)>3:
				if line[:3] == "id:":
					lspt = line.strip().split(':')
					go = "GO:"+lspt[2]
			if len(line)>5:
				if line[:5] == "name:":
					lspt = line.strip().split(':')
					name = lspt[1]
			if go != "" and name != "":
				go2name[go]=name
				go = ""
				name = ""
		except:
			print( line)
	ifile.close()
	#print "%i annotation loaded"%(len(go2name))
	return go2name


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
			geneid = lspt[0]
			term = lspt[1]
			try:
				geneid = Entrez2Sym[geneid]
			except KeyError:
				geneid = ''
			if geneid in geneset:
				if term not in GO_counter:
					GO_counter[term] = 0
				GO_counter[term] += 1
				if geneid not in gene2GO:
					gene2GO[geneid] = set()
				gene2GO[geneid].add(term)
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
	
	cls = []
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
			cls.append(i)
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
	for i in range(len(cls)):
		#print pvals[i], BHpvals[i]
		if BHpvals[i]<=0.05:
			BHcts += 1
			BH_enrichment[cls[i]].append([gos[i],BHpvals[i]])
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
# 	print(len(enriched_genes[0]), len(clusters))
# 	print(perc_genes)
    #print(perc_genes)
	perc_genes = 100.*perc_genes/float(len(gene2go))
	#print(perc_genes, float(len(gene2go)))    
	return [BH_enrichment, MNE, perc_genes, perc_enr_cluster, enriched_genes]


def parseKEGGnames(KEGGnames_file):
    KEGG2names = {}
    for line in KEGGnames_file.readlines():
        lspt=line.strip().split('\t')
        KEGG_name = lspt[0].split(':')[1]
        name = lspt[1].split('Homo')[0]
        KEGG2names[KEGG_name] = name
    return(KEGG2names)    



def Sort(sub_li):
  
    # reverse = None (Sorts in Ascending order)
    # key is set to sort using second element of 
    # sublist lambda has been used
    sub_li.sort(key = lambda x: len(x[1]), reverse = True)
    return sub_li

"""
	Main code starts here
"""
from goatools.obo_parser import GODag
import pickle

indir = 'input'
Entrez2Symbol_file = open('input/Homo_sapiens.gene_info', 'r')
Entrez2Sym = Entrez2Symbols(Entrez2Symbol_file)   
KEGGnames_file = open('input/hsa_pathway_names.lst', 'r')
KEGG2names = parseKEGGnames(KEGGnames_file)



All_genes = []
for root, dirs, files in os.walk('input/All_genes'):
    for file in files:
        if not file.endswith('.csv'):
            with open(f'{root}/{file}', 'rb') as handle:
                cc1cc2_allgenes = pickle.load(handle) 
        All_genes += cc1cc2_allgenes
All_genes = list(dict.fromkeys(All_genes))

fgoname_kp = f'{indir}/HSA_Kegg_Pathways.lst'
go_counts_kp,gene2go_kp = load_GO(fgoname_kp, All_genes)
print("%i KP annotated genes"%(len(gene2go_kp)))
   

     
topnCGs_files = []
for root, dirs, files in os.walk(indir):
    for file in files:
        if 'topnCGs' in file:# or 'PD_genes_DGN' in file:
            topnCGs_files.append(file) 
        
        
# for root, dirs, files in os.walk(indir):
#     for file in files:
for file in topnCGs_files:    
    # print(file)
    KP_terms_2sets = []      
    KP_termsmeaning_2sets = []    
         
    sd = save_file = file[:-4]
    print(sd)
    with open(f'input/{file}', 'rb') as handle:
        PD_genes = pickle.load(handle)  
    PD_genes = list(PD_genes.keys())
    PD_preds_common = PD_genes
        
    if not os.path.exists(f'output/{sd}'):
        os.makedirs(f'output/{sd}') 
    
    bh_kp, mne_kp, eg_kp, perc_cluster_kp, enriched_genes_kp = go_enrichment([PD_genes], gene2go_kp, go_counts_kp)      
    enriched_genes_kp = enriched_genes_kp[0]
         
    KP_terms = [[i[0] for i in nested] for nested in bh_kp]
    KP_terms_p = [[i[1] for i in nested] for nested in bh_kp][0]

  

    #write KEGG terms
    KP_terms_meaning = []
    for KP_term in KP_terms[0]:
        try:
            KP_terms_meaning.append(KEGG2names[KP_term])
        except KeyError:
            print(KP_term)
    KP_terms_2sets.append(KP_terms[0])
    KP_termsmeaning_2sets.append(KP_terms_meaning)

    KP_terms_meaning_p =  [list(i) for i in zip(KP_terms_meaning, KP_terms_p)]
    KP_terms_meaning_p.sort(key = lambda x: x[1])
    
    with open(f'output/{sd}/{save_file}_KPmeaning.txt', 'w') as f: 
        for i in range(len(KP_terms_meaning_p)):
            term = KP_terms_meaning_p[i][0]
            p = KP_terms_meaning_p[i][1]
            f.write(f'{term}\t{p}\n')
                      
    