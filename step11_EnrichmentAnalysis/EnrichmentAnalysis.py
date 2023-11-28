#import networkx as nx
import math, json
from scipy.stats import hypergeom
import numpy as np
import matplotlib.pyplot as plt
import os.path
import pandas as pd
from matplotlib_venn import venn2
from matplotlib.cm import ScalarMappable


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

def GO_meaning(clusters_GOs):
    GOterms_meaning = []
    for clust in clusters_GOs:
        GOterms_meaning_clst = []
        for GO_term in clust:
            try:
                GOterm_meaning = godag[GO_term].name
            except KeyError:
                GOterm_meaning = GO_term
            GOterms_meaning_clst.append(GOterm_meaning)
        GOterms_meaning.append(GOterms_meaning_clst)    
    return GOterms_meaning

def parseKEGGnames(KEGGnames_file):
    KEGG2names = {}
    for line in KEGGnames_file.readlines():
        lspt=line.strip().split('\t')
        KEGG_name = lspt[0].split(':')[1]
        name = lspt[1].split('Homo')[0]
        KEGG2names[KEGG_name] = name
    return(KEGG2names)    

def ReactomeNames(ReactomeNamesfile):
    Reactome2names = {}
    for line in ReactomeNamesfile.readlines():
        lspt=line.strip().split('\t')
        Reactome_name = lspt[0]
        name = lspt[1]
        Reactome2names[Reactome_name] = name
    return(Reactome2names) 

def SaveSharedORUniqueTerms(sd, termsmeaning_2sets, termtype, case = 'Shared'):
    if not os.path.exists(sd):
        os.makedirs(sd) 
    if case == 'Shared':
        terms = set.intersection(*map(set,termsmeaning_2sets))
        with open(f'{sd}/CP_PDG_{termtype}.txt', 'w') as f: 
            for term in terms:
                f.write(f'{term}\n')  
    elif case == 'Unique':
        terms = list(set(termsmeaning_2sets[0]) - set(termsmeaning_2sets[1]))  
        genes = sd.split('/')[1]
        with open(f'{sd}/{genes}_{termtype}Unique.txt', 'w') as f: 
            for term in terms:
                f.write(f'{term}\n')

def KEGGSubCats(sd, KP_termsmeaning_set):     
    PDs_KEGGSubGroups = {}
    for KEGG_path in KP_termsmeaning_set:
        try:
            Subgroup = KEGGPathwCats_df[KEGGPathwCats_df['PathwayName']==KEGG_path]['Subgroup'].values[0]   
            if Subgroup in PDs_KEGGSubGroups:
                PDs_KEGGSubGroups[Subgroup].append(KEGG_path)
            else:
                PDs_KEGGSubGroups[Subgroup] = [KEGG_path]   
        except IndexError:
           pass
    name = sd.split('/')[1]
    with open(f'{sd}/{name}_KEGGSubGroups.txt', 'w') as f:
        for key in PDs_KEGGSubGroups:
            f.write(f'{key}\t{len(PDs_KEGGSubGroups[key])}\t{PDs_KEGGSubGroups[key]}\n')

        
def VennDiagramTermss(terms_2sets, case): 
    sd = 'output/CP_PDG'       
    CGs_PDGs_inters= set.intersection(*map(set,terms_2sets))
    fig, ax = plt.subplots(figsize=(16, 9))
    len1 = len(terms_2sets[0]) - len(CGs_PDGs_inters)
    len2 = len(terms_2sets[1]) - len(CGs_PDGs_inters)
    len_inters =  len(CGs_PDGs_inters)
    out1 = venn2(subsets = (len1, len2, len_inters), alpha = 0.6, set_labels = ('Core PD Preds', 'PD Genes'), ax = ax)#, set_labels = (label1, label2))
    for text in out1.set_labels:  #change label size
     text.set_fontsize(26);
    for text in out1.subset_labels:  #change number size
     text.set_fontsize(24)
    plt.title(f'Enriched {case}',fontsize=28)#,pad=20)
    plt.savefig(f'{sd}/CP_PDG_{case}', dpi =600)
    plt.show()
    plt.close()

def Sort(sub_li): 
    # reverse = None (Sorts in Ascending order)
    # key is set to sort using second element of 
    # sublist lambda has been used
    sub_li.sort(key = lambda x: len(x[1]), reverse = True)
    return sub_li

def Sort_val(sub_li): 
    # reverse = None (Sorts in Ascending order)
    # key is set to sort using second element of 
    # sublist lambda has been used
    sub_li.sort(key = lambda x: x[1])
    return sub_li

def write_txt_enrRes(enrRes_dict, sd, save_file, save_ext):
    enr_terms = []
    with open(f'{sd}/{save_file}_{save_ext}.txt', 'w') as f: 
        for i in enrRes_dict:
            for pair in enrRes_dict[i]:                        
                f.write(f'{pair[0]}\t{pair[1]}\n')
                enr_terms.append(pair[0])
            f.write('\n')
    return enr_terms
 
    
"""
	Main code starts here
"""
from goatools.obo_parser import GODag
import pickle

indir = 'input'
godag = GODag(obo_file="input/go-basic.obo") 
Entrez2Symbol_file = open('input/Homo_sapiens.gene_info', 'r')
Entrez2Sym = Entrez2Symbols(Entrez2Symbol_file)   
KEGGnames_file = open('input/hsa_pathway_names.lst', 'r')
KEGG2names = parseKEGGnames(KEGGnames_file)
KEGGPathwCats_df = parseKEGGCats('input/hsa00001.json', KEGG2names)


with open('input/All_genes_Union.txt') as f:
    All_genes_Union = f.read().splitlines() 
       
        
fgoname_kp = f'{indir}/HSA_Kegg_Pathways.lst'
go_counts_kp,gene2go_kp = load_GO(fgoname_kp, All_genes_Union)
print("%i KP annotated genes"%(len(gene2go_kp)))

 
KP_terms_2sets = []      
KP_termsmeaning_2sets = []    
     
for root, dirs, files in os.walk(indir):
    for file in files:
        if 'CorePreds' in file or 'PD_genes_DGN' in file:
            print(file)
            if 'CorePreds' in file:           
                with open(f'{root}/{file}', 'rb') as handle:
                    PD_genes = pickle.load(handle)  
                PD_genes = list(PD_genes.keys())
                PD_preds_common = PD_genes
                save_file = sd = 'CorePreds'
            elif 'DGN' in file:
                with open(f'{root}/{file}', 'rb') as handel:
                    PD_genes = pickle.load(handel)
                PD_known = PD_genes
                save_file = sd = 'PDgenes'
                
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
                    pass
            # print(KP_terms)
            # print((KP_terms[0]))
            KP_terms_2sets.append(KP_terms[0])
            KP_termsmeaning_2sets.append(KP_terms_meaning)
        
            KP_terms_meaning_p =  [list(i) for i in zip(KP_terms_meaning, KP_terms_p)]
            KP_terms_meaning_p.sort(key = lambda x: x[1])
                
            with open(f'output/{sd}/{save_file}_KPmeaning.txt', 'w') as f: 
                for i in range(len(KP_terms_meaning_p)):
                    term = KP_terms_meaning_p[i][0]
                    p = KP_terms_meaning_p[i][1]
                    f.write(f'{term}\t{p}\n')
                              
            
            # save PD preds that are annotated with enriched Kegg pathways; save Pathways and Core PD preds that annotate them 
            if 'CorePreds' in file: 
                terms_genes = []
                terms_genes_count_p = []

                for term in KP_terms[0]:
                    term_genes = []
                    # print(term, KEGG2names[term])
                    # f.write(f'{KEGG2names[term]}\n')   
                    for gene, KP_terms_d in gene2go_kp.items():  # for name, age in dictionary.iteritems():  (for Python 2.x)
                        if term in KP_terms_d and gene in PD_preds_common:
                            term_genes.append(gene)                                # f.write(f'{gene}\n')   
                    terms_genes.append([term, term_genes])
                    for term_p in KP_terms_meaning_p:
                        # print(term_p, term)
                        try:
                            if KEGG2names[term] == term_p[0]:
                                terms_genes_count_p.append([KEGG2names[term][:-3], term_p[1], len(term_genes)])
                        except Exception:
                            pass

                                    

### Shared Terms
sd = 'output/CP_PDG'

SaveSharedORUniqueTerms(sd, KP_termsmeaning_2sets, 'KPmeaning')
# VennDiagramTermss(KP_terms_2sets, 'KEGG Pathways')

### Core PD preds
sd = 'output/CorePreds'
SaveSharedORUniqueTerms(sd, KP_termsmeaning_2sets, 'KPmeaning', case = 'Unique')
CorePreds_KEGGSubGroups = KEGGSubCats(sd, KP_termsmeaning_2sets[0])

### PD genes  
sd = 'output/PDgenes'
SaveSharedORUniqueTerms(sd, KP_termsmeaning_2sets, 'KPmeaning', case = 'Unique')
PDGenes_KEGGSubGroups = KEGGSubCats(sd, KP_termsmeaning_2sets[1])



# plot enrichments in a horizontal color bar

terms_genes_count_p.sort(key=lambda x: x[1], reverse = False)
terms_genes_count_p_t10 = terms_genes_count_p[:10]
y_pos = np.arange(len(terms_genes_count_p_t10))
data_color = [x[1] for x in terms_genes_count_p_t10]
data_color.reverse()
Keggs = [x[0] for x in terms_genes_count_p_t10]
Keggs.reverse()
GC = [x[2] for x in terms_genes_count_p_t10]
GC.reverse()
  

fig, ax = plt.subplots(figsize=(13, 9))
ax.set_yticks(y_pos, labels=Keggs, fontsize = 28)
ax.set_xlabel('Gene Count', fontsize = 32)
plt.xticks(fontsize=22)
ax.set_title('Top 10 most significantly enriched KEGG pathways', fontsize=34, fontweight = 'bold', pad=25, x=-0.3, y=1.02)
my_cmap = plt.colormaps.get_cmap("plasma")
rescale = lambda data_color: (data_color - np.min(data_color)) / (np.max(data_color) - np.min(data_color))
colors = my_cmap(rescale(data_color))
ax.barh(y_pos, GC, color=colors, align='center')
sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(min(data_color),max(data_color)))
sm.set_array([])
# fig.colorbar(sm)
cbaxes = fig.add_axes([1.05, 0.1, 0.03, 0.8])  
cbar = plt.colorbar(sm,ax = [ax], cax = cbaxes)
cbar.set_label('p.adjust', rotation=270,labelpad=30, fontsize = 32)
cbar.ax.tick_params(labelsize=22)
cbar.ax.yaxis.offsetText.set(size=22)
plt.tight_layout(pad=1)
plt.savefig('output/CorePreds/EnrichmentKEGG_hbar_CorePreds.jpg', dpi = 350, bbox_inches='tight', format='jpg')
plt.show()
plt.close()



### enrichment in PD map pathways
fgoname = 'input/PDmap_BasePathways_noSBUK.lst'
geneset=set(All_genes_Union)
PDmap_counter = {}
gene2PDmap = {}
#loading raw gene & go relationships
fRead = open(fgoname, "r")
for line in fRead.readlines():
	lspt = line.strip().split('\t')
	if len(lspt)>1:
		geneid = lspt[0]
		term = lspt[1]
		if geneid in geneset:
			if term not in PDmap_counter:
				PDmap_counter[term] = 0
			PDmap_counter[term] += 1
			if geneid not in gene2PDmap:
				gene2PDmap[geneid] = set()
			gene2PDmap[geneid].add(term)
fRead.close()
PDmap_counter2 = {}
gene2PDmap2 = {}
removed = set()
for term in PDmap_counter:
	if PDmap_counter[term] > 1:
		PDmap_counter2[term] = PDmap_counter[term]
	else:
		removed.add(term)
for gene in gene2PDmap:
	genego = gene2PDmap[gene].difference(removed)
	if len(genego)>0:
		gene2PDmap2[gene]=genego
        
        
file = 'PINK1_1neighCorePDPreds_PPI.txt'      
with open(f'input/{file}') as f:
    PD_genes = f.read().splitlines() 
    try:
        PD_genes.remove('')
    except Exception:
        pass        
save_file = sd = file[:-4]
sd = f'output/{sd}'
if not os.path.exists(sd):
        os.makedirs(sd) 
        
# PD map
bh_PDmap, mne_PDmap, eg_PDmap, perc_cluster_PDmap, enriched_genes_PDmap = go_enrichment([PD_genes], gene2PDmap2, PDmap_counter2)
# bh_PDmap = bh_PDmap[0]
bh_PDmap = [Sort_val(l) for l in bh_PDmap]    
EnrClsts_PDmap = {v: k for v, k in enumerate(bh_PDmap)}
enr_PDmap_terms = write_txt_enrRes(EnrClsts_PDmap, sd, save_file, save_ext='PDmap')        


        