import pandas as pd
import os
import networkx as nx
import numpy as np 


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


def CreateNetwork(GeneralNet, genes_SC, genes_PPI, gen_nets_dir, Entrez2Sym = Entrez2Sym):
    print(GeneralNet)
    if 'MI_General_KEGG' in GeneralNet:
        gennet = nx.read_edgelist(f'{gen_nets_dir}/{GeneralNet}', delimiter='\t')
    else:
        gennet = nx.read_edgelist(f'{gen_nets_dir}/{GeneralNet}')

    gennet = nx.relabel_nodes(gennet, Entrez2Sym)
    gennet = gennet.subgraph(genes_PPI)
    gennet = gennet.subgraph(genes_SC)
    

    print(GeneralNets[i].split('_')[0])
    return gennet

def cell_conds2singlecells(Metadata_file):
    cell_conds_2_singlecells = {}
    Metadata_file.readline()
    cell_conds = set()
    for line in Metadata_file.readlines():
        lspt = line.strip().split(',')
        key = lspt[1][1:-1]
        if key in cell_conds_2_singlecells:
            cell_conds_2_singlecells[key].append(lspt[0][1:-1])
        else:
            cell_conds_2_singlecells[key] = [lspt[0][1:-1]]
        cell_conds.add(lspt[1][1:-1])
    cell_conds = list(cell_conds)
    return(cell_conds_2_singlecells, cell_conds)

if not os.path.exists('output/MolecularEXPSpecificNetworks'):
    os.makedirs('output/MolecularEXPSpecificNetworks') 
if not os.path.exists('output/Expression_Matrix'):
    os.makedirs('output/Expression_Matrix') 
    

basefolderIN = 'input/'
basefolderOUT = 'output/'   


EXPs = pd.read_csv(f'{basefolderIN}/ExpressionMatrix/Normalized_data.zip', index_col=0)
genes_SC = EXPs.index.tolist()

Metadata_file = open(f'{basefolderIN}/ExpressionMatrix/Metadata.csv','r')
cell_conds_2_singlecells, cell_conds = cell_conds2singlecells(Metadata_file)

gen_nets_dir = 'input/GeneralNetworks'
GeneralNets = os.listdir(gen_nets_dir)
PPI = nx.read_edgelist(f'{gen_nets_dir}/PPI_General_Biogrid.edgelist')
PPI = nx.relabel_nodes(PPI, Entrez2Sym)
nx.write_edgelist(PPI, 'output/PPI_General_Biogrid_GeneSym.edgelist', data=False)



for cell_cond in cell_conds:
    if cell_cond != 'Control_D10':
        print(cell_cond)

        if not os.path.exists(f'output/MolecularEXPSpecificNetworks/{cell_cond}'):
            os.makedirs(f'output/MolecularEXPSpecificNetworks/{cell_cond}') 

        savedir_nets = f'{basefolderOUT}/MolecularEXPSpecificNetworks/{cell_cond}'
        
        EXP_cellcond = pd.DataFrame(index=EXPs.index)
        SCs = []
        data = []
        for SC in EXPs:
            if SC in cell_conds_2_singlecells[cell_cond]:
                SCs.append(SC)
                data.append(list(EXPs[SC]))
        data = np.array(data)
        data = data.T        
        EXP_cellcond = pd.DataFrame(data, index=EXPs.index, columns = SCs)
        EXP_values = EXP_cellcond.values.tolist()

        zero_genes = []
        for i in range(len(EXP_values)):
            row = EXP_values[i]
            sum_row = sum(row)
            if sum_row == 0:
                zero_genes.append(genes_SC[i])
        expressed_genes = [x for x in genes_SC if x not in zero_genes]   
        EXP_cellcond_unique = EXP_cellcond.loc[expressed_genes, :]

        PPI = nx.read_edgelist(f'{gen_nets_dir}/PPI_General_Biogrid.edgelist')
        PPI = nx.relabel_nodes(PPI, Entrez2Sym)
        genes_PPI = list(PPI.nodes())

        MolecularEXPSpecificNetworks_info = open(f'{savedir_nets}/MolecularEXPSpecificNetworks_Statistics_{cell_cond}.txt','w')
        
        for i in range(len(GeneralNets)):            
            SpecificNet = CreateNetwork(GeneralNets[i], expressed_genes, genes_PPI, gen_nets_dir)         
            net = GeneralNets[i].split('_')[0]
            net_name =  f'{net}_{cell_cond}.edgelist'
            nx.write_edgelist(SpecificNet, f'{savedir_nets}/{net_name}', data=False)
            SpecificNet = nx.read_edgelist(f'{savedir_nets}/{net_name}')
            print(f"genes: {SpecificNet.number_of_nodes()} \n interactions: {SpecificNet.number_of_edges()}\n density: {nx.density(SpecificNet)}\n")

            MolecularEXPSpecificNetworks_info.write(f'{net} \n genes: {SpecificNet.number_of_nodes()} \n interactions: {SpecificNet.number_of_edges()}\n density: {nx.density(SpecificNet)}\n\n')
        MolecularEXPSpecificNetworks_info.close()   
        
        PPI = nx.read_edgelist(f'{savedir_nets}/PPI_{cell_cond}.edgelist')
        # PPI = nx.relabel_nodes(PPI, Entrez2Sym)
        genes_PPI = list(PPI.nodes())
        EXP_cellcond_unique.to_csv(f'{basefolderOUT}/Expression_Matrix/E_{cell_cond}.csv' , header = True, index = True)
        
