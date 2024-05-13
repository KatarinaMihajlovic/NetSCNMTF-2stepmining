# -*- coding: UTF-8 -*-

''' Non-negative matrix tri-factorization (numpy)'''
# Author: kmihajlo
import os
os.environ['OMP_NUM_THREADS'] = '8' # export OMP_NUM_THREADS=4
os.environ['OPENBLAS_NUM_THREADS'] = '8' # export OPENBLAS_NUM_THREADS=4 
os.environ['MKL_NUM_THREADS'] = '8' # export MKL_NUM_THREADS=6
os.environ['VECLIB_MAXIMUM_THREADS'] = '8' # export VECLIB_MAXIMUM_THREADS=4
os.environ['NUMEXPR_NUM_THREADS'] = '8' # export NUMEXPR_NUM_THREADS=6


import networkx as nx
import numpy as np
import Network_Matrices as nm
import Matrix_Factorization as mf
import time
import sys
import pandas as pd
import pickle


print('\n\x1b[1;37;44m ############################################# \x1b[0m')
print('\x1b[1;37;44m #                                           # \x1b[0m')
print('\x1b[1;37;44m # ISF Framework                             # \x1b[0m')
print('\x1b[1;37;44m #                                           # \x1b[0m')
print('\x1b[1;37;44m ############################################# \x1b[0m')

case = 'Skupin'
epsilon = 1e-5

#Loading and preparing all the data
MOL_or_EXP = ['MOL','EXP'] # EXP has shape gene*single cells


try:
    os.makedirs('output/EXP_matrices')
except FileExistsError:
    # directory already exists
    pass
    
# Dataset in pairs, PPI + COEX
# nets_1 = ['PPI', 'COEX', 'GI', 'MI']
# nets_2 = ['PPI+COEX', 'PPI+GI', 'COEX+GI', 'PPI+MI', 'COEX+MI', 'GI+MI']
# nets_3 = ['PPI+GI+COEX', 'PPI+GI+MI', 'PPI+COEX+MI', 'GI+COEX+MI']
# nets_4 = ['ALL']
# nets_everything = ['EXP']  + nets_1 + nets_2 + nets_3 + nets_4
nets_everything = ['ALL']

pickle.dump(nets_everything, open('./output/Net_pairs_used.pkl', 'wb'))
'''
k1 = int( math.sqrt(nb_genes/2.) )
k2 = int( math.sqrt(nb_SCs/2.) )

k2s = [25, k2, 75, 100]
'''
cell_cond = str(sys.argv[1])
k1 = int(sys.argv[2])
k2 = int(sys.argv[3])

#ks = [10]
#repetitions = 2

file = 'E_' + cell_cond + '.csv'
print(cell_cond)        

MolNets_dir = 'input/MolecularEXPSpecificNetworks'
print('---- Loading PPI network')
PPI, nodes, node2ind = nm.Load_Network(MolNets_dir + '/' + cell_cond + '/PPI_' + cell_cond + '.edgelist') # nodes = gene nodes, node2ind = gene node to index
PPI_mat = nm.Make_Adj(PPI, nodes, node2ind, MOL_or_EXP[0])

print('---- Loading COEX network')
COEX = nx.read_edgelist(MolNets_dir + '/' + cell_cond + '/COEX_'+ cell_cond + '.edgelist')
COEX_mat = nm.Make_Adj(COEX, nodes, node2ind, MOL_or_EXP[0])

print('---- Loading GI network')
GI = nx.read_edgelist(MolNets_dir + '/' + cell_cond + '/GI_' + cell_cond + '.edgelist')
GI_mat = nm.Make_Adj(GI, nodes, node2ind, MOL_or_EXP[0])

print('---- Loading MI network')
MI = nx.read_edgelist(MolNets_dir + '/' + cell_cond + '/MI_' + cell_cond + '.edgelist')
MI_mat = nm.Make_Adj(MI, nodes, node2ind, MOL_or_EXP[0])

print('---- Loading Expression matrix')
flag = True
if flag == True:
    EXP = pd.read_csv('input/EXP_matrices/' + file, index_col=0)
    EXP_mat, i2SC = nm.Make_Adj(EXP, nodes, node2ind, MOL_or_EXP[1])
    #EXP_mat+=epsilon
    pickle.dump([EXP_mat, i2SC], open('output/EXP_matrices/EXPmat_' + cell_cond + '.pkl', 'wb'))
else:
	EXP_mat, i2SC = pickle.load(open('output/EXP_matrices/EXPmat_' + cell_cond + '.pkl', 'rb'))
#print(EXP_mat)

nb_genes = len(nodes)
nb_SCs = len(i2SC)

SCs = []
for key in i2SC:
   SCs.append(i2SC[key]) 

try:
    os.makedirs('output/' + cell_cond)
except FileExistsError:
    # directory already exists
    pass
    
#start NMTF    

start = time.time()
print(k1, k2)

k1k2_folder = 'output/' + cell_cond + '/k1_' + str(k1) + '_k2_' + str(k2) + '/' 
try:
    os.makedirs(k1k2_folder)
except FileExistsError:
    # directory already exists
    pass

featuresG1 = [str(i) for i in range(k1)]
featuresG2 = [str(i) for i in range(k2)]

for nets in nets_everything:
    try:
        os.makedirs(k1k2_folder + nets)
    except FileExistsError:
        # directory already exists
        pass

#do Data integration of all combinations in nets_everything
EXP = EXP_mat

for nets in nets_everything:
    MOLs = []
    used_nets = []
    if ('PPI' in nets):
        MOLs.append(PPI_mat)
        used_nets.append('PPI')
    if ('COEX' in nets):
        MOLs.append(COEX_mat)
        used_nets.append('COEX')
    if ('GI' in nets):
        MOLs.append(GI_mat)
        used_nets.append('GI')
    if ('MI' in nets):
        MOLs.append(MI_mat)
        used_nets.append('MI')
    if nets == 'ALL': 
        MOLs = [PPI_mat, COEX_mat, GI_mat, MI_mat]  
        used_nets = ['PPI','COEX','GI','MI']
    if nets == 'EXP':
        used_nets = ['noMOLs']
        MOLs = [np.zeros((nb_genes,nb_genes))+epsilon]
    if 'only' in nets:
        EXP = np.zeros((nb_genes,nb_SCs))

    print(nets)
    print(len(MOLs))

    Solver = mf.PD_SSNMTF(max_iter=500, verbose = 10)
    G1, G2, S_Mol, S_EXP = Solver.Solve_MUR(MOLs, EXP, k1, k2, init='SVD')
    # G1, G2, S_Mol, S_EXP = Solver.Solve_MUR(MOLs, EXP, k1, k2, init='rand')

    print('saving matrices')
    np.save(k1k2_folder + nets +'/%s_G1'%(nets), G1)
    #print(G1)
    np.save(k1k2_folder + nets +'/%s_G2'%(nets), G2)
    #print(G2)
    
    '''   
    # pickle save vars for future executions
    print('saving pickle with all vars')
    output_vars = {'G1' : G1, 'G2' : G2, 'nodes': nodes, 'SCs': SCs}
    pickle.dump(output_vars, open(k1k2_folder + nets + '/%s_output_vars.pkl'%nets, 'wb'))
    '''
end = time.time()
run_time = (end - start)/60
print('Runtime of the program is ' + str(run_time))

   
with open('output/' + cell_cond + '/Geneslist_' + cell_cond + '.csv', 'w') as f:
    f.write('genes' + '\n')
    for gene in nodes:
        f.write(gene + '\n')
        
with open('output/' + cell_cond + '/SCslist' + cell_cond + '.csv', 'w') as f:
    f.write('SCs' + '\n')
    for SC in SCs:
        f.write(SC + '\n')
