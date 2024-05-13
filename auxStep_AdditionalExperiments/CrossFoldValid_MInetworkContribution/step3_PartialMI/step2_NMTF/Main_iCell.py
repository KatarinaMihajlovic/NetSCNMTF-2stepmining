# -*- coding: UTF-8 -*-

''' Non-negative matrix tri-factorization (numpy)'''
# Author: kmihajlo
import os
os.environ['OMP_NUM_THREADS'] = '4' # export OMP_NUM_THREADS=4
os.environ['OPENBLAS_NUM_THREADS'] = '4' # export OPENBLAS_NUM_THREADS=4 
os.environ['MKL_NUM_THREADS'] = '4' # export MKL_NUM_THREADS=6
os.environ['VECLIB_MAXIMUM_THREADS'] = '4' # export VECLIB_MAXIMUM_THREADS=4
os.environ['NUMEXPR_NUM_THREADS'] = '4' # export NUMEXPR_NUM_THREADS=6


import networkx as nx
import numpy as np
import Network_Matrices as nm
import Matrix_Factorization as mf
import time
import sys
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

print('\n\x1b[1;37;44m ############################################# \x1b[0m')
print('\x1b[1;37;44m #                                           # \x1b[0m')
print('\x1b[1;37;44m # ISF Framework                             # \x1b[0m')
print('\x1b[1;37;44m #                                           # \x1b[0m')
print('\x1b[1;37;44m ############################################# \x1b[0m')

case = 'Skupin'
epsilon = 1e-5

#Loading and preparing all the data
MOL_or_EXP = ['MOL','EXP'] # EXP has shape gene*single cells

if not os.path.exists('output/EXP_matrices'):
    os.makedirs('output/EXP_matrices')

nets_everything = ['ALL']
# nets_everything = ['MI']

ks = {'Control_IPSCs' : [100, 50], 'Control_D06' : [100, 60], 'Control_D15' : [125, 50], 'Control_D21' : [125, 50], 
      'PINK1_IPSCs' : [75, 60], 'PINK1_D06' : [100, 50], 'PINK1_D15' : [100, 60], 'PINK1_D21' : [100, 40]}
legend = {'Control_D06':'C6','Control_D15':'C15', 'Control_D21':'C21', 'Control_IPSCs':'C0',
          'PINK1_D06':'PD6', 'PINK1_D15':'PD15', 'PINK1_D21':'PD21', 'PINK1_IPSCs':'PD0'}



for cc in ks.keys():
    if not os.path.exists(f'output/{cc}'):
        os.makedirs(f'output/{cc}') 
        
    MolNets_dir = 'input/MolecularEXPSpecificNetworks'
    print('---- Loading PPI network')
    PPI, nodes, node2ind = nm.Load_Network(f'{MolNets_dir}/{cc}/PPI_{cc}.edgelist') # nodes = gene nodes, node2ind = gene node to index
    PPI_mat = nm.Make_Adj(PPI, nodes, node2ind, MOL_or_EXP[0])
    
    print('---- Loading COEX network')
    COEX = nx.read_edgelist(f'{MolNets_dir}/{cc}/COEX_{cc}.edgelist')
    COEX_mat = nm.Make_Adj(COEX, nodes, node2ind, MOL_or_EXP[0])
    
    print('---- Loading GI network')
    GI = nx.read_edgelist(f'{MolNets_dir}/{cc}/GI_{cc}.edgelist')
    GI_mat = nm.Make_Adj(GI, nodes, node2ind, MOL_or_EXP[0])
  
    print('---- Loading Expression matrix')
    flag = True
    if flag == True:
        EXP = pd.read_csv(f'input/EXP_matrices/E_{cc}.csv', index_col=0)
        EXP_mat, i2SC = nm.Make_Adj(EXP, nodes, node2ind, MOL_or_EXP[1])
        #EXP_mat+=epsilon
        pickle.dump([EXP_mat, i2SC], open(f'output/EXP_matrices/EXPmat_{cc}.pkl', 'wb'))
    else:
    	EXP_mat, i2SC = pickle.load(open(f'output/EXP_matrices/EXPmat_{cc}.pkl', 'rb'))
    EXP = EXP_mat

    #print(EXP_mat)
    nb_genes = len(nodes)
    nb_SCs = len(i2SC)
    
    SCs = []
    for key in i2SC:
       SCs.append(i2SC[key]) 
       
    k1 = ks[cc][0]
    k2 = ks[cc][1]

    start = time.time()
    print(k1, k2)

    k1k2_folder = f'output/{cc}/k1_{k1}_k2_{k2}' 
    if not os.path.exists(k1k2_folder):
        os.makedirs(k1k2_folder)

    featuresG1 = [str(i) for i in range(k1)]
    featuresG2 = [str(i) for i in range(k2)]

  
    # removing Test genes interactions from MI matrix (5 fold)
    with open(f'input/Folds/{legend[cc]}_Train_Test5FOLD.pkl', 'rb') as handle:
        Train_Test5FOLD = pickle.load(handle)   
    
    # do 5 NMTFs, once for each Fold
    for i in range(len(Train_Test5FOLD)):
        Train = Train_Test5FOLD[i][0]
        Test = Train_Test5FOLD[i][1]
        print('---- Loading MI network')
        MI = nx.read_edgelist(f'{MolNets_dir}/{cc}/MI_{cc}.edgelist')
        MI_mat = nm.Make_Adj(MI, nodes, node2ind, MOL_or_EXP[0])
        
        MI_df =  pd.DataFrame(MI_mat, index=nodes, columns=nodes)
        MI_df.loc[Test, :] = 0
        MI_df.loc[:, Test] = 0
        
        MI_mat = MI_df.values
        
        MOLs = [PPI_mat, COEX_mat, GI_mat, MI_mat]  
        used_nets = ['PPI','COEX','GI','MI']
    
        foldn = f'FOLD{i}'
        #start NMTF    

        Solver = mf.PD_SSNMTF(max_iter=1000, verbose = 10)
        G1, G2, S_Mol, S_EXP, OBJ_fig = Solver.Solve_MUR(MOLs, EXP, k1, k2, init='SVD')
        
        #G1, G2, S_Mol, S_EXP = Solver.Solve_MUR(MOLs, EXP, k1, k2, init='rand')
        OBJ_fig.tight_layout()    
        OBJ_fig.savefig(k1k2_folder + '/' + foldn + '_OBJ.jpg', dpi = 350, format='jpg')
        plt.close(OBJ_fig)
        
        print('saving matrices')
        nm.Save_Matrix_Factor(G1, k1k2_folder  +'/%s_G1_with_headers.csv'%(foldn), nodes, featuresG1)
        nm.Save_Matrix_Factor_no_headers(G1, k1k2_folder  +'/%s_G1.csv'%(foldn), nodes, featuresG1)
        nm.Save_Matrix_Factor(G2, k1k2_folder  +'/%s_G2_with_headers.csv'%(foldn), nodes, featuresG2)
        nm.Save_Matrix_Factor_no_headers(G2, k1k2_folder  +'/%s_G2.csv'%(foldn), nodes, featuresG2)
        
        for i in range(len(used_nets)):
            nm.Save_Matrix_Factor(S_Mol[i], k1k2_folder  + '/'  + foldn + '_S_' + used_nets[i] + '.csv', featuresG1, featuresG1)	
        nm.Save_Matrix_Factor(S_EXP, k1k2_folder  + '/'  + foldn + '_S_EXP.csv', featuresG1, featuresG2)	

    end = time.time()
    run_time = (end - start)/60
    print('Runtime of the program is ' + str(run_time))

   
    with open(f'./output/{cc}/Geneslist_{cc}.csv', 'w') as f:
        f.write('genes' + '\n')
        for gene in nodes:
            f.write(f'{gene}\n')
            
    with open(f'./output/{cc}/SCslist_{cc}.csv', 'w') as f:
        f.write('SCs' + '\n')
        for SC in SCs:
            f.write(f'{SC}\n')
