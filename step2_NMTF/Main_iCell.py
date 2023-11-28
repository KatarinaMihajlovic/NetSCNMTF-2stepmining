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

# Dataset in pairs, PPI + COEX
nets_1 = ['PPI','GI','MI','COEX','EXP']
nets_2 = ['PPI+GI','PPI+MI','PPI+COEX','COEX+GI','COEX+MI','GI+MI']
nets_3 = ['GI+COEX+MI','PPI+GI+MI','PPI+COEX+MI','PPI+GI+COEX']
nets_everything = ['ALL'] + nets_1 + nets_2 + nets_3
# nets_everything = ['MI']

ks = {'Control_IPSCs' : [100, 50], 'Control_D06' : [100, 60], 'Control_D15' : [125, 50], 'Control_D21' : [125, 50], 
      'PINK1_IPSCs' : [75, 60], 'PINK1_D06' : [100, 50], 'PINK1_D15' : [100, 60], 'PINK1_D21' : [100, 40]}


for root, dirs, files in os.walk('input'):
    for file in files:
        if file.endswith('.csv'):
            print(file)
            lspt = file.split('.')
            lspt = lspt[0].split('_')
            cell_cond = lspt[1] + '_' + lspt[2]
            print(cell_cond)        
            
            MolNets_dir = 'input/MolecularEXPSpecificNetworks'
            print('---- Loading PPI network')
            PPI, nodes, node2ind = nm.Load_Network(f'{MolNets_dir}/{cell_cond}/PPI_{cell_cond}.edgelist') # nodes = gene nodes, node2ind = gene node to index
            PPI_mat = nm.Make_Adj(PPI, nodes, node2ind, MOL_or_EXP[0])
            
            print('---- Loading COEX network')
            COEX = nx.read_edgelist(f'{MolNets_dir}/{cell_cond}/COEX_{cell_cond}.edgelist')
            COEX_mat = nm.Make_Adj(COEX, nodes, node2ind, MOL_or_EXP[0])
            
            print('---- Loading GI network')
            GI = nx.read_edgelist(f'{MolNets_dir}/{cell_cond}/GI_{cell_cond}.edgelist')
            GI_mat = nm.Make_Adj(GI, nodes, node2ind, MOL_or_EXP[0])
            
            print('---- Loading MI network')
            MI = nx.read_edgelist(f'{MolNets_dir}/{cell_cond}/MI_{cell_cond}.edgelist')
            MI_mat = nm.Make_Adj(MI, nodes, node2ind, MOL_or_EXP[0])
            
            print('---- Loading Expression matrix')
            flag = True
            if flag == True:
                EXP = pd.read_csv(f'input/EXP_matrices/{file}', index_col=0)
                EXP_mat, i2SC = nm.Make_Adj(EXP, nodes, node2ind, MOL_or_EXP[1])
                #EXP_mat+=epsilon
                pickle.dump([EXP_mat, i2SC], open(f'output/EXP_matrices/EXPmat_{cell_cond}.pkl', 'wb'))
            else:
            	EXP_mat, i2SC = pickle.load(open(f'output/EXP_matrices/EXPmat_{cell_cond}.pkl', 'rb'))
            #print(EXP_mat)
            
            nb_genes = len(nodes)
            print(nb_genes)
            nb_SCs = len(i2SC)
            
            SCs = []
            for key in i2SC:
               SCs.append(i2SC[key]) 
            
            if not os.path.exists(f'output/{cell_cond}'):
                os.makedirs(f'output/{cell_cond}') 
                
            #start NMTF    
            k1 = ks[cell_cond][0]
            k2 = ks[cell_cond][1]

            start = time.time()
            print(k1, k2)

            k1k2_folder = f'output/{cell_cond}/k1_{k1}_k2_{k2}' 
            if not os.path.exists(k1k2_folder):
                os.makedirs(k1k2_folder)

            featuresG1 = [str(i) for i in range(k1)]
            featuresG2 = [str(i) for i in range(k2)]

            for nets in nets_everything:
                if not os.path.exists(f'{k1k2_folder}/{nets}'):
                    os.makedirs(f'{k1k2_folder}/{nets}')


            #do Data integration of all combinations in nets_everything
            EXP = EXP_mat

            for nets in nets_everything:
                k1k2_folder = f'output/{cell_cond}/k1_{k1}_k2_{k2}/{nets}' 
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
                
                if nets != 'EXP':
                    Solver = mf.PD_SSNMTF(max_iter=1000, verbose = 10)
                    G1, G2, S_Mol, S_EXP, OBJ_fig = Solver.Solve_MUR(MOLs, EXP, k1, k2, init='SVD')
                else:
                    Solver = mf.NMTF_basic(max_iter=1000, verbose = 10)
                    G1, G2, S_EXP, OBJ_fig = Solver.Solve_MUR(EXP, k1, k2, init='SVD')
                    
                #G1, G2, S_Mol, S_EXP = Solver.Solve_MUR(MOLs, EXP, k1, k2, init='rand')
                OBJ_fig.tight_layout()    
                OBJ_fig.savefig(k1k2_folder + '/' + nets + '_OBJ.jpg', dpi = 350, format='jpg')
                plt.close(OBJ_fig)
                
                print('saving matrices')
                nm.Save_Matrix_Factor(G1, k1k2_folder  +'/%s_G1_with_headers.csv'%(nets), nodes, featuresG1)
                nm.Save_Matrix_Factor_no_headers(G1, k1k2_folder  +'/%s_G1.csv'%(nets), nodes, featuresG1)
                nm.Save_Matrix_Factor(G2, k1k2_folder  +'/%s_G2_with_headers.csv'%(nets), nodes, featuresG2)
                nm.Save_Matrix_Factor_no_headers(G2, k1k2_folder  +'/%s_G2.csv'%(nets), nodes, featuresG2)
                
                if nets != 'EXP':
                    for i in range(len(used_nets)):
                        nm.Save_Matrix_Factor(S_Mol[i], k1k2_folder  + '/'  + nets + '_S_' + used_nets[i] + '.csv', featuresG1, featuresG1)	
                nm.Save_Matrix_Factor(S_EXP, k1k2_folder  + '/'  + nets + '_S_EXP.csv', featuresG1, featuresG2)	

            end = time.time()
            run_time = (end - start)/60
            print('Runtime of the program is ' + str(run_time))

           
            with open(f'./output/{cell_cond}/Geneslist_{cell_cond}.csv', 'w') as f:
                f.write('genes' + '\n')
                for gene in nodes:
                    f.write(f'{gene}\n')
                    
            with open(f'./output/{cell_cond}/SCslist_{cell_cond}.csv', 'w') as f:
                f.write('SCs' + '\n')
                for SC in SCs:
                    f.write(f'{SC}\n')
