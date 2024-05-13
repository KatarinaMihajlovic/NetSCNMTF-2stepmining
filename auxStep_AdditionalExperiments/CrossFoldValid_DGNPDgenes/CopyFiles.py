# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 10:49:56 2021

@author: kmihajlo
"""


from shutil import copyfile
import os

work_dir = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


copy_from = 'step3_ClusterG1NMTF/output'
copy_to = f'{work_dir}/input/Clusters'         

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if 'PPI+GI+COEX' in root or 'ALL' in root:
            cell_cond = root.split('\\')[1]
            net = root.split('\\')[2]
            if not os.path.exists( f'{copy_to}/{cell_cond}/{net}'):
                os.makedirs(f'{copy_to}/{cell_cond}/{net}')
            copyfile(f'{root}/{file}', f'{copy_to}/{cell_cond}/{net}/{file}')


copy_to = f'{work_dir}/input'         

copyfile('CommonData/PD_genes_DGN.pkl', f'{copy_to}/PD_genes_DGN.pkl') 

os.chdir(work_dir) 