# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 12:32:21 2021

@author: kmihajlo
"""

from shutil import copyfile, copytree, rmtree
import os

workdir = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

if not os.path.exists(f'{workdir}/input'):
    os.makedirs(f'{workdir}/input') 
if not os.path.exists(f'{workdir}/output'):
    os.makedirs(f'{workdir}/output')

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

copy_from = 'step4_PDgenesClustEnirch/output/best_runs'
copy_to = f'{workdir}/input'         

for root, dirs, files in os.walk(copy_from):
    for file in files:
        # if '_PDEnrichClusts.pkl' in file and 'ALL' in root and best_clustMeth in file and 'PINK1' in root:
        if 'kMeans' in file:
            if ('PPI+GI+COEX' in root) and ('PINK1' in root) and 'PDEnrich' in file:
                print(root)
                cell_cond = root.split('\\')[1]
                save_dir = f'{copy_to}/Clusters/{cell_cond}'
                if not os.path.exists(save_dir):
                    os.makedirs(save_dir)
                copyfile(f'{root}/{file}', f'{save_dir}/{file}')


copyfile('CommonData/PD_genes_DGN.pkl', f'{copy_to}/PD_genes_DGN.pkl') 

copy_from = 'step2_NMTF/output'
copy_to = f'{workdir}/input/Genelists'         

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if 'Geneslist' in file:
            cell_cond = root.split('\\')[1]
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)
            copyfile(f'{root}/{file}', f'{copy_to}/{file}')

os.chdir(workdir) 
