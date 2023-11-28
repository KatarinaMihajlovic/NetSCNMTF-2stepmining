# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 10:49:56 2021

@author: kmihajlo
"""

from shutil import copyfile, copytree, rmtree
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

copy_from = 'step2_NMTF/output'
copy_to = f'{work_dir}/input'         

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file == 'ALL_G1_with_headers.csv':
            cell_cond = root.split('/')[2]
            if not os.path.exists( f'{copy_to}/NMTF_G1s/{cell_cond}'):
                os.makedirs(f'{copy_to}/NMTF_G1s/{cell_cond}')
            copyfile(f'{root}/{file}', f'{copy_to}/NMTF_G1s/{cell_cond}/{file}')

copy_from = 'step9_Literature_Validation/output/StageSpecPreds'
copy_to = f'{work_dir}/input/StageSpecPreds'         

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file.endswith('.pkl') and 'Unique' not in file:
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)
            copyfile(f'{root}/{file}', f'{copy_to}/{file}')

copy_from = 'step9_Literature_Validation/output/StageSpecPreds/Other_genes'

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file.endswith('.pkl'):
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)
            copyfile(f'{root}/{file}', f'{copy_to}/{file}')


copy_to = f'{work_dir}/input'  
copyfile('CommonData/PD_genes_DGN.pkl', f'{copy_to}/PD_genes_DGN.pkl')
copyfile('CommonData/Homo_sapiens.gene_info', f'{copy_to}/Homo_sapiens.gene_info') 
copyfile('CommonData/LitValid_AllGenes.pkl', f'{copy_to}/LitValid_AllGenes.pkl') 
copyfile('step9_Literature_Validation/output/CorePreds/All_CommonGenes_LitValid.pkl', f'{copy_to}/All_CommonGenes_LitValid.pkl') 

os.chdir(work_dir) 
