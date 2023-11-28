# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 12:50:26 2023

@author: kmihajlo
"""

from shutil import copyfile, rmtree, copytree
import os


work_dir = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')



copy_to = f'{work_dir}/input/LIGERPreds'  
if not os.path.exists(copy_to):
    os.makedirs(copy_to)
       
copy_from = 'PD/output'
for root, dirs, files in os.walk(copy_from):
    for file in files:
        copyfile(f'{root}/{file}', f'{copy_to}/{file}')
        
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

copy_dir = 'step8_StageSpec_CorePreds/output/All_genes'
destination_dir = f'{work_dir}/input/All_genes'
if os.path.exists(destination_dir):
    rmtree(destination_dir) 
copytree(copy_dir, destination_dir) 

copy_from = 'step9_Literature_Validation/output/StageSpecPreds'
copy_to = f'{work_dir}/input/StageSpecPreds'         

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file.endswith('.pkl') and 'OGs' not in file:
            print(file)
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)
            copyfile(f'{root}/{file}', f'{copy_to}/{file}')

copy_to = f'{work_dir}/input'
copy_from = 'step2_NMTF/output'
save_dir = f'{copy_to}/All_genes'
for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file.endswith('.csv') and 'Geneslist' in file:
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)  
            copyfile(f'{root}/{file}', f'{save_dir}/{file}')    

copy_to = f'{work_dir}/input'
copyfile('CommonData/Homo_sapiens.gene_info', f'{copy_to}/Homo_sapiens.gene_info') 
copyfile('CommonData/PD_genes_DGN.pkl', f'{copy_to}/PD_genes_DGN.pkl')             
copyfile('CommonData/LitValid_AllGenes.pkl', f'{copy_to}/LitValid_AllGenes.pkl')

os.chdir(work_dir) 
