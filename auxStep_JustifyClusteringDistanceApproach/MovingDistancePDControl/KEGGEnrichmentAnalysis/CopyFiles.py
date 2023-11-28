# -*- coding: utf-8 -*-
"""
Created on Thu Oct  14 12:48:56 2021

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

copy_to = f'{work_dir}/input' 
copy_from = 'output' 
for root, dirs, files in os.walk(copy_from):
    for file in files:
        if 'topnCGs' in file:
            copyfile(f'{root}/{file}', f'{copy_to}/{file}')    

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

copy_dir = 'step9_Literature_Validation/input/All_genes'
destination_dir = f'{copy_to}/All_genes'
if os.path.exists(destination_dir):
    rmtree(destination_dir) 
copytree(copy_dir, destination_dir)


# copy_from = 'CommonData'
# copy_to = f'{work_dir}/input'  

# for root, dirs, files in os.walk(copy_from):
#     for file in files:
#         if not file.endswith('.txt'):
#             copyfile(f'{copy_from}/{file}', f'{copy_to}/{file}')    
copyfile('CommonData/hsa_pathway_names.lst', f'{copy_to}/hsa_pathway_names.lst')
copyfile('CommonData/HSA_Kegg_Pathways.lst', f'{copy_to}/HSA_Kegg_Pathways.lst')
copyfile('CommonData/Homo_sapiens.gene_info', f'{copy_to}/Homo_sapiens.gene_info') 


os.chdir(work_dir)
