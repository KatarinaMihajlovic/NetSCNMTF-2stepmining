# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 16:57:53 2021

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

prev_step = 'step8_StageSpec_CorePreds/output'
copy_from = f'{prev_step}/StageSpecPreds'
copy_to = f'{workdir}/input'         


for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file.endswith('.pkl'):
            print(root, file)
            sd = copy_from.split('/')[2]
            save_dir = f'{copy_to}/{sd}'
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)             
            copyfile(f'{root}/{file}', f'{save_dir}/{file}') 

copy_from = f'{prev_step}/StageSpecPreds/EuclideanDistance_kMeans/Eucl_dist'
copy_to = f'{workdir}/input'   
for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file.endswith('.csv') and '' in root:
            print(root, file)
            sd = copy_from.split('/')[2] + '/' + copy_from.split('/')[3] + '_' + copy_from.split('/')[4]
            save_dir = f'{copy_to}/{sd}'
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)             
            copyfile(f'{root}/{file}', f'{save_dir}/{file}')    
            


copy_dir = f'{prev_step}/All_genes'
destination_dir = f'{copy_to}/All_genes'
if os.path.exists(destination_dir):
    rmtree(destination_dir) 
copytree(copy_dir, destination_dir)

copy_from = 'step2_NMTF/output'
save_dir = f'{copy_to}/All_genes'
for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file.endswith('.csv') and 'Geneslist' in file:
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)  
            copyfile(f'{root}/{file}', f'{save_dir}/{file}')                       

copyfile(f'{prev_step}/CorePreds.pkl', f'{copy_to}/CorePreds.pkl') 
copyfile('CommonData/Homo_sapiens.gene_info', f'{copy_to}/Homo_sapiens.gene_info') 

copyfile('CommonData/PD_genes_DGN.pkl', f'{copy_to}/PD_genes_DGN.pkl')             
copyfile('CommonData/LitValid_AllGenes.pkl', f'{copy_to}/LitValid_AllGenes.pkl')
copyfile('CommonData/DEGs_Skupin.txt', f'{copy_to}/DEGs_Skupin.txt')


os.chdir(workdir)
