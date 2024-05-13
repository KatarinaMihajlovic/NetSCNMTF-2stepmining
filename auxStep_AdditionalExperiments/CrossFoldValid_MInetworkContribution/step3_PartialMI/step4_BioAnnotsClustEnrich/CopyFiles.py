# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 15:29:39 2024

@author: Katarina
"""


from shutil import copyfile, copytree, rmtree
import os

work_dir = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

legend = {'Control_D06':'C6','Control_D15':'C15', 'Control_D21':'C21', 'Control_IPSCs':'C0',
          'PINK1_D06':'PD6', 'PINK1_D15':'PD15', 'PINK1_D21':'PD21', 'PINK1_IPSCs':'PD0'}


if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')

copy_from = 'step2_NMTF/output'
copy_to = f'{work_dir}/input' 

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if 'Geneslist' in file:
            sd = f'{copy_to}/Geneslist'
            if not os.path.exists(sd):
                os.makedirs(sd)
            copyfile(f'{root}/{file}', f'{sd}/{file}')


dir2copy = 'step3_ClusterG1NMTF/output'
copy_to = f'{work_dir}/input/Clusters'         

if os.path.exists(copy_to): 
    rmtree(copy_to) 
copytree(dir2copy, copy_to) 


path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 
copy_to = f'{work_dir}/input/Folds'
dir2copy = 'Folds/output'
if os.path.exists(copy_to): 
    rmtree(copy_to) 
copytree(dir2copy, copy_to) 


path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 



copy_to = f'{work_dir}/input'
copyfile('CommonData/HSA_GO-BP.lst', f'{copy_to}/HSA_GO-BP.lst')
copyfile('CommonData/HSA_Kegg_Pathways.lst', f'{copy_to}/HSA_Kegg_Pathways.lst')
copyfile('CommonData/HSA_Reactome_Pathways.lst', f'{copy_to}/HSA_Reactome_Pathways.lst')
copyfile('CommonData/Homo_sapiens.gene_info', f'{copy_to}/Homo_sapiens.gene_info')

os.chdir(work_dir) 
