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


copy_from = 'step2_NMTF/output'
copy_to = f'{work_dir}/input/NMTF_G1'         

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if 'G1_with_headers' in file:
            cell_cond = root.split('\\')[1]
            if not os.path.exists( f'{copy_to}/{cell_cond}'):
                os.makedirs(f'{copy_to}/{cell_cond}')
            copyfile(f'{root}/{file}', f'{copy_to}/{cell_cond}/{file}')


path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 
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

