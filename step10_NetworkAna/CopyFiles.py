# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 14:12:56 2023

@author: kmihajlo
"""

from shutil import copyfile, copytree, rmtree
import os
import pandas as pd

work_dir = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')


copy_from = 'step1_CreateInputNetworksMatrices/output'
copy_to = f'{work_dir}/input'         
copyfile(f'{copy_from}/PPI_General_Biogrid_GeneSym.edgelist', f'{copy_to}/PPI_General_Biogrid_GeneSym.edgelist')

copy_from = 'step8_StageSpec_CorePreds/output'
copyfile(f'{copy_from}/CorePreds.txt', f'{copy_to}/CorePreds.txt')

copyfile('step9_Literature_Validation/output/All_genes_Union.txt', f'{copy_to}/All_genes_Union.txt')
copyfile('CommonData/DEGs_Skupin.txt', f'{copy_to}/DEGs_Skupin.txt')

os.chdir(work_dir) 
