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


copy_from = 'step6_InitPreds/output'
copy_to = f'{work_dir}/input/InitPreds'     


for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file.endswith('.pkl'):
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)  
            copyfile(f'{root}/{file}', f'{copy_to}/{file}')    


copy_to = f'{work_dir}/input'                 
dir2copy = 'step7_ComputeGMM/output'
if os.path.exists(f'{copy_to}/SimMeasuresG'): 
    rmtree(f'{copy_to}/SimMeasuresG') 
copytree(dir2copy, f'{copy_to}/SimMeasuresG') 


copyfile('CommonData/PD_genes_DGN.pkl', f'{copy_to}/PD_genes_DGN.pkl') 

os.chdir(work_dir) 
