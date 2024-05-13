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

MolNets_dir = 'MolecularEXPSpecificNetworks'


if not os.path.exists(f'{work_dir}/input/EXP_matrices'):
    os.makedirs(f'{work_dir}/input/EXP_matrices') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')

if not os.path.exists(f'{work_dir}/input/{MolNets_dir}'):
    os.makedirs(f'{work_dir}/input/{MolNets_dir}') 


copy_from = 'PD-Genes-main/step1_CreateInputNetworksMatrices/output'
copy_to = f'{work_dir}/input'

for root, dirs, files in os.walk(f'{copy_from}/Expression_Matrix'):
    for file in files:
        if file.endswith('.csv'):
            print(file)
            copyfile(f'{copy_from}/Expression_Matrix/{file}', f'{copy_to}/EXP_matrices/{file}')
            

if os.path.exists(f'{copy_to}/{MolNets_dir}'): 
    rmtree(f'{copy_to}/{MolNets_dir}') 
copytree(f'{copy_from}/{MolNets_dir}', f'{copy_to}/{MolNets_dir}')


os.chdir(work_dir) 
