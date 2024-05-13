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
copy_to = f'{work_dir}/input'         

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if 'G1' in file:
            cell_cond = root.split('\\')[1]
            k1k2 = root.split('\\')[2]
            print(root, file)
            sd = f'{copy_to}/NMTF_G1/{cell_cond}'
            if not os.path.exists(sd):
                os.makedirs(sd)
            copyfile(f'{root}/{file}', f'{sd}/{k1k2}_{file}')
        elif 'Geneslist' in file:
            sd = f'{copy_to}/Geneslist'
            if not os.path.exists(sd):
                os.makedirs(sd)
            copyfile(f'{root}/{file}', f'{sd}/{file}')
            

os.chdir(work_dir) 
