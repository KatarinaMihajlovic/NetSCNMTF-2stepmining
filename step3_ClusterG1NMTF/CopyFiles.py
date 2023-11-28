# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 10:49:56 2021

@author: kmihajlo
"""
from shutil import copyfile
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
cell_conds = os.listdir('step2_NMTF/input/MolecularEXPSpecificNetworks')

for cell_cond in cell_conds:
    if not os.path.exists( f'{copy_to}/{cell_cond}'):
        os.makedirs(f'{copy_to}/{cell_cond}')
    for root, dirs, files in os.walk(f'{copy_from}/{cell_cond}'):
        for file in files:
            if 'G1_with_headers' in file:
                copyfile(f'{root}/{file}', f'{copy_to}/{cell_cond}/{file}')


os.chdir(work_dir) 
