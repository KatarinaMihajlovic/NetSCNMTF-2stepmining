# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 15:14:45 2024

@author: Katarina
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

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


copy_from = 'step2_NMTF/output'
copy_to = f'{work_dir}/input' 

for root, dirs, files in os.walk(copy_from):
	for file in files:
		if 'Geneslist' in file:
			cell_cond = root.split('\\')[1]
			if not os.path.exists(f'{copy_to}/Geneslist'):
				os.makedirs(f'{copy_to}/Geneslist')
			copyfile(f'{root}/{file}', f'{copy_to}/Geneslist/{file}')

os.chdir(work_dir) 
