# -*- coding: utf-8 -*-
"""
Created on Mon May 13 15:53:30 2024

@author: Katarina
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

copy_from = 'step2_BioAnnotsClustEnrich_MI_noMI/output'
copy_to = f'{work_dir}/input'

copyfile(f'{copy_from}/ALL_BPperc_TestCorrectClst.npy', f'{copy_to}/ALL_BPperc_TestCorrectClst.npy')
copyfile(f'{copy_from}/ALL_KPperc_TestCorrectClst.npy', f'{copy_to}/ALL_KPperc_TestCorrectClst.npy')
copyfile(f'{copy_from}/ALL_RPperc_TestCorrectClst.npy', f'{copy_to}/ALL_RPperc_TestCorrectClst.npy')

copyfile(f'{copy_from}/PPI+GI+COEX_BPperc_TestCorrectClst.npy', f'{copy_to}/PPI+GI+COEX_BPperc_TestCorrectClst.npy')
copyfile(f'{copy_from}/PPI+GI+COEX_KPperc_TestCorrectClst.npy', f'{copy_to}/PPI+GI+COEX_KPperc_TestCorrectClst.npy')
copyfile(f'{copy_from}/PPI+GI+COEX_RPperc_TestCorrectClst.npy', f'{copy_to}/PPI+GI+COEX_RPperc_TestCorrectClst.npy')

copy_from = 'step3_PartialMI/step4_BioAnnotsClustEnrich/output'
copyfile(f'{copy_from}/PartMImissing_BPperc_TestCorrectClst.npy', f'{copy_to}/PartMImissing_BPperc_TestCorrectClst.npy')
copyfile(f'{copy_from}/PartMImissing_KPperc_TestCorrectClst.npy', f'{copy_to}/PartMImissing_KPperc_TestCorrectClst.npy')
copyfile(f'{copy_from}/PartMImissing_RPperc_TestCorrectClst.npy', f'{copy_to}/PartMImissing_RPperc_TestCorrectClst.npy')

os.chdir(work_dir) 
