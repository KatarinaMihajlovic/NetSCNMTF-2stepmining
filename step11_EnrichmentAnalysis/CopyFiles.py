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

copy_to = f'{work_dir}/input'  
copyfile('step9_Literature_Validation/output/CorePreds/CorePreds_LitValid.pkl', f'{copy_to}/CorePreds_LitValid.pkl') 
copyfile('step9_Literature_Validation/output/All_genes_Union.txt', f'{copy_to}/All_genes_Union.txt')


copy_from = 'CommonData'
copy_to = f'{work_dir}/input'  

copyfile(f'{copy_from}/go-basic.obo', f'{copy_to}/go-basic.obo')
copyfile(f'{copy_from}/Homo_sapiens.gene_info', f'{copy_to}/Homo_sapiens.gene_info')
copyfile(f'{copy_from}/HSA_Kegg_Pathways.lst', f'{copy_to}/HSA_Kegg_Pathways.lst')
copyfile(f'{copy_from}/hsa_pathway_names.lst', f'{copy_to}/hsa_pathway_names.lst')
copyfile(f'{copy_from}/hsa00001.json', f'{copy_to}/hsa00001.json')
copyfile(f'{copy_from}/PD_genes_DGN.pkl', f'{copy_to}/PD_genes_DGN.pkl')
copyfile(f'{copy_from}/PDmap_BasePathways_noSBUK.lst', f'{copy_to}/PDmap_BasePathways_noSBUK.lst') 
copyfile('step10_NetworkAna/output/PINK1_1neighCorePDPreds_PPI.txt', f'{copy_to}/PINK1_1neighCorePDPreds_PPI.txt')


os.chdir(work_dir)
