from shutil import copyfile, copytree, rmtree
import os

work_dir = os.getcwd()

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

MolNets_dir = 'MolecularEXPSpecificNetworks'


if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/input/EXP_matrices'):
    os.makedirs(f'{work_dir}/input/EXP_matrices') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')

if not os.path.exists(f'{work_dir}/input/{MolNets_dir}'):
    os.makedirs(f'{work_dir}/input/{MolNets_dir}') 


copy_from = 'step1_CreateInputNetworksMatrices/output'
copy_to = f'{work_dir}/input'

cluster_folders = set()
for root, dirs, files in os.walk(f'{copy_from}/Expression_Matrix'):
    for file in files:
        if file.endswith('.csv'):
            print(file)
            copyfile(f'{copy_from}/Expression_Matrix/{file}', f'{copy_to}/EXP_matrices/{file}')
            

if os.path.exists(f'{copy_to}/{MolNets_dir}'): 
    rmtree(f'{copy_to}/{MolNets_dir}') 
copytree(f'{copy_from}/{MolNets_dir}', f'{copy_to}/{MolNets_dir}')

os.chdir(work_dir) 
