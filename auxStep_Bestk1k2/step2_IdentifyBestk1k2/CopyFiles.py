from shutil import copyfile, copytree, rmtree
import os

work_dir = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 



if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 

if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')

copy_from = 'step1_ComputeDispersionCoefficient/output'
copy_to = f'{work_dir}/input'

cluster_folders = set()
for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file.endswith('.txt'):
            outdir = root.split('\\')[1]
            if not os.path.exists(f'{copy_to}/{outdir}'):
                os.makedirs(f'{copy_to}/{outdir}')
            copyfile(f'{root}/{file}', f'{copy_to}/{outdir}/{file}')
            
os.chdir(work_dir) 
