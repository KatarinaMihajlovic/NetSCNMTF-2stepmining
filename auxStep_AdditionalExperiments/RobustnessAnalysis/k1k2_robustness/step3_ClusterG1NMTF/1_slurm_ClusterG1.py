# -*- coding: UTF-8 -*-

# Author: kmihajlo

import os
import sys
import shutil

tasks_dir = 'tasks'
jobs_dir = 'jobs'
out_dir = 'Nord3_output'

if os.path.exists('output'):
 	shutil.rmtree('output')
 	os.makedirs('output')  
else:
    os.makedirs('output') 

if os.path.exists(tasks_dir):
 	shutil.rmtree(tasks_dir)
 	os.makedirs(tasks_dir)  
else:
    os.makedirs(tasks_dir) 

if os.path.exists(jobs_dir):
 	shutil.rmtree(jobs_dir) 
 	os.makedirs(jobs_dir) 
else:
    os.makedirs(jobs_dir) 

if os.path.exists(out_dir):
 	shutil.rmtree(out_dir) 
 	os.makedirs(out_dir) 
else:
    os.makedirs(out_dir) 


for cell_cond in os.listdir('input/NMTF_G1'):
    task_filename = cell_cond + '.txt'
    task_file = open(tasks_dir + '/' + task_filename, 'w')
    for file in os.listdir(f'input/NMTF_G1/{cell_cond}'):
        k1k2 = file[:-11]
        task_file.write('python ClusterNMTFEmbeddings.py' + ' ' + cell_cond + ' ' + k1k2 + '\n')
    task_file.close()

    job_name = cell_cond  
    job_filename = job_name + '.sh'
    job_file = open(jobs_dir + '/' + job_filename,'w')
		
    job_file.write('#!/bin/bash\n')
    job_file.write('#SBATCH --ntasks=10\n')
    job_file.write('#SBATCH --cpus-per-task=4\n')

    job_file.write('#SBATCH --output=' + out_dir + '/output_%J.out\n')
    job_file.write('#SBATCH --error=' + out_dir + '/Err_output_%J.err\n')
    job_file.write('#SBATCH --job-name=' + job_name + '\n')
    job_file.write('#SBATCH --qos=bsc_ls\n')
    job_file.write('#SBATCH --time=4:00:00\n')
    job_file.write('module purge && module load intel impi greasy mkl python/3.9.10\n')
    job_file.write('FILE=' + tasks_dir + '/' + task_filename + '\n')
    job_file.write('export GREASY_LOGFILE=Nord3_output/nord3_%j.log\n')
    job_file.write('/apps/GREASY/2.2.3/INTEL/IMPI/bin/greasy $FILE\n')
    job_file.close()
    cmd = "sbatch " + jobs_dir + '/' + job_filename 
    os.system(cmd)
