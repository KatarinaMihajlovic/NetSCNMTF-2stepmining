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

# k1s = [25, 50, 75, 100, 125, 150, 175, 200, 250, 275]
# k2s = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

cell_conds = ['PINK1_D06', 'Control_D06', 'PINK1_IPSCs', 'PINK1_D15', 'Control_D21', 'Control_D15', 'PINK1_D21', 'Control_IPSCs']

k1s_cc = {'Control_IPSCs' : [50, 75, 100, 125, 150], 'Control_D06' : [50, 75, 100, 125, 150,], 
       'Control_D15' : [75, 100, 125, 150, 175], 'Control_D21' : [75, 100, 125, 150, 175], 
      'PINK1_IPSCs' : [25, 50, 75, 100, 125], 'PINK1_D06' : [50, 75, 100, 125, 150], 
      'PINK1_D15' : [50, 75, 100, 125, 150], 'PINK1_D21' : [50, 75, 100, 125, 150]}
k2s_cc = {'Control_IPSCs' : [30, 40, 50, 60, 70], 'Control_D06' : [40, 50, 60, 70, 80], 
       'Control_D15' : [30, 40, 50, 60, 70], 'Control_D21' : [30, 40, 50, 60, 70], 
      'PINK1_IPSCs' : [40, 50, 60, 70, 80], 'PINK1_D06' : [30, 40, 50, 60, 70], 
      'PINK1_D15' : [40, 50, 60, 70, 80], 'PINK1_D21' : [20, 30, 40, 50, 60]}


# k1s = [275]
# k2s = [50]
# cell_conds = ['Control_IPSCs']

for cell_cond in cell_conds:
    task_filename = cell_cond + '.txt'
    task_file = open(tasks_dir + '/' + task_filename, 'w')
    for k1 in k1s_cc[cell_cond]:
        for k2 in k2s_cc[cell_cond]:
            task_file.write('python Main_iCell.py' + ' ' + cell_cond + ' ' + str(k1) + ' ' + str(k2) + '\n')
    task_file.close()

    job_name = cell_cond  
    job_filename = job_name + '.sh'
    job_file = open(jobs_dir + '/' + job_filename,'w')
		
    job_file.write('#!/bin/bash\n')
    job_file.write('#SBATCH --ntasks=5\n')
    job_file.write('#SBATCH --cpus-per-task=8\n')

    job_file.write('#SBATCH --output=' + out_dir + '/output_%J.out\n')
    job_file.write('#SBATCH --error=' + out_dir + '/Err_output_%J.err\n')
    job_file.write('#SBATCH --job-name=' + job_name + '\n')
    job_file.write('#SBATCH --qos=bsc_ls\n')
    job_file.write('#SBATCH --time=16:00:00\n')
    job_file.write('module purge && module load intel impi greasy mkl python/3.9.10\n')
    job_file.write('FILE=' + tasks_dir + '/' + task_filename + '\n')
    job_file.write('export GREASY_LOGFILE=Nord3_output/nord3_%j.log\n')
    job_file.write('/apps/GREASY/2.2.3/INTEL/IMPI/bin/greasy $FILE\n')
    job_file.close()
    cmd = "sbatch " + jobs_dir + '/' + job_filename 
    os.system(cmd)
