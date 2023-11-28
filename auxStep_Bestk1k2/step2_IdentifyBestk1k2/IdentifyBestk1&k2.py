# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 13:58:07 2021

@author: kmihajlo
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
from kneed import KneeLocator


in_dir = 'input'
out_dir = 'output'
cell_conds = os.listdir(in_dir)
columns = ['Net_comb', 'k1', 'k2', 'dispersionCoefG1', 'dispersionCoefG2', 'DispCoeff_avg']

labels_key = {'Control_IPSCs':'Control$_{D0}$', 'Control_D06':'Control$_{D6}$', 'Control_D15':'Control$_{D15}$', 'Control_D21':'Control$_{D21}$', 
              'PINK1_IPSCs':'PD$_{D0}$', 'PINK1_D06':'PD$_{D6}$', 'PINK1_D15':'PD$_{D15}$', 'PINK1_D21':'PD$_{D21}$'}
bestk1 = {'Control_IPSCs':100, 'Control_D06':125, 'Control_D15':125, 'Control_D21':125, 
              'PINK1_IPSCs':75, 'PINK1_D06':100, 'PINK1_D15':100, 'PINK1_D21':100}
bestk2 = {'Control_IPSCs':50, 'Control_D06':60, 'Control_D15':50, 'Control_D21':50, 
              'PINK1_IPSCs':60, 'PINK1_D06':50, 'PINK1_D15':60, 'PINK1_D21':40}

cm = 'kmeans'
            
for cell_cond in cell_conds:
    print(cell_cond)
    DispCoeff_df = pd.DataFrame(columns=columns)
    cnt = 0
    for file in os.listdir(f'{in_dir}/{cell_cond}'):
        #print(file)
        if cm in file:
            with open(f'{in_dir}/{cell_cond}/{file}') as f:
                lines = f.readlines()
                if lines:
                    lines.pop(0)
                    for line in lines:
                        lspt = line.split('\t')
                        Net_comb = lspt[0]
                        k1 = lspt[1]
                        k2 = lspt[2]
                        dispersionCoefG1 = float(lspt[3])
                        dispersionCoefG2 = float(lspt[4])
                        DispCoeff_avg = (float(dispersionCoefG1) + float(dispersionCoefG2))/2
                        DispCoeff_df.loc[cnt] = [Net_comb, k1, k2, dispersionCoefG1, dispersionCoefG2, DispCoeff_avg] 
                        cnt+=1
    
    DispCoeff_df = DispCoeff_df.sort_values(by = 'DispCoeff_avg', ascending=False)  
    # print(DispCoeff_df)      
    if not os.path.exists(f'{out_dir}/{cell_cond}'):
        os.makedirs(f'{out_dir}/{cell_cond}')
    DispCoeff_df.to_csv(f'{out_dir}/{cell_cond}/DispCoeff_{cell_cond}_{cm}.csv')
    
    # extract values for ALL comb and sort for k1 and k2
    DispCoeffALL_df = DispCoeff_df.loc[DispCoeff_df['Net_comb'] == 'ALL']
    DispCoeffALL_df = DispCoeffALL_df.sort_values(by = 'DispCoeff_avg', ascending=False)  
    
    DispCoeffALL_df_k1sorted = DispCoeffALL_df.sort_values(by = 'k1', ascending=False)  
    DispCoeffALL_df_k1sorted['k1'] = DispCoeffALL_df['k1'].astype(int)
    DispCoeffALL_df_k1sorted = DispCoeffALL_df_k1sorted.sort_values(by = 'k1', ascending=True)  
    DispCoeffALL_df_k2sorted = DispCoeffALL_df.sort_values(by = 'k2', ascending=False)  
    DispCoeffALL_df_k2sorted['k2'] = DispCoeffALL_df['k2'].astype(int)
    DispCoeffALL_df_k2sorted = DispCoeffALL_df_k2sorted.sort_values(by = 'k2', ascending=True)     
    
    DispCoeff_avglistk1 = DispCoeffALL_df_k1sorted.loc[:,'DispCoeff_avg'].tolist()
    k1s = DispCoeffALL_df_k1sorted.loc[:,'k1'].tolist()   
    k2s = DispCoeffALL_df_k2sorted.loc[:,'k2'].tolist()

    
    k1s = set(k1s)
    k2s = set(k2s)
    DispCoeff_df_k1k2 = pd.DataFrame(columns=list(k2s),index=list(k1s))

    for k1 in k1s:
        for k2 in k2s:
            try:
                DispCoeff = DispCoeffALL_df.loc[(DispCoeffALL_df['k1'] == str(k1)) & (DispCoeffALL_df['k2'] == str(k2)), 'DispCoeff_avg'].values[0]
                DispCoeff_df_k1k2.loc[k1,k2] = DispCoeff
            except IndexError:
                DispCoeff_df_k1k2.loc[k1,k2] = None
    
    DispCoeff_df_k1k2 = DispCoeff_df_k1k2.sort_index(axis=1)
    DispCoeff_df_k1k2 = DispCoeff_df_k1k2.sort_index(axis=0)
    DispCoeff_df_k1k2 = DispCoeff_df_k1k2.astype(float)

 
    
    fig, ax = plt.subplots(figsize=(11, 9)) 
    sb.lineplot(data=DispCoeff_df_k1k2[DispCoeff_df_k1k2.columns],legend='brief',linewidth = 4)
    # plt.axvline(x = bestk1[cell_cond], color = 'r')
    plt.plot(bestk1[cell_cond],DispCoeff_df_k1k2.loc[bestk1[cell_cond],bestk2[cell_cond]], color = 'r',marker='X',markersize=14)
    ax.tick_params(axis='x', labelsize= 22)
    ax.tick_params(axis='y', labelsize= 22)
    ax.set_title(f'{labels_key[cell_cond]}', fontsize = 28)
    ax.set_xlabel('k$_{1}$', fontsize = 24)
    ax.set_ylabel(r'$\mathrm{mean}_{\rho_{k_1},\rho_{k_2}}$', fontsize = 24)
    leg = ax.legend(title='k$_{2}$',title_fontsize=24, fontsize=20, loc='lower right',ncol=2)
    # change the line width for the legend
    for line in leg.get_lines():
        line.set_linewidth(4.0)
    # ax.legend(title='k$_{2}$',title_fontsize=18, fontsize=16, loc='lower right')
    ax.xaxis.grid(color='gray', alpha=0.3)
    ax.yaxis.grid(color='gray', alpha=0.3)
    
    plt.savefig(f'{out_dir}/{cell_cond}/DispCoeffML_{cell_cond}_{cm}.jpg',  dpi = 350, format='jpg')
    plt.close()
    









    
        