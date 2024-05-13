# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 23:28:36 2024

@author: Katarina
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


nets = ['ALL','PartMImissing', 'PPI+GI+COEX']
nets_l = {'ALL':'E+PPI+MI+COEX+GI','PartMImissing':'E+PPI+MI$_{CF}$+COEX+GI', 'PPI+GI+COEX':'E+PPI+COEX+GI'}

AnnotPerc_tot = []
case_tot = []
net_type_tot = []
for net in nets:    
    print(net)
    KPperc = np.load(f'input/{net}_KPperc_TestCorrectClst.npy').flatten()
    RPperc = np.load(f'input/{net}_RPperc_TestCorrectClst.npy').flatten()
    BPperc = np.load(f'input/{net}_BPperc_TestCorrectClst.npy').flatten()
    
    print(np.mean(KPperc))
    print(np.mean(RPperc))
    print(np.mean(BPperc))

    
    AnnotPerc = np.concatenate((KPperc, RPperc, BPperc))
    case = ['KP']*len(KPperc) + ['RP']*len(RPperc) + ['GO']*len(BPperc)
    net_type = [nets_l[net]]*len(AnnotPerc)
    
    AnnotPerc_tot.append(AnnotPerc)
    case_tot.append(case)
    net_type_tot.append(net_type)
    
AnnotPerc_tot = [item for items in AnnotPerc_tot for item in items]
case_tot = [item for items in case_tot for item in items]
net_type_tot = [item for items in net_type_tot for item in items]

AnnotPerc_df = pd.DataFrame({'AnnotPerc': AnnotPerc_tot, 'case': case_tot,'Net': net_type_tot})

fig, ax = plt.subplots(figsize=(12, 8))
# plt.grid(axis='y', alpha = 0.5)
ax = sns.boxplot(data=AnnotPerc_df, x='case', y='AnnotPerc', hue = 'Net', palette=sns.color_palette('pastel'))
# annotator = Annotator(ax, [['FoldTest','FoldBG']], data=avg_enrichments_df, x='case', y='value', order=['FoldTest','FoldBG'])
# annotator.configure(test='Mann-Whitney-gt', text_format='full',fontsize=22, loc='outside')
# annotator.apply_and_annotate()
ax.set_ylabel('% of genes with at least one of their\nannotations enriched in their clusters', fontsize = 24)#, pad = 24)
# correcr cluster - cluster enriched in at least one annotation that annotates the Test gene
ax.set_xlabel('')
ax.legend(fontsize=16, title='NetSC-NMTF',title_fontsize=20)
# ax.set_xticklabels(['$PD_{{DGN}_{Test}}$','Background'])
plt.xticks(fontsize=24)
plt.yticks(fontsize=20)
plt.savefig('output/5FoldTestgenes_CorrectClst_PathAnnot.jpg',  dpi = 350, format='jpg')  

plt.show()
plt.close()     
