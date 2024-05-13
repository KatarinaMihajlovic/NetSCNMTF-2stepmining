# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 20:00:30 2021

@author: kmihajlo
"""

import pickle, os
import pandas as pd 


labels_key = {'Control_IPSCs':'C0', 'Control_D06':'C6', 'Control_D15':'C15', 'Control_D21':'C21', 
              'PINK1_IPSCs':'PD0', 'PINK1_D06':'PD6', 'PINK1_D15':'PD15', 'PINK1_D21':'PD21'}



def Sort(sub_li):
  
    sub_li.sort(key = lambda x: x[1], reverse=True)
    return sub_li



### Main Code

with open('input/PD_genes_DGN.pkl', 'rb') as handle:
    PD_genes = pickle.load(handle)  

All_genes = pd.read_csv('input/All_genes.csv', header = None)      
All_genes = All_genes[0].tolist()  
All_genes = [str(x) for x in All_genes]


for root, dirs, files in os.walk('input/Clusters'):
    for file in files:
        print(file)
        Init_preds = []
        cell_cond = root.split('\\')[1]
        print(cell_cond)
        cm = file.split('_PD')[0]

        with open(f'{root}/{file}', 'rb') as handle:
            ClustsPDEnrich = pickle.load(handle)
        
        Geneslist = pd.read_csv(f'input/Genelists/Geneslist_{cell_cond}.csv')
        Geneslist = Geneslist['genes'].tolist()

        i = 0
        for clst in range(len(ClustsPDEnrich)):
            PD_gene_Predictions = []
            EnrClust = ClustsPDEnrich[clst]
            Litvalid_genes = []
            for gene in EnrClust:
                if gene not in PD_genes:
                    PD_gene_Predictions.append(gene)
                    i+=1
            Init_preds.append(PD_gene_Predictions)
        Init_preds = [j for i in Init_preds for j in i]
        print(len(Init_preds))
            
        with open(f'output/{labels_key[cell_cond]}_{cm}_InitPreds.txt', 'w') as f:
            for gene in Init_preds:
                f.write(f'{gene}\n')
        with open(f'output/{labels_key[cell_cond]}_{cm}_InitPreds.pkl', 'wb') as handle:
            pickle.dump(Init_preds, handle)


