# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:36:08 2023

@author: Administrator
"""

import os 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

os.chdir('/home/4vt/Documents/data/SLT01_PUFs/')
# rng = np.random.default_rng(1)

with open('presentations_and_reports/Ikenna_paper/Supp data file_PraI-PobA proteomics paper.xlsx', 'rb') as xlsx:
    abundance = pd.read_excel(xlsx, 'S1', skiprows=2).fillna(0)
    stats = pd.read_excel(xlsx, 'S3', skiprows=1)

bad_aroe = [1554, 1857]

stats['Name'] = [n[0].upper() + n[1:-1] + n[-1].upper() for n in stats['Name']]
stats = stats.drop(bad_aroe[0])
abundance['Name'] = [n[0].upper() + n[1:-1] + n[-1].upper() for n in abundance['Name']]
abundance = abundance.drop(bad_aroe[1])

code_map = [('CJ475', 'AD'), ('CJ680', 'BE'), ('CJ781', 'CF')]

def plot_prots(proteins, bh_scale, name, lloc = 'upper center'):
    sig_comps = {p:(s.split(';') if type(s) == str else '') for p,s in zip(stats['Name'], stats['Sample Differences'])}
    # sig_comps['fpr-II'] = sig_comps['Fpr-II']
    bar_edges = {'A_B':[-(1/6),0], 'A_C':[-(1/6),(1/6)], 'B_A':[-(1/6),0], 
                 'B_C':[0,(1/6)], 'C_A':[-(1/6),(1/6)], 'C_B':[0,(1/6)]}
    
    bar_heights = {'A_B':1*bh_scale, 'A_C':3*bh_scale, 'B_A':1*bh_scale, 
                 'B_C':2*bh_scale, 'C_A':3*bh_scale, 'C_B':2*bh_scale}
    i_cols = ['A1','A2','A3','B1','B2','B3','C1','C2','C3']
    prot_max = {p:np.max(abundance[abundance['Name'] == p][i_cols].to_numpy()) for p in proteins}

    prot_idx = {p:abundance[abundance['Name'] == p].index[-1] for p in proteins}
    cond_colors = {('CJ475', 'AD'):'brown', 
                   ('CJ680', 'BE'):'gold', 
                   ('CJ781', 'CF'):'paleturquoise'}
    cond_objs = {}
    fig, ax = plt.subplots(figsize = (4,4))
    for protein in proteins:
        for cond in code_map:
            cond_objs[cond] = ax.scatter([proteins.index(protein) - (1/6) + code_map.index(cond)/6]*3,
                                         abundance.loc[[prot_idx[protein]],
                                                       [c for c in abundance.columns if c.startswith(cond[1][0]) and len(c) == 2]],
                                         s = 7, c = cond_colors[cond])
        for comp in sig_comps[protein]:
            sigbar = ax.plot(np.asarray(bar_edges[comp]) + proteins.index(protein),
                             [prot_max[protein] + bar_heights[comp]]*2, '-k', linewidth = 1)[0]
            ax.scatter(np.asarray(bar_edges[comp]) + proteins.index(protein),
                       [prot_max[protein] + bar_heights[comp]]*2, s = 5, c = 'k', marker = '|')
    sigbar.set_label('Tukey’s HSD')
    _=[cond_objs[c].set_label(c[0]) for c in cond_objs.keys()]
    ax.legend(facecolor = 'whitesmoke', edgecolor = 'whitesmoke', loc = lloc)
    ax.set_xticks(range(len(proteins)),
                  [p[0].upper() + p[1:] for p in proteins],
                  rotation = 45, ha='right')
    ax.set_xlim(-(2/6), len(proteins) -1 + (2/6))
    ax.ticklabel_format(axis = 'y', style = 'sci', useMathText = True)
    ax.set_box_aspect(1)
    ax.set_ylabel('Summed Intensity')
    y0, y1 = ax.get_ylim()
    ax.set_ylim(0,y1)
    fig.savefig(f'presentations_and_reports/Ikenna_paper/{name}.png', 
                bbox_inches = 'tight', dpi = 1000)

#FIGURE 3B
proteins = ['ZwfB', 'Fpr-II', 'NadK', 'AroE', 'PntAA', 'SthA', 'GntZ']
bh_scale = 1e8
name = 'Figure3B'
plot_prots(proteins, bh_scale, name)

#FIGURE S1
proteins = ['PP_3775', 'PP_3776','PP_3777','ProC']
bh_scale = 5e8
name = 'FigureS1'
plot_prots(proteins, bh_scale, name)

#FIGURE S2
proteins = ['PP_3784', 'PP_3785','AspC','PP_3787','PP_3788']
bh_scale = 4e8
name = 'FigureS2'
plot_prots(proteins, bh_scale, name)

#FIGURE S3
proteins = ['PobA', 'PP_3536']
bh_scale = 2e9
name = 'FigureS3'
plot_prots(proteins, bh_scale, name, lloc = 'upper right')

#figure S5
proteins = ['AacS','TesA','EstP','MccB','FabZ','FabV','PP_4379']
bh_scale = 8e7
name = 'FigureS5'
plot_prots(proteins, bh_scale, name, lloc = 'upper left')


# for p in proteins:
#     np.max(abundance[abundance['Name'] == p][i_cols].to_numpy())
