#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 17:11:14 2023

@author: 4vt
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

stats['Name'] = [n[0].upper() + n[1:] for n in stats['Name']]
stats = stats.drop(bad_aroe[0])
abundance['Name'] = [n[0].upper() + n[1:] for n in abundance['Name']]
abundance = abundance.drop(bad_aroe[1])

i_cols = ['A1','A2','A3','B1','B2','B3','C1','C2','C3']
code_map = [('CJ475', 'AD'), ('CJ680', 'BE'), ('CJ781', 'CF')]

sig_comps = {p:(s.split(';') if type(s) == str else '') for p,s in zip(stats['Name'], stats['Sample Differences'])}
# sig_comps['fpr-II'] = sig_comps['Fpr-II']
bar_edges = {'A_B':[-(1/6),0], 'A_C':[-(1/6),(1/6)], 'B_A':[-(1/6),0], 
             'B_C':[0,(1/6)], 'C_A':[-(1/6),(1/6)], 'C_B':[0,(1/6)]}
bh_scales = [1e7, 1.1e9]
rescale = [5e-8, 5e-10]
protlists = [['CspD','OmpR','AhpF','PP_3156','OxyR','AhpC','KatA','Ohr','SrkA','CorC','IbpA','PP_3234'],
              ['PP_0189','HemA','HemE','HemF','PP_3781','PP_1358','Bfr-I','FbpA','PP_3155','PP_3330']]
protlists = [sorted(p, key = lambda x: np.max(abundance[abundance['Name'] == x][i_cols])) for p in protlists]
llocs = ['upper left']*2

# all_prots = ['CspD','OmpR','AhpF','PP_3156','OxyR','AhpC','KatA','Ohr','SrkA','CorC','IbpA','PP_3234',
#              'PP_0189','HemA','HemE','HemF','HemN','PP_1358','Bfr-I','FbpA','PP_3155','PP_3330']
# all_prots = sorted(all_prots, key = lambda x: np.max(abundance[abundance['Name'] == x][i_cols]))
# protlists = [all_prots[:int(len(all_prots)/2)],all_prots[int(len(all_prots)/2):]]

cond_colors = {('CJ475', 'AD'):'brown', 
               ('CJ680', 'BE'):'gold', 
               ('CJ781', 'CF'):'paleturquoise'}
cond_objs = {}
size_scale = 0.7
fig, axes = plt.subplots(2,1, layout ='constrained', figsize = (8*size_scale,6*size_scale))
for i,ax in enumerate(axes):
    bh_scale = bh_scales[i]
    proteins = protlists[i]
    bar_heights = {'A_B':1*bh_scale, 'A_C':3*bh_scale, 'B_A':1*bh_scale, 
                 'B_C':2*bh_scale, 'C_A':3*bh_scale, 'C_B':2*bh_scale}
    prot_max = {p:np.max(abundance[abundance['Name'] == p][i_cols].to_numpy()) for p in proteins}

    prot_idx = {p:abundance[abundance['Name'] == p].index[-1] for p in proteins}
    
    for protein in proteins:
        for cond in code_map:
            cond_objs[cond] = ax.scatter([proteins.index(protein) - (1/6) + code_map.index(cond)/6]*3,
                                         abundance.loc[[prot_idx[protein]],
                                                       [c for c in abundance.columns if c.startswith(cond[1][0]) and len(c) == 2]],
                                         s = 7, c = cond_colors[cond])
        for comp in sig_comps[protein]:
            sigbar = ax.plot(np.asarray(bar_edges[comp]) + proteins.index(protein),
                             [prot_max[protein] + bar_heights[comp]*prot_max[protein]*rescale[i]]*2,
                             '-k', linewidth = 1)[0]
            ax.scatter(np.asarray(bar_edges[comp]) + proteins.index(protein),
                       [prot_max[protein] + bar_heights[comp]*prot_max[protein]*rescale[i]]*2,
                       s = 5, c = 'k', marker = '|')
    sigbar.set_label('Tukey’s HSD')
    _=[cond_objs[c].set_label(c[0]) for c in cond_objs.keys()]
    if i == 0:
        ax.legend(facecolor = 'whitesmoke', edgecolor = 'whitesmoke', loc = llocs[i],
                  fontsize = 7)
    ax.set_xticks(range(len(proteins)),
                  [p[0].upper() + p[1:] for p in proteins],
                  fontsize = 7)
    ax.set_xlim(-(2/6), len(proteins) -1 + (2/6))
    ax.ticklabel_format(axis = 'y', style = 'sci', useMathText = True)
    ax.set_box_aspect(1/3)
    ax.set_yscale('log')
    ax.yaxis.set_label_position("right")
    ax.set_ylabel('AB'[i], ha='left', y=1, rotation=0)
fig.supylabel('Summed Intensity')
fig.savefig('presentations_and_reports/Ikenna_paper/Figure6.png', 
            bbox_inches = 'tight', dpi = 1000)
    