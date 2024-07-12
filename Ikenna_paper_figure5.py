#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 15:01:23 2023

@author: 4vt
"""



import os 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec

os.chdir('/home/4vt/Documents/data/SLT01_PUFs/')
# rng = np.random.default_rng(1)

### protein data
with open('presentations_and_reports/Ikenna_paper/Supp data file_PraI-PobA proteomics paper.xlsx', 'rb') as xlsx:
    abundance = pd.read_excel(xlsx, 'S1', skiprows=2).fillna(0)
    stats = pd.read_excel(xlsx, 'S3', skiprows=1)

bad_aroe = [1554, 1857]

stats['Name'] = [n[0].upper() + n[1:] for n in stats['Name']]
stats = stats.drop(bad_aroe[0])
abundance['Name'] = [n[0].upper() + n[1:] for n in abundance['Name']]
abundance = abundance.drop(bad_aroe[1])

code_map = [('CJ475', 'AD'), ('CJ680', 'BE'), ('CJ781', 'CF')]

proteins = ['ZwfB', 'Fpr-II', 'NadK', 'AroE', 'PntAA', 'SthA', 'GntZ']

sig_comps = {p:(s.split(';') if type(s) == str else '') for p,s in zip(stats['Name'], stats['Sample Differences'])}
bar_edges = {'A_B':[-(1/6),0], 'A_C':[-(1/6),(1/6)], 'B_A':[-(1/6),0], 
             'B_C':[0,(1/6)], 'C_A':[-(1/6),(1/6)], 'C_B':[0,(1/6)]}

bh_scale = 1e8
bar_heights = {'A_B':1*bh_scale, 'A_C':3*bh_scale, 'B_A':1*bh_scale, 
             'B_C':2*bh_scale, 'C_A':3*bh_scale, 'C_B':2*bh_scale}
i_cols = ['A1','A2','A3','B1','B2','B3','C1','C2','C3']
prot_max = {p:np.max(abundance[abundance['Name'] == p][i_cols].to_numpy()) for p in proteins}

prot_idx = {p:abundance[abundance['Name'] == p].index[-1] for p in proteins}
cond_colors = {('CJ475', 'AD'):'brown', 
               ('CJ680', 'BE'):'gold', 
               ('CJ781', 'CF'):'paleturquoise'}

### pathway image
with mpl.cbook.get_sample_data('/home/4vt/Documents/data/SLT01_PUFs/presentations_and_reports/Ikenna_paper/old_files/figure3A.png') as file:
    arr_image = plt.imread(file, format='png')

### plot figure
cond_objs = {}
gs = gridspec.GridSpec(1,2, width_ratios= [arr_image.shape[1]/arr_image.shape[0],1])
fig = plt.figure()
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

#panel B
for protein in proteins:
    for cond in code_map:
        cond_objs[cond] = ax2.scatter([proteins.index(protein) - (1/6) + code_map.index(cond)/6]*3,
                                     abundance.loc[[prot_idx[protein]],
                                                   [c for c in abundance.columns if c.startswith(cond[1][0]) and len(c) == 2]],
                                     s = 7, c = cond_colors[cond])
    for comp in sig_comps[protein]:
        sigbar = ax2.plot(np.asarray(bar_edges[comp]) + proteins.index(protein),
                         [prot_max[protein] + bar_heights[comp]]*2, '-k', linewidth = 1)[0]
        ax2.scatter(np.asarray(bar_edges[comp]) + proteins.index(protein),
                   [prot_max[protein] + bar_heights[comp]]*2, s = 5, c = 'k', marker = '|')
sigbar.set_label('Tukey’s HSD')
_=[cond_objs[c].set_label(c[0]) for c in cond_objs.keys()]
ax2.legend(facecolor = 'whitesmoke', edgecolor = 'whitesmoke', 
           loc = 'upper left', borderpad=0.2)
ax2.set_xticks(range(len(proteins)),
              [p[0].upper() + p[1:] for p in proteins],
              fontsize = 9)
ax2.set_xlim(-(2/6), len(proteins) -1 + (2/6))
ax2.ticklabel_format(axis = 'y', style = 'sci', useMathText = True)
ax2.set_box_aspect(1)
ax2.set_ylabel('Summed Intensity')
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
ax2.set_title('B', loc = 'left')
y0,y1 = ax2.get_ylim()
ax2.set_ylim((0,y1))

#panel A
ax1.imshow(arr_image)
for spine in ['top','right','left','bottom']:
    ax1.spines[spine].set_visible(False)
ax1.set_title('A', loc = 'left')
ax1.set_yticks([],[])
ax1.tick_params(axis='y',
               which='both',
               left=False,
               right=False,
               labelbottom=False)
ax1.tick_params(axis='x',
               which='both',
               bottom=False,
               top=False,
               labelbottom=False)

gs.tight_layout(fig)
fig.savefig('presentations_and_reports/Ikenna_paper/Figure5.png', 
            bbox_inches = 'tight', dpi = 1000)



