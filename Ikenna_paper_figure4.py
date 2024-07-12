#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 16:41:31 2023

@author: 4vt
"""


import os 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from adjustText import adjust_text

wd = '/home/4vt/Documents/data/SLT01_PUFs/presentations_and_reports/Ikenna_paper/'
os.chdir(wd)

# panel A data
with open('Supp data file_PraI-PobA proteomics paper.xlsx','rb') as xlsx:
    top50 = pd.read_excel(xlsx, 'S2', skiprows= 1)
data_cols = ['CJ475', 'CJ781', 'CJ680']
top50[data_cols] = top50[data_cols].fillna(0)
annotate = ['PraI', 'PP_4981', 'PobA', 'VanB', 'RplL', 'PP_3358', 'GroL', 'OprF', 'TufB', 'VanA']
label_pos = {p:np.nansum([i6,i7])/2 for p,i7,i6 in zip(top50['Proteins'], top50['CJ781'], top50['CJ680'])}
cj475 = {p:i if np.isfinite(i) else min(top50['CJ475']) for p,i in zip(top50['Proteins'], top50['CJ475'])}

# panel B data
with mpl.cbook.get_sample_data(f'{wd}old_files/Figure4B.png') as file:
    panelB = plt.imread(file, format='png')

# panel C data
with mpl.cbook.get_sample_data(f'{wd}old_files/Figure4C.png') as file:
    panelC = plt.imread(file, format='png')

### plot figure

def no_borders(ax):
    for spine in ['top','right','left','bottom']:
        ax.spines[spine].set_visible(False)
    ax.set_yticks([],[])
    ax.tick_params(axis='y',
                   which='both',
                   left=False,
                   right=False,
                   labelbottom=False)
    ax.tick_params(axis='x',
                   which='both',
                   bottom=False,
                   top=False,
                   labelbottom=False)

fig = plt.figure(layout="constrained", figsize = (6.5, 4))
gs = GridSpec(2, 2, figure=fig)

axA = fig.add_subplot(gs[0,:-1])
axB = fig.add_subplot(gs[0,-1])
axC = fig.add_subplot(gs[1,:])

#panel A
axA.scatter(top50['CJ475'], top50['CJ781'], s = 5, c = 'k', marker = '.', label = 'CJ781')
axA.scatter(top50['CJ475'], top50['CJ680'], s = 5, c = 'g', marker = '.', label = 'CJ680')

labels = [axA.text(cj475[p], label_pos[p], p, fontsize = 5) for p in annotate]
adjust_text(labels, ax=axA)

axA.plot([np.min(top50[data_cols]),np.max(top50['CJ475'])],
         [np.min(top50[data_cols]),np.max(top50['CJ475'])],
         '--r', linewidth = 0.5)
axA.legend(fontsize = 7)
axA.set_title('A', loc = 'left')

pad = 3e9
axA.set_xlim(min(top50['CJ475']) - pad,max(top50['CJ475']) + pad)
axA.ticklabel_format(axis = 'both', style = 'sci', useMathText = True)
axA.tick_params(axis='both', labelsize=7)

#panel B
axB.imshow(panelB)
no_borders(axB)
axB.set_title('B', loc = 'left')

#panel C
axC.imshow(panelC)
no_borders(axC)
axC.set_title('C', loc = 'left')


