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
from matplotlib import colormaps as cmaps
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

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
with open('Supp data file_PraI-PobA proteomics paper.xlsx','rb') as xlsx:
    raw_data = pd.read_excel(xlsx, 'S1', skiprows= 2)
raw_data.index = raw_data['Name']
#MISSING: Ech, AroY, EcdBD
queries = ['fcs', 'vdh', 'pobA', 'PraI', 'vanA', 'vanB', 'catA-I', 'catA-II']
icols = [f'{l}{n}' for l in 'ABC' for n in range(1,4)]
panel_B = raw_data.loc[queries, ['Name'] + icols]
vals = panel_B[icols].to_numpy()
means = np.nanmean(vals, axis = 1)
stds = np.nanstd(vals, axis = 1)
vals = (vals - means[:, np.newaxis])/stds[:, np.newaxis]
panel_B[icols] = vals
icols = ['CJ475 1', 'CJ475 2', 'CJ475 3',
         'CJ680 1', 'CJ680 2', 'CJ680 3',
         'CJ781 1', 'CJ781 2', 'CJ781 3']
panel_B.columns = ['Name'] + icols

absmax = np.nanmax(np.abs(panel_B[icols].to_numpy()))
norm = Normalize(vmin = -absmax, vmax = absmax)

# panel C data
with mpl.cbook.get_sample_data(f'{wd}Figure3_pathway.png') as file:
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

scale = 0.9
fig = plt.figure(layout="constrained", figsize = (6*scale, 5*scale))
gs = GridSpec(2, 4, figure=fig)

axA = fig.add_subplot(gs[0,:-1])
axB = fig.add_subplot(gs[:,-1])
axC = fig.add_subplot(gs[1,:-1])

#panel A
axA.scatter(top50['CJ475'], top50['CJ781'], s = 5, c = 'k', marker = '.', label = 'CJ781')
axA.scatter(top50['CJ475'], top50['CJ680'], s = 5, c = 'g', marker = '.', label = 'CJ680')

labels = [axA.text(cj475[p], label_pos[p], p, fontsize = 7) for p in annotate]
adjust_text(labels, 
            x = list(top50['CJ475'])*2, 
            y = list(top50['CJ781']) + list(top50['CJ680']), 
            ax=axA)

axA.plot([np.min(top50[data_cols]),np.max(top50['CJ475'])],
         [np.min(top50[data_cols]),np.max(top50['CJ475'])],
         '--r', linewidth = 0.5)
axA.legend(fontsize = 7, loc = 'upper center')
axA.set_title('A', loc = 'left')

pad = 3e9
axA.set_xlim(min(top50['CJ475']) - pad,max(top50['CJ475']) + pad)
axA.ticklabel_format(axis = 'both', style = 'sci', useMathText = True)
axA.tick_params(axis='both', labelsize=7)

# #panel B
axB.imshow(panel_B[icols], 
           aspect = 'auto', 
           interpolation = 'nearest', 
           cmap = cmaps['coolwarm'], 
           norm = norm)
axB.set_yticks(range(len(queries)), queries, fontsize = 8)
axB.set_xticks(range(len(icols)), icols, fontsize = 8, rotation = 90)
axB.yaxis.tick_right()

sm = ScalarMappable(norm = norm, cmap = cmaps['coolwarm'])
clb = plt.colorbar(sm, 
                   ax = axB, 
                   shrink = .7, 
                   aspect = 15,
                   orientation="horizontal",
                   location = 'top',
                   pad = -0.08)
clb.ax.set_title('Z-score', fontsize = 8)
clb.ax.tick_params(labelsize=8)
axB.set_title('B', loc = 'left')

#panel C
axC.imshow(panelC, aspect = 'equal')
no_borders(axC)
axC.set_title('C', loc = 'left')

fig.savefig('Figure4.svg', dpi = 900, bbox_inches = 'tight')
