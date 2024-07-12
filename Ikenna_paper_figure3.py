#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 18:25:15 2023

@author: 4vt
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text

os.chdir('/home/4vt/Documents/data/SLT01_PUFs/')


with open('presentations_and_reports/Ikenna_paper/Supp data file_PraI-PobA proteomics paper.xlsx', 'rb') as xlsx:
    abundance = pd.read_excel(xlsx, 'S1', skiprows=2).fillna(0)
    stats = pd.read_excel(xlsx, 'S3', skiprows=1)

stats['Name'] = [n[0].upper() + n[1:] for n in stats['Name']]
stats['Sample Differences'] = [s if type(s) == str else '' for s in stats['Sample Differences']]
abundance['Name'] = [n[0].upper() + n[1:] for n in abundance['Name']]

anova_p = {n:p for n,p in zip(stats['Name'], stats['minus log10 pval'])}
tukey_hsd = {n:[''.join(sorted(i.split('_'))) for i in t.split(';')] for n,t in zip(stats['Name'],stats['Sample Differences'])}
code_map = {'A':'CJ475', 'B':'CJ680', 'C':'CJ781'}
color_map = {'A':'brown', 'B':'gold', 'C':'paleturquoise'}

comps = ['AC','BC','AB']
cond_cols = {i:[c for c in abundance.columns if c.startswith(i) and len(c) == 2] for i in 'ABC'}
for comp in comps:
    l2fc = np.log2(np.nanmean(abundance[cond_cols[comp[0]]], axis = 1)/np.nanmean(abundance[cond_cols[comp[1]]], axis = 1))
    l2fc = [n if np.isfinite(n) else np.sign(n) * 9 for n in l2fc]
    abundance[f'l2fc_{comp}'] = l2fc

up_color = '#BB5566'
down_color = '#004488'
def color(proteins, l2fc, comp):
    colors = []
    for p, l in zip(proteins, l2fc):
        if comp in tukey_hsd[p]:
            colors.append(down_color if l < 0 else up_color)
        else:
            colors.append('k')
    return colors

queries = {'AC':['PobA','PP_3784','PP_3781','PP_3785','PP_3775','ProC','PP_3536','PP_3787',
                 'PP_3777','AspC','PP_3788','CycA','PP_3776','PP_2099','PP_1146','PraI','PP_0218',
                 'CspD','HemE','OmpR'],
           'BC':['PobA','PP_3784','PP_3775','PP_3536','PP_0218','PcaY','FecA','PP_3550','PP_4296',
                 'Ffh','YijP','PP_1937','PP_1289'],
           'AB':['PP_3784','SyrB','PP_3777','PP_3785','PP_3775','ProC','PP_3781','PP_3787',
                 'PP_1104','AspC','PP_3788','FpvA','CycA','PP_3776','PP_0614','PP_1146','KdsA2',
                 'PraI','FecA','PP_3316','PP_0317','CspD']}

size = 8
fig, axes = plt.subplots(nrows = 1, ncols = 3, sharey = True, figsize = (6.5, 3),
                         layout = 'constrained') #(original figure (2.59, 6.5))
for comp,ax in zip(comps, axes):
    colors = color(abundance['Name'], abundance[f'l2fc_{comp}'], comp)
    ndown = len([c for c in colors if c == down_color])
    nup = len([c for c in colors if c == up_color])
    ax.scatter(abundance[f'l2fc_{comp}'],
               [anova_p[n] for n in abundance['Name']],
               s = 2, c = colors, marker = '.')
    query_idx = [abundance[abundance['Name'] == q].index[0] for q in queries[comp]]
    abu_dict = {n:l for n,l in zip(abundance['Name'], abundance[f'l2fc_{comp}'])}
    labels = [ax.text(abu_dict[q], anova_p[q], q,
                      fontsize = 5) for q in queries[comp] if np.isfinite(abu_dict[q]) and np.isfinite(anova_p[q])]
    adjust_text(labels, ax=ax, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    ax.set_xticks([-9,-5,0,5,9],
                  [code_map[comp[1]],'-5','0','5',code_map[comp[0]]],
                  rotation = 60, fontsize = size)
    xshift = 6
    ax.text(-xshift, 1, str(ndown), color = down_color, ha = 'center', va = 'center')
    ax.text(xshift, 1, str(nup), color = up_color, ha = 'center', va = 'center')
    if comps.index(comp):
        ax.tick_params(axis='y',
                       which='both',
                       left=False,
                       right=False,
                       labelbottom=False)
    else:
        ax.set_ylabel('Negative Log 2 ANOVA p-value', fontsize = size)
        ax.tick_params(axis='y', labelsize = size)
    ax.set_title('ABC'[comps.index(comp)], loc = 'left', fontsize = size)
    y0,y1 = ax.get_ylim()
    ax.set_ylim(0,y1)
fig.supxlabel('Log 2 Fold Change', fontsize = size)
fig.savefig('presentations_and_reports/Ikenna_paper/Figure3.png',
            bbox_inches = 'tight', dpi = 1000)

