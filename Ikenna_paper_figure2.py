#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 13:59:11 2023

@author: 4vt
"""
import os
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
from matplotlib import colormaps as cmaps
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


os.chdir('/home/4vt/Documents/data/SLT01_PUFs/')

conditions = ['A','B','C','D','E','F']
cond_map = {'A':'CJ475 p-CA','B':'CJ680 p-CA','C':'CJ781 p-CA',
            'D':'CJ475 glu', 'E':'CJ680 glu', 'F':'CJ781 glu'}

with open('presentations_and_reports/Ikenna_paper/filtered raw data_A-F.xlsx','rb') as xlsx:
    excel = {c:pd.read_excel(xlsx, c) for c in conditions}

# frames = []
# for cond in conditions:
#     c_frame = pd.DataFrame({f'{cond}{i}':excel[cond][[c for c in excel[cond].columns if c.startswith('Abundance')][i]] for i in range(3)})
#     c_frame.index = excel[cond]['Accession']
#     frames.append(c_frame)

with open('presentations_and_reports/Ikenna_paper/filtered raw data_A-F.xlsx','rb') as xlsx:
    a_f = pd.read_excel(xlsx, 'A-F')
intensities = a_f[list(a_f.columns)[3:]]

cols = intensities.columns
rows = intensities.index
intensities = intensities.to_numpy()
medians = np.nanmedian(intensities, axis = 0)
intensities = intensities - [medians]*intensities.shape[0]
intensities = intensities + np.mean(medians)
zscores = (intensities - np.asarray([np.nanmean(intensities, axis = 1)]*intensities.shape[1]).T)/np.asarray([np.nanstd(intensities, axis = 1)]*intensities.shape[1]).T

dendro_metric = 'cosine'
prot_cluster = AgglomerativeClustering(n_clusters=None, 
                                       distance_threshold= np.nansum(abs(zscores)),
                                       compute_full_tree= True,
                                       linkage = 'average',
                                       metric = dendro_metric).fit(np.nan_to_num(zscores,  copy = True, nan = -4))
dendro_out = dendrogram(np.column_stack([prot_cluster.children_, prot_cluster.distances_, list(range(len(prot_cluster.children_)))]))
prot_order = dendro_out['leaves']


samp_cluster = AgglomerativeClustering(n_clusters=None, 
                                       distance_threshold= np.nansum(abs(zscores)),
                                       compute_full_tree= True,
                                       linkage = 'average',
                                       metric = dendro_metric).fit(np.nan_to_num(zscores,  copy = True, nan = -4).T)
dendro_out = dendrogram(np.column_stack([samp_cluster.children_, samp_cluster.distances_, list(range(len(samp_cluster.children_)))]))
samp_order = dendro_out['leaves']

zscores = pd.DataFrame(zscores, columns = cols, index = rows)
zscores = zscores[zscores.columns[samp_order]]
zscores = zscores.loc[zscores.index[prot_order]]


fig = plt.figure(figsize = (6,10), layout = 'constrained')
gs = fig.add_gridspec(5, 1)
ax1 = fig.add_subplot(gs[0])
dendrogram(np.column_stack([samp_cluster.children_, samp_cluster.distances_, list(range(len(samp_cluster.children_)))]),
           ax = ax1, labels = [f'{cond_map[cols[i][0]]} {cols[i]}' for i in samp_order], leaf_rotation=90, 
           color_threshold=0, link_color_func = lambda k: 'k')
ax1.set_yticks([])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(False)

absmax = np.max(np.abs(zscores.to_numpy()))
norm = Normalize(vmin = -absmax, vmax = absmax)

ax2 = fig.add_subplot(gs[1:])
ax2.imshow(zscores, aspect = 'auto', interpolation = 'nearest', cmap = cmaps['coolwarm'], norm = norm)
ax2.set_yticks([])
ax2.set_xticks([])

sm = ScalarMappable(norm = norm, cmap = cmaps['coolwarm'])
clb = plt.colorbar(sm, ax = [ax1,ax2], shrink = .37, aspect = 20)
clb.ax.set_title('Z-score')

fig.savefig('presentations_and_reports/Ikenna_paper/Figure2.svg', bbox_inches = 'tight')
fig.savefig('presentations_and_reports/Ikenna_paper/Figure2.png', bbox_inches = 'tight', dpi = 1000)
