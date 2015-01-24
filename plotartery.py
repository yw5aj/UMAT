# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 16:09:19 2015

@author: Administrator
"""

# %% Load data
import numpy as np, pandas as pd, matplotlib.pyplot as plt
df = pd.read_excel('./csvs/radius_pressure.xlsx')
# Experiment data
press_exp = df['Pressure (kPa)'].dropna()
radius_exp = df['Radius (mm)'].dropna()
# Modeling data
press_num = df['Unnamed: 8'][2:10] * 1e-3
radius_num = df['Unnamed: 13'][2:10] * 1e3
press_ana = df['Unnamed: 8'][13:20] * 1e-3
radius_ana = df['Unnamed: 13'][13:20] * 1e3
# Stress distribution
radius_dist = np.linspace(1., 0., 7)
press_dist = np.linspace(0, 25000, 5) * 1e-3
press_num_raw = df.iloc[24:31, 6] * 1e-3
stress_num_raw = df.iloc[24:31, 8:15] * 1e-3
press_ana_raw = df.iloc[35:41, 6] * 1e-3
stress_ana_raw = df.iloc[35:41, 8:15] * 1e-3
stress_num = np.empty((press_dist.size, radius_dist.size))
stress_ana = np.empty((press_dist.size, radius_dist.size))
for i, key in enumerate(stress_num_raw):
    stress_num[:, i] = np.interp(press_dist, press_num_raw.astype(float),
        stress_num_raw[key].astype(float))
for i, key in enumerate(stress_ana_raw):
    stress_ana[:, i] = np.interp(press_dist, press_ana_raw.astype(float),
        stress_ana_raw[key].astype(float))
# %% Start plotting
fig, axs = plt.subplots(3, 1, figsize=(3.27, 8.75))
# Screenshot
im = plt.imread('./plots/screenshot_artery.png')
axs[0].imshow(im)
axs[0].axis('off')
axs[0].axvline(x=3, ls='--', lw=1.5, c='k', dashes=(8, 3, 2, 3))
axs[0].axhline(y=im.shape[0]-3, ls='--', lw=1.5, c='k', dashes=(8, 3, 2, 3))
axs[0].text(-0.01, 0.5, 'Symmetric axis', va='center', ha='right',
    rotation='vertical', transform=axs[0].transAxes, size=8)
axs[0].text(0.5, -0.01, 'Symmetric axis', va='top', ha='center',
    rotation='horizontal', transform=axs[0].transAxes, size=8)
# Pressure-radius
axs[1].plot(press_exp, radius_exp, '^', mfc='none', label='Experiment')
axs[1].plot(press_num, radius_num, '-ok', mfc='none', label='Numerical FE')
axs[1].plot(press_ana, radius_ana, '--xk', ms=5, label='Analytical FE')
axs[1].legend(loc=4)
axs[1].set_xlabel('Pressure (kPa)')
axs[1].set_ylabel('Radius (mm)')
# Distributed stress
for i in range(stress_num.shape[0]):
    axs[2].plot(radius_dist, stress_num[i], '-o', mfc='none', color=str(i*.15),
        label='Numerical FE')
    axs[2].plot(radius_dist, stress_ana[i], '--x', ms=5, color=str(i*.15),
        label='Analytical FE')
handles, labels = axs[2].get_legend_handles_labels()
axs[2].legend(handles[:2], labels[:2])
axs[2].set_xlabel('Normalized distance')
axs[2].set_ylabel('Max. principal stress (kPa)')
# Organize and save
fig.tight_layout()
for axes_id, axes in enumerate(axs.ravel()):
    axes.text(-.15, 1.05, chr(65+axes_id), transform=axes.transAxes,
        fontsize=12, fontweight='bold', va='top')
fig.savefig('./plots/arteryinfl.png')
plt.close(fig)