# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:07:13 2022

@author: goneill
"""


import numpy as np
import os.path as op
import matplotlib.pyplot as plt

import sys
sys.path.append('..//') # Use this is running tests without a pip install!

import opym
import mne


data_root = op.abspath('D:\\20220623_095714_cMEG_Data')
data_bin = op.join(data_root,'20220623_095714.cMEG')

rawC = opym.io.read_raw_cerca(data_bin,slotnames=True)
a = rawC[76,:]
plt.plot(a[1],a[0].transpose())



#%% Empty room SSP

empty_room_projs = mne.compute_proj_raw(rawC, n_mag=3)
rawC.add_proj(empty_room_projs,remove_existing=True)


# for proj in (False, True):
#     with mne.viz.use_browser_backend('matplotlib'):
#         fig = rawC.plot(butterfly=True, proj=proj)


#%%


rawHFC = opym.preproc.denoise_hfc(rawC)
a = rawHFC[76,:]
plt.plot(a[1],a[0].transpose())


rawHFC2 = opym.preproc.denoise_hfc(rawC,method='ssp')
a = rawHFC2[76,:]
plt.plot(a[1],a[0].transpose())


r2 = rawHFC2.get_data(76)
plt.plot(a[1],r2)