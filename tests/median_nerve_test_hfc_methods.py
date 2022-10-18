# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 13:32:23 2022

@author: goneill
"""

import os.path as op
import sys
sys.path.append("D:\\Documents\\GitHub\\OPyM") # Use this is running tests without a pip install!

import opym
import mne

#%% the main import and coreg

# import the data
data_root = op.abspath('D:\old_mn_data')
data_bin = op.join(data_root,'sub-001_ses-001_task-mediannerve_run-001_meg.bin')
raw = opym.io.read_raw_ucl(data_bin)

#%% no HFC

raw4 = raw.copy()
raw4.load_data()
raw4.filter(1, 80)

events = mne.find_events(raw4, stim_channel='TRIG1')
epochs = mne.Epochs(raw4, events, tmin=-0.1, tmax=0.3, baseline=(-0.1, 0.), preload=True)
evoked = epochs.average()
evoked.plot()


#%% HFC denoising - just the homogenous bit and filter - using old method!

raw2 = opym.preproc.denoise_hfc(raw)
# raw2 = raw.copy()
raw2.filter(1, 80)


events = mne.find_events(raw2, stim_channel='TRIG1')
epochs = mne.Epochs(raw2, events, tmin=-0.1, tmax=0.3, baseline=(-0.1, 0.), preload=True)
evoked = epochs.average()
evoked.plot()


#%% HFC denoising - just the homogenous bit and filter - using SSP as required!

raw3 = opym.preproc.denoise_hfc(raw,update_proj=True)
raw3.load_data()
raw3.filter(1, 80)

events = mne.find_events(raw3, stim_channel='TRIG1')
epochs = mne.Epochs(raw3, events, tmin=-0.1, tmax=0.3, baseline=(-0.1, 0.), preload=True, proj=False)
evoked = epochs.average()
evoked.plot()




