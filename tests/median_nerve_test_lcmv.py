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
import numpy as np
import matplotlib.pyplot as plt
from mne.beamformer import make_lcmv, apply_lcmv
#%% the main import and coreg

# import the data
data_root = op.abspath('D:\old_mn_data')
data_bin = op.join(data_root,'sub-001_ses-001_task-mediannerve_run-001_meg.bin')
raw = opym.io.read_raw_ucl(data_bin)

# Lets have a look at trying to align data in the source space
subjects_dir = op.abspath('D:\\recons')
subject = 'masters-22-001'

coreg = opym.coreg.coreg_headcast(raw.info,subject=subject, subjects_dir=subjects_dir)

# View to check the registration looks correct
opym.coreg.check_datareg(raw.info, coreg=coreg, subject=subject, subjects_dir=subjects_dir)

#%% HFC denoising - just the homogenous bit and filter

raw2 = opym.preproc.denoise_hfc(raw)
# raw2 = raw.copy()
raw2.filter(1, 80)


#%% Epoch, average etc

events = mne.find_events(raw2, stim_channel='TRIG1')
epochs = mne.Epochs(raw2, events, tmin=-0.1, tmax=0.3, baseline=(-0.1, 0.), preload=True)
evoked = epochs.average()
evoked.plot()


#%% Covariance matricies

data_cov = mne.compute_covariance(epochs, tmin=0.01, tmax=0.3,
                                  method='empirical')

noise_cov = mne.compute_covariance(epochs,tmin=-0.1, tmax=0.0,
                                  method='empirical')



data_cov.plot(epochs.info)

#%% Source space and BEM construction

src = mne.setup_source_space(subject, spacing='oct6', add_dist='patch',
                             subjects_dir=subjects_dir)
# 
conductivity = (0.3,)  # for single layer
# conductivity = (0.3, 0.006, 0.3)  # for three layers
model = mne.make_bem_model(subject=subject, ico=4,
                           conductivity=conductivity,
                           subjects_dir=subjects_dir)
bem = mne.make_bem_solution(model)

#%% solve forwards

# coil_def_fname = op.abspath('C:\\Users\\goneill\\mne_data\\MNE-OPM-data\\MEG\\OPM\\coil_def.dat')


fwd = mne.make_forward_solution(epochs.info, trans=coreg.trans, src=src, bem=bem,
                    meg=True, eeg=False, mindist=5.0, verbose=False)
    
#%% Inverse time

filters = make_lcmv(evoked.info, fwd, data_cov, reg=0.05,
                    pick_ori='max-power',
                    noise_cov=noise_cov,weight_norm='unit-noise-gain', rank=None)


stc = apply_lcmv(evoked, filters)

#%% Visualise result (hopefully!)

vertno_max, time_max = stc.get_peak(hemi='rh')

surfer_kwargs = dict(
    hemi='rh', subjects_dir=subjects_dir, surface='inflated',
    clim=dict(kind='value', lims=[0.7*stc.data.max(), 0.85*stc.data.max(), stc.data.max()]), views='lateral',
    initial_time=time_max, time_unit='s', size=(800, 800), smoothing_steps=10)
brain = stc.plot(**surfer_kwargs)

#%%


# dip_opm, _ = mne.fit_dipole(evoked.copy().crop(0.040, 0.080),
#                             noise_cov, bem, coreg.trans, verbose=True)
# idx = np.argmax(dip_opm.gof)
# print('Best dipole at t=%0.1f ms with %0.1f%% GOF'
#       % (1000 * dip_opm.times[idx], dip_opm.gof[idx]))

# # Plot N20m dipole as an example
# dip_opm.plot_locations(coreg.trans, subject, subjects_dir,
#                        mode='orthoview', idx=idx)