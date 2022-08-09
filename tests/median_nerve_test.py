# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 13:32:23 2022

@author: goneill
"""

import os.path as op
import sys
sys.path.append('../') # Use this is running tests without a pip install!

import opym
import mne
import numpy as np
import matplotlib.pyplot as plt
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
raw2.filter(1, 80)


#%% Epoch, average etc

events = mne.find_events(raw2, stim_channel='TRIG1')
epochs = mne.Epochs(raw2, events, tmin=-0.1, tmax=0.3, baseline=(-0.1, 0.), preload=True)
evoked = epochs.average()
evoked.plot()

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

fwd = mne.make_forward_solution(epochs.info, trans=coreg.trans, src=src, bem=bem,
                                meg=True, eeg=False, mindist=5.0, verbose=True)

#%% MNE?

noise_cov = mne.compute_covariance(
    epochs, tmax=0., method='auto', rank=None, verbose=True)

fig_cov, fig_spectra = mne.viz.plot_cov(noise_cov, epochs.info)

evoked.plot_white(noise_cov)

inverse_operator = mne.minimum_norm.make_inverse_operator(
    evoked.info, fwd, noise_cov, loose=0.2, depth=0.8)


method = "dSPM"
snr = 3.
lambda2 = 1. / snr ** 2
stc, residuals = mne.minimum_norm.apply_inverse(evoked, inverse_operator, lambda2,
                              method=method, pick_ori=None,
                              return_residual=True, verbose=True)

#%% Visualise result (hopefully!)

vertno_max, time_max = stc.get_peak(hemi='rh')

surfer_kwargs = dict(
    hemi='rh', subjects_dir=subjects_dir,
    clim=dict(kind='value', lims=[8, 12, 15]), views='lateral',
    initial_time=time_max, time_unit='s', size=(800, 800), smoothing_steps=10)
brain = stc.plot(**surfer_kwargs)
brain.add_foci(vertno_max, coords_as_verts=True, hemi='rh', color='blue',
                scale_factor=0.6, alpha=0.5)

#%%


# dip_opm, _ = mne.fit_dipole(evoked.copy().crop(0.040, 0.080),
#                             noise_cov, bem, coreg.trans, verbose=True)
# idx = np.argmax(dip_opm.gof)
# print('Best dipole at t=%0.1f ms with %0.1f%% GOF'
#       % (1000 * dip_opm.times[idx], dip_opm.gof[idx]))

# # Plot N20m dipole as an example
# dip_opm.plot_locations(coreg.trans, subject, subjects_dir,
#                        mode='orthoview', idx=idx)