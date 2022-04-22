# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 15:59:52 2022

@author: goneill


"""
import mne
from opym import read_raw_ucl
import os.path as op
from mne.datasets import fetch_fsaverage

# Download fsaverage files
fs_dir = fetch_fsaverage(verbose=True)
subjects_dir = op.dirname(fs_dir)

# The files live in:
subject = 'fsaverage'
trans = 'fsaverage'  # MNE has a built-in fsaverage transformation
src = op.join(fs_dir, 'bem', 'fsaverage-ico-5-src.fif')
bem = op.join(fs_dir, 'bem', 'fsaverage-5120-5120-5120-bem-sol.fif')

#read in OPM dataset
data_root = op.abspath('D:\masters_example_data\data')
data_bin = op.join(data_root,'sub-001_ses-001_task-motor4way_run-001_meg.bin')

raw = read_raw_ucl(data_bin)

# Check that the locations of EEG electrodes is correct with respect to MRI
mne.viz.plot_alignment(
    raw.info, src=src, eeg=['original', 'projected'], trans=trans,
    show_axes=True, mri_fiducials=True, dig='fiducials')