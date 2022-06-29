# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 15:59:52 2022

@author: goneill


"""
import numpy as np

import mne
from mne.coreg import Coregistration
from mne.io import read_info
from opym import read_raw_ucl

import os.path as op
from mne.datasets import fetch_fsaverage


# Download fsaverage files
fs_dir = fetch_fsaverage(verbose=True)
subjects_dir = op.dirname(fs_dir)

# The files live in:
subject = 'gb'

# read in data
data_root = op.abspath('D:\masters_example_data\data')
data_bin = op.join(data_root,'sub-001_ses-001_task-motor4way_run-001_meg.bin')
raw = read_raw_ucl(data_bin)
info = raw.info
plot_kwargs = dict(subject=subject, subjects_dir=subjects_dir,
                   dig=True,
                   meg='sensors', show_axes=True,
                   coord_frame='meg')
view_kwargs = dict(azimuth=45, elevation=90, distance=0.6,
                   focalpoint=(0., 0., 0.))

fiducials = "estimated"  # get fiducials from fsaverage
coreg = Coregistration(info, subject, subjects_dir, fiducials=fiducials)
fig = mne.viz.plot_alignment(info, trans=coreg.trans, **plot_kwargs)


coreg.fit_fiducials(verbose=True)
fig = mne.viz.plot_alignment(info, trans=coreg.trans, **plot_kwargs)
