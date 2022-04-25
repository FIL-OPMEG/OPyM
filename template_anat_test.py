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
