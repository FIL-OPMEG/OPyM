# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 11:34:57 2022

@author: goneill
"""


import os.path as op
import sys
sys.path.append('..\\') # Use this is running tests without a pip install!

import opym
import mne
import numpy as np
import matplotlib.pyplot as plt
from mne.beamformer import make_lcmv, apply_lcmv
#%% the main import and coreg

# import the data
data_root = op.abspath('D:\mixed_gen_noise')
data_bin = op.join(data_root,'sub-OP00042_ses-001_task-noise_run-001_meg.bin')
raw = opym.io.read_raw_ucl(data_bin)