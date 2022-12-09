# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 11:34:57 2022

@author: goneill
"""


import os.path as op
import sys

import opym
import mne
import numpy as np
import matplotlib.pyplot as plt

#%% the main import and coreg

# import the data
data_root = op.abspath('D:/triaxial_empty_room')
data_bin = op.join(data_root,'sub-P01_ses-002_task-noise_run-001_meg.bin')
raw = opym.io.read_raw_ucl(data_bin)