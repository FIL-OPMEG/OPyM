# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 13:32:23 2022

@author: goneill
"""

import os.path as op
import sys
sys.path.append('../') # Use this is running tests without a pip install!

import opym

data_root = op.abspath('D:\old_mn_data')
data_bin = op.join(data_root,'sub-001_ses-001_task-mediannerve_run-001_meg.bin')


raw = opym.io.read_raw_ucl(data_bin)
info = raw.info;

# Lets have a look at trying to align data in the source space
subjects_dir = op.abspath('D:\\recons')
subject = 'masters-22-001'


coreg = opym.coreg.coreg_headcast(info,subject=subject, subjects_dir=subjects_dir)
opym.coreg.check_datareg(info, coreg, subject=subject, subjects_dir=subjects_dir)