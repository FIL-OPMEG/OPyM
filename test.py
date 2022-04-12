# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 13:32:23 2022

@author: goneill
"""

import os.path as op
import tsv_handler as tsv
import json

data_root = op.abspath('C:/Users/goneill/Downloads/masters_22/data');
bids_root = 'sub-001_ses-001_task-motor4way_run-001'

files = dict();
files['bin'] = op.join(data_root,bids_root + '_meg.bin')
files['meg'] = op.join(data_root,bids_root + '_meg.json')
files['chans'] = op.join(data_root,bids_root + '_channels.tsv')
files['positions'] = op.join(data_root,bids_root + '_positions.tsv')
files['coordsystem'] = op.join(data_root,bids_root + '_coordsystem.json')

