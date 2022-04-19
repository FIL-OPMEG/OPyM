# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 13:32:23 2022

@author: goneill
"""

from opym import _compose_meas_info
import numpy as np
import matplotlib.pyplot as plt
import os.path as op
import tsv_handler as tsv
import json

precision = 'single';

match precision:
    case 'single':
        dt = '>f'
    case 'double':
        dt = '>d'

data_root = op.abspath('C:/Users/goneill/Downloads/masters_22/data')
bids_root = 'sub-001_ses-001_task-motor4way_run-001'

files = dict();
files['bin'] = op.join(data_root,bids_root + '_meg.bin')
files['meg'] = op.join(data_root,bids_root + '_meg.json')
files['chans'] = op.join(data_root,bids_root + '_channels.tsv')
files['positions'] = op.join(data_root,bids_root + '_positions.tsv')
files['coordsystem'] = op.join(data_root,bids_root + '_coordsystem.json')

data = np.fromfile(files['bin'],dt)
chans = tsv._from_tsv(files['chans'])
chanpos = tsv._from_tsv(files['positions'])
nchans = len(chans['name'])
nlocs = len(chanpos['name'])
nsamples = len(data)

data.shape = (int(nsamples/nchans),nchans)

# add some additional channel information
# chans['type_mne'] = list()
# for ii in range(0,nchans):
#     match chans['type'][ii]:
#         case 'MEGMAG':
#             chans['type_mne'].append('mag')
#         case 'MEGREFMAG':
#             chans['type_mne'].append('ref_meg')
#         case 'TRIG':
#             chans['type_mne'].append('misc')
#         case _:
#             chans['type_mne'].append('misc')

chans['pos'] = [None] * nchans
chans['ori'] = [None] * nchans
for ii in range(0,nlocs):
    idx = chans['name'].index(chanpos['name'][ii])
    tmp = np.array([chanpos['Px'][ii], chanpos['Py'][ii], chanpos['Pz'][ii]])
    chans['pos'][idx] = tmp.astype(np.float64)/100
    tmp = np.array([chanpos['Ox'][ii], chanpos['Oy'][ii], chanpos['Oz'][ii]])
    chans['ori'][idx] = tmp.astype(np.float64)
    
fid = open(files['meg'],'r')
meg = json.load(fid)
#fid = open(files['coordsystem'],'r')
#coord = json.load(fid)
# chans = 

info = _compose_meas_info(meg,chans)
