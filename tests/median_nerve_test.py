# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 13:32:23 2022

@author: goneill
"""

import mne
import opym
import os.path as op
import copy


data_root = op.abspath('D:\old_mn_data')
data_bin = op.join(data_root,'sub-001_ses-001_task-mediannerve_run-001_meg.bin')


raw = opym.io.read_raw_ucl(data_bin)
info = raw.info;

# Lets have a look at trying to align data in the source space

subjects_dir = op.abspath('D:\\recons')
subject = 'masters-22-001'
# fs_dir = mne.datasets.fetch_fsaverage(verbose=True)
# subjects_dir = op.dirname(fs_dir)
# subject = 'fsaverage'

# set up the coreg object (this doesnt initialise it)
dig = copy.deepcopy(info['dig'])

# hack to shift by CRAS
lol = mne.read_surface(op.join(subjects_dir,subject,'bem','watershed',subject+'_brain_surface'),read_metadata=True)
cras = lol[2]['cras']/1000
for ii in range(3):
    dig[ii]['r']=dig[ii]['r']-cras

coreg =  mne.coreg.Coregistration(info, subject, subjects_dir,fiducials=dig)

# first alignment
plot_kwargs = dict(subject=subject, subjects_dir=subjects_dir, dig=True, eeg=[],
                   meg='sensors', show_axes=True,
                   coord_frame='mri',mri_fiducials=dig) #debugging with estimated
view_kwargs = dict(azimuth=45, elevation=90, distance=0.6,
                   focalpoint=(0., 0., 0.))


# do the initial fit with the fiducials
coreg.fit_fiducials(verbose=True)
# do the initial fit with the 


fig = mne.viz.plot_alignment(info, trans=coreg.trans, **plot_kwargs)
print(coreg.trans)