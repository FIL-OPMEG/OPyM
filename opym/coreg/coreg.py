#!/usr/bin/python
"""
@author: goneill
"""

from copy import deepcopy

from numpy import eye
from numpy.linalg import inv

from mne import read_surface
from mne.coreg import Coregistration
from mne.viz import plot_alignment, set_3d_view
from mne.transforms import apply_trans, Transform
from mne.utils import warn
from mne._freesurfer import get_mni_fiducials

import os

def coreg_identity(info,subject=None,subjects_dir=None,cras_shift=True):
    
    # if the dataset has no fiducials in it, we will assume the traformation
    # head -> meg is identity (and maybe a cras shift). This function will also
    # return a corrected info with the fiducials in it
    
    # try and autoget the SUBJECTS_DIR, otherwise warn
    if subjects_dir == None:
        subjects_dir = _determine_subjects_dir()
        
    cras = _get_cras(subjects_dir, subject)
    xform = _cras2xform(cras)
    
    return xform
    

def coreg_headcast(info,subject=None,subjects_dir=None,cras_shift=True):
    
    # simple coreg script assuming a custom heacast was a made for paticipant
    # using the same MRI processed in freesurfer
    
    # try and autoget the SUBJECTS_DIR, otherwise warn
    if subjects_dir == None:
        subjects_dir = _determine_subjects_dir()
    
    dig = deepcopy(info['dig'])
        
    if dig == None:
        print('no fiducials found, tryng fallback options')
        xform = coreg_identity(info,subject=subject,
                               subjects_dir=subjects_dir,
                               cras_shift=cras_shift)
        
        print('estimating fiducials')
        dig = get_mni_fiducials(subject,subjects_dir)
        dig_meg = deepcopy(dig)
        it = inv(xform['trans']);
        
        for ii in range(len(dig)):
            dig_meg[ii]['r'] = apply_trans(it, dig[ii]['r'])
            dig_meg[ii]['coord_frame'] = 4; # FIFFV_COORD_HEAD
            
        with info._unlock():
             info['dig'] = dig_meg
        
    else:
        
        # if there are digitsation points, inverse transform from device to head
        # coordinate frame
        
        t = info['dev_head_t']['trans']
        it = inv(t);
        
        for ii in range(len(dig)):
            dig[ii]['r'] = apply_trans(it, dig[ii]['r'])
            dig[ii]['coord_frame'] = 4; # FIFFV_COORD_HEAD
            
        if cras_shift == True:
            cras = _get_cras(subjects_dir, subject)
            dig = _shift_dig_by_cras(dig, cras)
        
    coreg =  Coregistration(info, subject, subjects_dir, fiducials=dig)
    coreg.fit_fiducials(verbose=True)
    
    return coreg
    

def check_datareg(info,coreg=None,subject=None,subjects_dir=None,trans=None):
    
    # raise warning if no subject name specified
    if subject == None:
        warn('No subject specified, prepare for error!')
    
    # try and autoget the SUBJECTS_DIR, otherwise warn
    if subjects_dir == None:
        subjects_dir = _determine_subjects_dir()
    
    # arguments for visualisation
    if coreg == None:
        plot_kwargs = dict(subject=subject, subjects_dir=subjects_dir, dig=False,
                           eeg=[], meg='sensors', show_axes=True,
                           coord_frame='mri') 
    else:
        plot_kwargs = dict(subject=subject, subjects_dir=subjects_dir, dig=True,
                           eeg=[], meg='sensors', show_axes=True,
                           coord_frame='mri',mri_fiducials=coreg.fiducials.dig)
        
    view_kwargs = dict(azimuth=45, elevation=90, distance=0.6,
                       focalpoint=(0., 0., 0.))
    
    if coreg == None:
        fig = plot_alignment(_switch_coildef(info), trans=trans, **plot_kwargs)
    else:
        fig = plot_alignment(_switch_coildef(info), trans=coreg.trans, **plot_kwargs)
    set_3d_view(fig, **view_kwargs)
    
    
def plot_sens(info):
        
    # arguments for visualisation
    plot_kwargs = dict(subject=None,  dig=False,
                           eeg=[], meg='sensors', show_axes=True,
                           coord_frame='meg')
   
        
    view_kwargs = dict(azimuth=45, elevation=90, distance=0.6,
                       focalpoint=(0., 0., 0.))
    
    fig = plot_alignment(_switch_coildef(info), **plot_kwargs)
    set_3d_view(fig, **view_kwargs)
        
        

def _switch_coildef(info,old=8002,new=7002):
    # update the info structure for a coil definition which is better
    # for visualisation (default is the BabyMEG magnetometers, aww!)
    new_info = deepcopy(info)
    for ii in range(len(new_info['chs'])):
        if info['chs'][ii]['coil_type'] == old:
            new_info['chs'][ii]['coil_type'] = new
    
    return new_info

def _determine_subjects_dir():
    subjects_dir = os.getenv('SUBJECTS_DIR')
    if subjects_dir == None:
        raise ValueError('subjects_dir not specified!')
        
    return subjects_dir

def _shift_dig_by_cras(dig,cras):
    for ii in range(len(dig)):
        dig[ii]['r']=dig[ii]['r']-cras
    return dig

def _get_cras(subjects_dir,subject):
    tmp = read_surface(os.path.join(subjects_dir,subject,
                                   'bem','watershed',
                                   subject+'_brain_surface'),
                           read_metadata=True)
    cras = tmp[2]['cras']/1000
    return cras

def _cras2xform(cras):
    t = eye(4)
    # assuming cras is in m, no scaling needed
    t[0:3,-1] = -cras
    xform = Transform(5,4,t)
    return xform
