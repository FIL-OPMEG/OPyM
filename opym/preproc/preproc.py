# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 14:44:00 2022

@author: goneill
"""

from mne.io.pick import pick_types
from mne.io.proj import Projection, activate_proj
from mne.io.compensator import get_current_comp
from numpy import zeros, isnan, roll, eye, arange, shape, array, newaxis
from scipy import linalg


def denoise_hfc(raw,L=1,copy=True,method='old'):
    
    if L > 1:
        raise ValueError('harmonic order higher than (L=1) not yet supported')
    else:
        basis, mag_inds = _get_homogenous_field_basis(raw.info)



    if method == 'old':
        # make the projector
        nchans = len(mag_inds)
        M = eye(nchans) - basis @ linalg.pinv(basis)
           
        if copy:
            raw = raw.copy()
            if not raw.preload:
                print('HFC: Loading raw data from disk')
                raw.load_data(verbose=False)
            else:
                print('HFC: Using loaded raw data')
      
        starts, stops = _get_chunks(raw)
        for start, stop in zip(starts, stops):
            orig_data = raw._data[mag_inds, start:stop]
            raw._data[mag_inds, start:stop] = M @ orig_data
    else:
        
        
        if copy:
            raw = raw.copy()
            
        l = (1,1,1)
        m = (-1,0,1)
        names = [];
        projs = [];
        for ii in mag_inds:
            names.append(raw.info['chs'][ii]['ch_name'])
        for ii in range(len(l)):
            data = basis[:,ii]
            proj_data = dict(col_names=names, row_names=None,
                             data=data[newaxis,:], ncol =len(names), nrow=1)
            this_desc = "HFC-l=%02d-m=%02d" % (l[ii], m[ii])
            # logger.info("Adding projection: %s" % this_desc)
            proj = Projection(
                active=True, data=proj_data, desc=this_desc)
            projs.append(proj)

        raw.add_proj(projs, remove_existing=True)
        # projs = raw.info['projs']
        # projs = activate_proj(projs)
        
        # with raw.info._unlock():
        #     raw.info['projs'] = activate_proj(raw.info['projs'])
        raw.info.normalize_proj()
        
        
    return raw
       
    
def _get_homogenous_field_basis(info):
    
    # easily defined basis set which is based on just the orientation of
    # the sensor (rather than fitting spherical harmonics)
    
    mag_ind = pick_types(info, meg='mag', ref_meg=False, exclude='bads')
    
    basis = zeros((len(mag_ind),3))
    for ii in range(len(mag_ind)):
        loc = info['chs'][mag_ind[ii]]['loc']
        o = loc[-3:]
        if isnan(sum(o)) == False:
            basis[ii,:] = o
            
    # roll the array to be in y-z-x which is the order of
    # m = [-1, 0, 1] when l = 1;
    basis = roll(basis,-1,axis=1)
    
    return basis, mag_ind

def _get_chunks(raw,size=512):
    # get blocks of data in a specific size (default 512 MB), return the beginning
    # and end samples to do this
    chunk_samples = round(size*1e6/(8*len(raw.info['chs'])))
    
    start = arange(0,len(raw),chunk_samples)
    stop = zeros(shape(start),dtype=int)
    if len(start) == 1:
        pass
    else:
        for ii in range(1,len(start)):
            stop[ii-1] = start[ii]-1
    stop[-1] = len(raw)-1   
    
    return start, stop
        