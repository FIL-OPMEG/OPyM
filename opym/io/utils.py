# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 14:15:21 2022

@author: goneill
"""

from copy import deepcopy
import numpy as np

def _refine_sensor_orientation(chanin):
    # the ex and ey elements from _convert_channel_info were randomly oriented
    # but it doesnt have to be this way, we can use (if available) the
    # orientation information from mulit-axis recordings to refine this.
    # THIS NEEDS SOME WORK TO MAKE BULLETPROOF
    print('refining sensor orientations')
    chanout = deepcopy(chanin)
    tmpname = list()
    for ii in range(len(chanin)):
        tmpname.append(chanin[ii]['ch_name'])
    
    for ii in range(len(chanin)):
        tmploc = deepcopy(chanin[ii]['loc']);
        tmploc = tmploc.reshape(3,4,order='F');
        if np.isnan(tmploc.sum()) == False:
            target = _guess_other_chan_axis(tmpname, ii)
            if np.isnan(target) == False:
                targetloc = deepcopy(chanin[target]['loc']);
                if np.isnan(targetloc.sum()) == False:
                    targetloc = targetloc.reshape(3,4,order='F');
                    tmploc[:,2] = targetloc[:,3]
                    tmploc[:,1] = np.cross(tmploc[:,2],tmploc[:,3])
                    chanout[ii]['loc'] = tmploc.reshape(12,order='F')
                    
        
    return chanout

def _guess_other_chan_axis(tmpname,seedID):
    # mad script which tries to guess the name of another channel which is
    # from the same sensor, but another axis
    targetID = np.NAN
    
    # see if its using the old RAD/TAN convention first
    if tmpname[seedID][-3:] == 'RAD':
        prefix1 = 'RAD'
        prefix2 = 'TAN'
    elif tmpname[seedID][-3:] == 'TAN':
        prefix1 = 'TAN'
        prefix2 = 'RAD'
    elif tmpname[seedID][-1:] == 'Z' or tmpname[seedID][-3:] == '[Z]':
        prefix1 = 'Z'
        prefix2 = 'Y'
    elif tmpname[seedID][-1:] == 'Y' or tmpname[seedID][-3:] == '[Y]':
        prefix1 = 'Y'
        prefix2 = 'Z'
    elif tmpname[seedID][-1:] == 'X' or tmpname[seedID][-3:] == '[X]':
        prefix1 = 'X'
        prefix2 = 'Y'
    else:
        prefix1 = '?'
        prefix2 = '?'
        
    targetName = tmpname[seedID][:-len(prefix1)] + prefix2
    
    # iterate through loop to find target, stop when found
    hit = False;
    ii = 0;
    while (ii < len(tmpname)) and (hit == False):
        if targetName == tmpname[ii]:
            targetID = ii;
            hit = True
        ii += 1
    
    return targetID

def _determine_position_units(pos):
    
    # get rid of None elements
    nppos = np.empty((0,3))
    for ii in range(0,len(pos)):
        if pos[ii] is not None and sum(np.isnan(pos[ii])) == 0:
            nppos = np.vstack((nppos, pos[ii]))
    
    idrange = np.empty(shape=(0,3))        
    for ii in range(0,3):
        q90, q10 = np.percentile(nppos[:,ii], [90 ,10])
        idrange= np.append(idrange,q90 - q10)
        
    siz = np.linalg.norm(idrange)
    
    unit, sf = _size2units(siz)
    
    return unit, sf
   
    
def _size2units(siz):
    
    if siz >= 0.050 and siz < 0.500:
        unit = 'm'
        sf = 1;
    elif siz >= 0.50 and siz < 5:
        unit = 'dm'
        sf = 10;
    elif siz >= 5 and siz < 50:
        unit = 'cm'
        sf = 100;
    elif siz >= 50 and siz < 500:
        unit = 'mm'
        sf = 1000;
    else:
        unit = 'unknown'
        sf = 1;
        
    return unit, sf

def _get_plane_vectors(ez):
    """Get two orthogonal vectors orthogonal to ez (ez will be modified)."""
    assert ez.shape == (3,)
    ez_len = np.sqrt(np.sum(ez * ez))
    if ez_len == 0:
        raise RuntimeError('Zero length normal. Cannot proceed.')
    if np.abs(ez_len - np.abs(ez[2])) < 1e-5:  # ez already in z-direction
        ex = np.array([1., 0., 0.])
    else:
        ex = np.zeros(3)
        if ez[1] < ez[2]:
            ex[0 if ez[0] < ez[1] else 1] = 1.
        else:
            ex[0 if ez[0] < ez[2] else 2] = 1.
    ez /= ez_len
    ex -= np.dot(ez, ex) * ez
    ex /= np.sqrt(np.sum(ex * ex))
    ey = np.cross(ez, ex)
    return ex, ey