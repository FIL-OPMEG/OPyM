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
    elif tmpname[seedID][-3:] == 'Z':
        prefix1 = 'Z'
        prefix2 = 'Y'
    elif tmpname[seedID][-3:] == 'Y':
        prefix1 = 'Y'
        prefix2 = 'Z'
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