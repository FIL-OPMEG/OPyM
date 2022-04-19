# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 13:29:53 2022

@author: goneill
"""

import numpy as np
from collections import OrderedDict

from mne.io.constants import FIFF
from mne.io.meas_info import _empty_info
from mne.io.write import get_new_file_id



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

def _convert_channel_info(chans):
    nmeg = nstim = nmisc = nref = 0
    chs = list()
    for ii in range(0,len(chans['name'])):
        ch = dict(scanno=ii + 1, range=1., cal=1., loc=np.full(12, np.nan),
                  unit_mul=FIFF.FIFF_UNITM_NONE, ch_name=chans['name'][ii],
                  coil_type=FIFF.FIFFV_COIL_NONE)    
       
        chs.append(ch)
       
        # create the channel information
        if chans['pos'][ii] is not None:
            r0 = chans['pos'][ii].copy() # mm to m
            ez = chans['ori'][ii].copy()
            ex, ey = _get_plane_vectors(ez)
            ch['loc'] = np.concatenate([r0, ex, ey, ez])
            
        if chans['type'][ii] == 'MEGMAG':
            nmeg += 1
            ch.update(logno=nmeg, coord_frame=FIFF.FIFFV_COORD_DEVICE,
                      kind=FIFF.FIFFV_MEG_CH, unit=FIFF.FIFF_UNIT_T,
                      coil_type=FIFF.FIFFV_COIL_QUSPIN_ZFOPM_MAG2)
        elif chans['type'][ii] == 'MEGREFMAG':
            nref += 1
            ch.update(logno=nref, coord_frame=FIFF.FIFFV_COORD_UNKNOWN,
                      kind=FIFF.FIFFV_REF_MEG_CH, unit=FIFF.FIFF_UNIT_T,
                      coil_type=FIFF.FIFFV_COIL_QUSPIN_ZFOPM_MAG2)
        elif chans['type'][ii] == 'TRIG':
            nstim += 1
            ch.update(logno=nstim, coord_frame=FIFF.FIFFV_COORD_UNKNOWN,
                      kind=FIFF.FIFFV_STIM_CH, unit=FIFF.FIFF_UNIT_V)
        else:
            nmisc += 1
            ch.update(logno=nmisc, coord_frame=FIFF.FIFFV_COORD_UNKNOWN,
                      kind=FIFF.FIFFV_MISC_CH, unit=FIFF.FIFF_UNIT_NONE)
    return chs

def _compose_meas_info(meg,chans):
    """Create info structure"""
    info = _empty_info(meg['SamplingFrequency'])
    
    # Collect all the necessary data from the structures read
    info['meas_id'] = get_new_file_id()
    info['chs'] = _convert_channel_info(chans)
    info['line_freq'] = meg['PowerLineFrequency']
    info['bads'] = _read_bad_channels(chans)
    info._unlocked = False
    info._update_redundant()
    return info

def _read_bad_channels(chans):
    bads = list()
    for ii in range(0,len(chans['status'])):
        if chans['status'][ii] == 'bad':
            bads.append(chans['name'][ii])
    return bads
                    
def _from_tsv(fname, dtypes=None):
    """Read a tsv file into an OrderedDict.
    Parameters
    ----------
    fname : str
        Path to the file being loaded.
    dtypes : list, optional
        List of types to cast the values loaded as. This is specified column by
        column.
        Defaults to None. In this case all the data is loaded as strings.
    Returns
    -------
    data_dict : collections.OrderedDict
        Keys are the column names, and values are the column data.
    """
    data = np.loadtxt(fname, dtype=str, delimiter='\t', ndmin=2,
                      comments=None, encoding='utf-8-sig')
    column_names = data[0, :]
    info = data[1:, :]
    data_dict = OrderedDict()
    if dtypes is None:
        dtypes = [str] * info.shape[1]
    if not isinstance(dtypes, (list, tuple)):
        dtypes = [dtypes] * info.shape[1]
    if not len(dtypes) == info.shape[1]:
        raise ValueError('dtypes length mismatch. Provided: {0}, '
                         'Expected: {1}'.format(len(dtypes), info.shape[1]))
    for i, name in enumerate(column_names):
        data_dict[name] = info[:, i].astype(dtypes[i]).tolist()
    return data_dict