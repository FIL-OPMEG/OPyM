# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 13:29:53 2022

@author: goneill
"""

import json
import numpy as np
import os.path as op
from collections import OrderedDict

from mne.io.constants import FIFF
from mne.io.meas_info import _empty_info
from mne.io.write import get_new_file_id
from mne.io.base import BaseRaw
from mne.io.utils import _read_segments_file
from mne.io._digitization import _make_dig_points
from mne.transforms import get_ras_to_neuromag_trans, apply_trans, Transform


def read_raw_ucl(binfile, precision='single', preload=False):

    return RawUCL(binfile, precision=precision, preload=preload)

class RawUCL(BaseRaw):
    def __init__(self, binfile, precision='single', preload=False):
        
        if precision == 'single':
            dt = '>f'
            bps = 4
        else:
            dt = '>d'
            bps = 8
            
       
        sample_info = dict()
        sample_info['dt'] = dt
        sample_info['bps'] = bps
        
        files = _get_file_names(binfile)
    
        chans = _from_tsv(files['chans'])
        chanpos = _from_tsv(files['positions'])
        nchans = len(chans['name'])
        nlocs = len(chanpos['name'])
        nsamples = _determine_nsamples(files['bin'], nchans, precision) - 1
        sample_info['nsamples'] = nsamples
        
        raw_extras = list()
        raw_extras.append(sample_info)
        
        chans['pos'] = [None] * nchans
        chans['ori'] = [None] * nchans       
        
        for ii in range(0,nlocs):
            idx = chans['name'].index(chanpos['name'][ii])
            tmp = np.array([chanpos['Px'][ii], chanpos['Py'][ii], chanpos['Pz'][ii]])
            chans['pos'][idx] = tmp.astype(np.float64)
            tmp = np.array([chanpos['Ox'][ii], chanpos['Oy'][ii], chanpos['Oz'][ii]])
            chans['ori'][idx] = tmp.astype(np.float64)
            
        fid = open(files['meg'],'r')
        meg = json.load(fid)
        fid.close()
        info = _compose_meas_info(meg, chans)
        
        super(RawUCL, self).__init__(
            info, preload, filenames=[files['bin']],raw_extras=raw_extras,
            last_samps=[nsamples], orig_format=dt)
        
        if op.exists(files['coordsystem']):
            fid = open(files['coordsystem'],'r')
            csys = json.load(fid)
            fid.close()
            hc = csys['HeadCoilCoordinates']
            
            for key in hc:
                if key == 'lpa' or key == 'LPA':
                    lpa = np.asarray(hc[key])
                elif key == 'rpa' or key == 'RPA':
                    rpa = np.asarray(hc[key])
                elif key == 'nas' or key == 'NAS' or key == 'nasion':
                    nas = np.asarray(hc[key])
                    
            siz = np.linalg.norm(nas - rpa)
            unit, sf = _size2units(siz) 
            lpa/=sf
            rpa/=sf
            nas/=sf
            
            # t = get_ras_to_neuromag_trans(nas, lpa, rpa)
            # with self.info._unlock():
            #     self.info['dev_head_t'] = \
            #         Transform(FIFF.FIFFV_COORD_DEVICE,
            #                   FIFF.FIFFV_COORD_HEAD, t)

            # # transform fiducial points
            # nas = apply_trans(t, nas)
            # lpa = apply_trans(t, lpa)
            # rpa = apply_trans(t, rpa)
            
            with self.info._unlock():
                self.info['dig'] = _make_dig_points(nasion=nas,
                                                    lpa=lpa,
                                                    rpa=rpa)
                
                
            
            
                
        
    def _read_segment_file(self, data, idx, fi, start, stop, cals, mult):
        """Read a chunk of raw data."""
        si = self._raw_extras[fi]
        _read_segments_file(
            self, data, idx, fi, start, stop, cals, mult, dtype=si['dt'])
    


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
    
    units, sf = _determine_position_units(chans['pos'])
    
    chs = list()
    for ii in range(0,len(chans['name'])):
        ch = dict(scanno=ii + 1, range=1., cal=1., loc=np.full(12, np.nan),
                  unit_mul=FIFF.FIFF_UNITM_NONE, ch_name=chans['name'][ii],
                  coil_type=FIFF.FIFFV_COIL_NONE)    
       
        chs.append(ch)
       
        # create the channel information
        if chans['pos'][ii] is not None:
            r0 = chans['pos'][ii].copy()/sf # mm to m
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
            
            
        # set the calibration based on the units - MNE expects T units for meg
        # and V for eeg
        if chans['units'][ii] == 'fT':
            ch.update(cal=1e-15)
        elif chans['units'][ii] == 'pT':
            ch.update(cal=1e-12)
        elif chans['units'][ii] == 'nT':
            ch.update(cal=1e-9)
        elif chans['units'][ii] == 'mV':
            ch.update(cal=1e3)
        elif chans['units'][ii] == 'uV':
                ch.update(cal=1e6)

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

def _determine_nsamples(bin_fname,nchans,precision):
    bsize = op.getsize(bin_fname)
    if precision == 'single':
        bps = 4
    else:
        bps = 8
    nsamples = int(bsize/(nchans*bps))
    return nsamples

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

def _get_file_names(binfile):
    
    files = dict();
    files['dir'] = op.dirname(binfile)

    tmp = op.basename(binfile)
    tmp = str.split(tmp,'_meg.bin')

    files['root'] = tmp[0];
    files['bin'] = op.join(files['dir'],files['root'] + '_meg.bin')
    files['meg'] = op.join(files['dir'],files['root'] + '_meg.json')
    files['chans'] = op.join(files['dir'],files['root'] + '_channels.tsv')
    files['positions'] = op.join(files['dir'],files['root'] + '_positions.tsv')
    files['coordsystem'] = op.join(files['dir'],files['root'] + '_coordsystem.json')
    
    return files

def _determine_position_units(pos):
    
    # get rid of None elements
    nppos = np.empty((0,3))
    for ii in range(0,len(pos)):
        if pos[ii] is not None:
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
    
    