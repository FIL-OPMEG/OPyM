# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 14:21:18 2022

@author: goneill
"""
from copy import deepcopy

import numpy as np
import os.path as op
import scipy.io 

from .utils import _refine_sensor_orientation, _determine_position_units, _get_plane_vectors
from ..utils import root

from mne.io.constants import FIFF
from mne.io.meas_info import _empty_info
from mne.io.write import get_new_file_id
from mne.io.base import BaseRaw

def read_raw_cerca(binfile,slotnames=False,line_freq=50):
    return RawCerca(binfile,slotnames)

class RawCerca(BaseRaw):
    def __init__(self, binfile,slotnames=False,line_freq=50):
                
        Adim_conv = np.array((2**32, 2**16, 2**8, 1))

        nbytes = op.getsize(binfile)

        # determine size of each block and how many there might be
        blockhdr = np.fromfile(binfile,dtype='>u1',count=8)
        blockhdr = blockhdr.reshape((2,4))
        blockdim = Adim_conv @ blockhdr.T
        nelements = 1 + blockdim.prod()
        bytesperblock = 8*nelements
        nblocks = nbytes/bytesperblock

        # magic number = blockhdr bytes cast as a float
        magic = np.fromfile(binfile,dtype='>d',offset=0,count=1)

        dat = np.zeros((blockdim[0],blockdim[1],int(nblocks)))

        for ii in range(int(nblocks)):
            tmp = np.fromfile(binfile,dtype='>d',offset=ii*bytesperblock,count=nelements)
            assert tmp[0]==magic
            tmp = np.delete(tmp,0,0)
            tmp = tmp.reshape(blockdim)
            dat[:,:,ii] = tmp.copy();

        dat = dat.reshape(blockdim[0],-1,order='F')
        
        infofile = op.splitext(binfile)[0] + '_SessionInfo.txt'
        posfile = op.splitext(binfile)[0] + '_sensor_order.mat'
        
        sens = _cerca_get_sensor_names(infofile)
        cal = _cerca_get_cal(infofile)
        
        """Determine which helmet used"""
        if op.exists(posfile):
            helmet = _cerca_determine_helmet(infofile)
            pos, ori, sens = _cerca_get_pos(posfile,sens,helmet,update_names=slotnames)
        else:
            print('No helmet definition found!')
            pos = [];
            ori = [];

        """Create info structure"""
        info = _empty_info(_cerca_get_fsample(dat))

        # Collect all the necessary data from the structures read
        info['meas_id'] = get_new_file_id()
        tmp = _cerca_convert_channel_info(sens, cal, pos, ori)
        info['chs'] = _refine_sensor_orientation(tmp)
        info['line_freq'] = line_freq
        info._unlocked = False
        info._update_redundant()
        
        # purge time channel at this point
        dat = np.delete(dat,0,0)
        # rescale as this isnt working automatically ATM
        dat = _cerca_rescale_data(dat, info['chs'])
        
        
        super(RawCerca, self).__init__(
            info, preload=dat, last_samps=dat.shape[1]-1)     
    
def _cerca_get_sensor_names(data_info):
    assert op.exists(data_info)
    with open(data_info) as f:
       txt_info = f.readlines()
        
    # find line with sensor names
    for text in txt_info:
        if 'Sensor Names' in text:
            txt_sensors = text
    
    sensor_list = txt_sensors.split(',')
    tmp = sensor_list[0]
    tmp = tmp.split(':')
    sensor_list[0] = tmp[1];
    
    for ii in range(len(sensor_list)):
        sensor_list[ii] = sensor_list[ii].strip()
    
    return sensor_list

def _cerca_get_fsample(dat):
    fsample = round(1/(dat[0,1] - dat[0,0]))
    return fsample

def _cerca_get_cal(data_info):
    assert op.exists(data_info)
    with open(data_info) as f:
        txt_info = f.readlines()
        
    # find line with sensor names
    for text in txt_info:
        if 'OPM Gain' in text:
            txt_gain = text.split(':')[1]
    gain = float(txt_gain.split('x')[0])    

    # assuming data is converted from V to nT to T!
    cal = 1e-9/(2.7*gain)
    return cal

def _cerca_convert_channel_info(chans,cal,pos,ori):
    nmeg = nstim = nref = 0
    
    if pos is not None:
        units, sf = _determine_position_units(pos)
        
    chs = list()
    for ii in range(0,len(chans)):
        ch = dict(scanno=ii + 1, range=1., cal=1., loc=np.full(12, np.nan),
                  unit_mul=FIFF.FIFF_UNITM_NONE, ch_name=chans[ii],
                  coil_type=FIFF.FIFFV_COIL_NONE)    
       
        chs.append(ch)
        
        # create the channel information
        if sum(np.isnan(pos[ii])) == 0:
            r0 = pos[ii].copy()/sf # mm to m
            ez = ori[ii].copy()
            ez = ez/np.linalg.norm(ez)
            ex, ey = _get_plane_vectors(ez)
            ch['loc'] = np.concatenate([r0, ex, ey, ez])
                    
        if chans[ii][0:4] == 'Trig': # its a trigger!
            nstim += 1
            ch.update(logno=nstim, coord_frame=FIFF.FIFFV_COORD_UNKNOWN,
                      kind=FIFF.FIFFV_STIM_CH, unit=FIFF.FIFF_UNIT_V)
        elif sum(np.isnan(pos[ii])) == 3:       # its a sensor with no location info
            nref += 1
            ch.update(logno=nref, coord_frame=FIFF.FIFFV_COORD_UNKNOWN,
                      kind=FIFF.FIFFV_REF_MEG_CH, unit=FIFF.FIFF_UNIT_T,
                      coil_type=FIFF.FIFFV_COIL_QUSPIN_ZFOPM_MAG2,
                      cal=cal)
        else:                                   # its a sensor!
            nmeg += 1
            ch.update(logno=nmeg, coord_frame=FIFF.FIFFV_COORD_DEVICE,
                      kind=FIFF.FIFFV_MEG_CH, unit=FIFF.FIFF_UNIT_T,
                      coil_type=FIFF.FIFFV_COIL_QUSPIN_ZFOPM_MAG2,
                      cal=cal)

    return chs
    
def _cerca_rescale_data(dat,chs):
    
    calmat = np.zeros((dat.shape[0],1))
    
    for ii in range(dat.shape[0]):
        calmat[ii] = chs[ii]['cal']
        
    dat = calmat * dat
    return dat

def _cerca_determine_helmet(data_info):
    assert op.exists(data_info)
    with open(data_info) as f:
       txt_info = f.readlines()
        
    # find line with sensor names
    for text in txt_info:
        if 'Helmet' in text:
            txt_helmet = text.lower()
    
    if txt_helmet.find('orange') > -1:
        helmet = 'orange_adult_S'
    elif txt_helmet.find('purple') > -1:
        helmet = 'purple_adult_L'
    elif txt_helmet.find('light') > -1 and txt_helmet.find('blue'):
        helmet = 'light_blue_4YO'
    else:
        helmet = 'unknown'
    
    print('Helmet used: {0}'.format(helmet))
    return helmet

def _cerca_get_pos(posfile,sens,helmet,update_names=False):
    
    helmetinfo = _cerca_helmet_grabber(helmet)
    mat = scipy.io.loadmat(posfile,simplify_cells=True)
    garbage = mat['__function_workspace__'][0]
    nchans = len(helmetinfo['sens_labels'])
    chan_list = _cerca_garbage_compactor(garbage, nchans)
    
    pos = np.empty((len(sens),3))
    ori = np.empty((len(sens),3))
    slotno = []
    
    # identify which slot number each channel was in and get pos info
    for ii in range(len(sens)):
        prefix = sens[ii].split(' ')[0]
        pos[ii] = np.nan
        slotno.append(None)
        for jj, slot in enumerate(chan_list):
            if prefix in slot.split(' ')[0]:
                slotno[ii] = int(jj)
                pos[ii] = helmetinfo['sens_pos'][jj]
            
    # determine which axis and grab the ori informaiton
    for ii in range(len(sens)):
        ori[ii] = np.nan
        if slotno[ii] is not None:
            suffix = sens[ii].split(' ')[1]
            if suffix == '[X]':
                ori[ii] = helmetinfo['sens_ors_X'][slotno[ii]]
            elif suffix == '[Y]':
                ori[ii] = helmetinfo['sens_ors_Y'][slotno[ii]]
            elif suffix == '[Z]':
                ori[ii] = helmetinfo['sens_ors_Z'][slotno[ii]]
                
    # update the names to those of the slots if required
    channames = deepcopy(sens)
    if update_names is True:
        for ii, chan in enumerate(channames):
            if sum(np.isnan(pos[ii])) == 0:
                string = helmetinfo['sens_labels'][slotno[ii]] + ' ' + channames[ii].split(' ')[1]
                channames[ii] = string
                
    return pos, ori, channames
            


def _cerca_helmet_grabber(helmet):
    helmet_path = op.join(root(),'io','cerca_helmets',helmet+'.mat')
    mat = scipy.io.loadmat(helmet_path,simplify_cells=True)
    helmet_info = mat['Helmet_info']
    return helmet_info

def _cerca_garbage_compactor(data,nchans):
    # Somone decided to save out the slot2sensors information as a MATLAB
    # TABLE object, a close source class and no python compatibillity :(
    # HOWEVER
    # we know we are looking for a cell which is nslots x 1 in size,
    # so this means we can look for a magic number in the data stream (ugh).
    # see https://www.mathworks.com/help/pdf_doc/matlab/matfile_format.pdf
    # this is really bad code, but it works.
    
    target = np.array((0,0,0,1,
                       0,0,0,0,
                       0,0,0,5,
                       0,0,0,8,
                       0,0,0,nchans,
                       0,0,0,1),dtype='uint8')
    hits = []
    for ii in range(len(data)):
        if np.array_equal(data[ii:ii+len(target)], target):
            hits.append(ii)
    if len(hits) != 1:
        raise ValueError('couldnt locate the cell array which maps the sensors '
                         'to the helmet!')
    start = hits[0] + len(target)
    
    # rewind to the start of the first cell defintion ->  look for miMATRIX tag
    target = np.array((0,0,0,14),dtype='uint8')
    hits = -1
    while hits < 0:
        start += 1
        if np.array_equal(data[start:start+len(target)], target):
            hits = start
    start += 4;
    slot2sens = [];
    for ii in range(nchans):
        
        # roll forward 12 bytes
        start += 12;
        
        isChar = False
        
        if _int8tointX(data[start:start+4]) == 4:
            isChar = True
        
        start += 8
        string = ''
        if isChar is True:
            start += 4
            
            ndims = int(_int8tointX(data[start:start+4])/4)
        
            nels = 1;
            for jj in range(ndims):
                start += 4
                nels *= int(_int8tointX(data[start:start+4]))
            
            start += 8
            len_name = int(_int8tointX(data[start:start+4]))
            start += 8 + len_name + np.mod(len_name,4)
            assert(_int8tointX(data[start:start+4])==nels)
            start += 4
            
            while data[start] != 14:
                if data[start] > 31:
                    string += chr(data[start])
                start+=1
            start +=1
        else:
            while data[start] != 14:
                start+=1
            start +=1
        slot2sens.append(string)
    return slot2sens
    
    
def _int8tointX(arr):
    binary = np.unpackbits(arr)   
    return sum(val*(2**idx) for idx, val in enumerate(reversed(binary)))