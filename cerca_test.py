# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:07:13 2022

@author: goneill
"""


import numpy as np
import os.path as op
import matplotlib.pyplot as plt

from opym import read_raw_cerca
# from collections import OrderedDict

# from mne.io.constants import FIFF
# from mne.io.meas_info import _empty_info
# from mne.io.write import get_new_file_id
# from mne.io.base import BaseRaw
# from mne.io.utils import _read_segments_file
# from mne.io._digitization import _make_dig_points
# from mne.transforms import get_ras_to_neuromag_trans, apply_trans, Transform

# def _cerca_get_sensor_names(data_info):
#     assert op.exists(data_info)
#     with open(data_info) as f:
#        txt_info = f.readlines()
        
#     # find line with sensor names
#     for text in txt_info:
#         if 'Sensor Names' in text:
#             txt_sensors = text
    
#     sensor_list = txt_sensors.split(',')
#     tmp = sensor_list[0]
#     tmp = tmp.split(':')
#     sensor_list[0] = tmp[1];
    
#     for ii in range(len(sensor_list)):
#         sensor_list[ii] = sensor_list[ii].strip()
    
#     return sensor_list

# def _cerca_get_fsample(dat):
#     fsample = round(1/(dat[0,1] - dat[0,0]))
#     return fsample

# def _cerca_get_cal(data_info):
#     assert op.exists(data_info)
#     with open(data_info) as f:
#         txt_info = f.readlines()
        
#     # find line with sensor names
#     for text in txt_info:
#         if 'OPM Gain' in text:
#             txt_gain = text.split(':')[1]
#     gain = float(txt_gain.split('x')[0])    

#     # assuming data is converted from V to nT to T!
#     cal = 1e-9/(2.7*gain)
#     return cal

# def _cerca_convert_channel_info(chans,cal):
#     nmeg = nstim = 0
        
#     chs = list()
#     for ii in range(0,len(chans)):
#         ch = dict(scanno=ii + 1, range=1., cal=1., loc=np.full(12, np.nan),
#                   unit_mul=FIFF.FIFF_UNITM_NONE, ch_name=chans[ii],
#                   coil_type=FIFF.FIFFV_COIL_NONE)    
       
#         chs.append(ch)
                    
#         if chans[ii][0:4] == 'Trig': # its a trigger!
#             nstim += 1
#             ch.update(logno=nstim, coord_frame=FIFF.FIFFV_COORD_UNKNOWN,
#                       kind=FIFF.FIFFV_STIM_CH, unit=FIFF.FIFF_UNIT_V)
#         else:                      # its a sensor!     
#             nmeg += 1
#             ch.update(logno=nmeg, coord_frame=FIFF.FIFFV_COORD_DEVICE,
#                       kind=FIFF.FIFFV_MEG_CH, unit=FIFF.FIFF_UNIT_T,
#                       coil_type=FIFF.FIFFV_COIL_QUSPIN_ZFOPM_MAG2,
#                       cal=cal)
      
#     return chs
    

# precision = 'single';


data_root = op.abspath('D:\\20220623_095714_cMEG_Data')
data_bin = op.join(data_root,'truncated.cMEG')

rawC = read_raw_cerca(data_bin)
a = rawC[76,:]
plt.plot(a[1],a[0].transpose())

# Adim_conv = np.array((2**32, 2**16, 2**8, 1))

# nbytes = op.getsize(data_bin)

# # determine size of each block and how many there might be
# blockhdr = np.fromfile(data_bin,dtype='>u1',count=8)
# blockhdr = blockhdr.reshape((2,4))
# blockdim = Adim_conv @ blockhdr.T
# nelements = 1 + blockdim.prod()
# bytesperblock = 8*nelements
# nblocks = nbytes/bytesperblock

# # magic number = blockhdr bytes cast as a float
# magic = np.fromfile(data_bin,dtype='>d',offset=0,count=1)

# dat = np.zeros((blockdim[0],blockdim[1],int(nblocks)))

# for ii in range(int(nblocks)):
#     tmp = np.fromfile(data_bin,dtype='>d',offset=ii*bytesperblock,count=nelements)
#     assert tmp[0]==magic
#     tmp = np.delete(tmp,0,0)
#     tmp = tmp.reshape(blockdim)
#     dat[:,:,ii] = tmp.copy();

# dat = dat.reshape(blockdim[0],-1,order='F')

# sens = _cerca_get_sensor_names(data_info)
# cal = _cerca_get_cal(data_info)

# """Create info structure"""
# info = _empty_info(_cerca_get_fsample(dat))

# # Collect all the necessary data from the structures read
# info['meas_id'] = get_new_file_id()
# info['chs'] = _cerca_convert_channel_info(sens, cal)
# info['line_freq'] = 50.0
# info._unlocked = False
# info._update_redundant()


# # remove the time channel (ch 1)
# raww = mne.io.RawArray(np.delete(dat,0,0),info)