# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 20:17:57 2022

@author: smellor
"""

import numpy as np
import pandas
import warnings
import copy
from scipy.signal import butter,filtfilt, find_peaks
from scipy.interpolate import interp1d
from scipy.spatial.transform import Rotation as R
import mne

# Classes for storing optitrack data
class Marker:
    """
    Class for storing single OptiTrack markers
    """
    def __init__(self, Position = [], Name = None):
        self.__position__ = Position
        self.__name__ = Name
        
    def appendPosition(self, newPos):
        if newPos.shape == (1, 3):
            newPos = np.transpose(newPos)
        if np.size(newPos,0) != 3:
            raise ValueError('Size of new position is incompatible. The row indicates x, y and z so should have 3 elements')
        self.__position__ = np.append(self.__position__, newPos, axis=1)
        
    def setPosition(self, newPos):
        self.__position__ = newPos
        
    def deletePosition(self):
        self.__position__ = []
        
    def getPosition(self):
        return self.__position__
    
    def getName(self):
        return self.__name__
        
class RigidBody(Marker):
    """
    Class for storing OptiTrack rigid bodies
    """
    def __init__(self, Position = [], Rotation = [], Name = None):
        self.__position__ = Position
        self.__rotation__ = Rotation
        self.__name__ = Name
        
    def appendRotation(self, newRot):
        if newRot.shape == (1,4):
            newRot = np.transpose(newRot)
        self.__rotation__ = np.append(self.__rotation__, newRot, axis=1)
        
    def setRotation(self, newRot):
        self.__rotation__ = newRot
        
    def deleteRotation(self):
        self.__rotation__ = [];
        
    def getRotation(self):
        return self.__rotation__
    
class OTobj:
    """
    Class for storing OptiTrack data.
    """
    
    def __init__(self, filename = None, rigidBodies = {}, markers = {}, time = None, recordingProperties = None):
        if filename is None and recordingProperties is None:
            self.rigidBodies = {}
            self.markers = {}
            self.time = []
            self.recordingProperties = {
                    'frameRate' : 120,
                    'startTime' : 0,
                    'totalFrames' : 0,
                    'rotationType' : 'Quaternion',
                    'lengthUnits' : 'm'
                    }
        elif filename is None and recordingProperties != None:
            self.rigidBodies = rigidBodies
            self.markers = markers
            self.recordingProperties = recordingProperties
            self.time = time
        elif filename!= None and recordingProperties is None:
            # Start by importing data
            data = pandas.read_csv(filename, header=[1,2,4,5])
            
            # Create rigid bodies
            self.rigidBodies = {}
            if any('Rigid Body' in s for s in data.columns.levels):
                RBidx = [i for i, x in enumerate(data.columns) if x[0] == 'Rigid Body']
                RBcount = int(np.round((len(RBidx))/8))
                
                for i in range(0, RBcount):
                    # Find Name
                    RBlabpos = RBidx[i*8]
                    name = data.columns[RBlabpos][1]
                    
                    # Find Position
                    pos = np.zeros((3, len(data['Rigid Body', name, 'Position', 'X'])))
                    pos[0,:] = data['Rigid Body', name, 'Position', 'X'].values
                    pos[1,:] = data['Rigid Body', name, 'Position', 'Y'].values
                    pos[2,:] = data['Rigid Body', name, 'Position', 'Z'].values
                    
                    # Find Rotation
                    rot = np.zeros((4, len(data['Rigid Body', name, 'Rotation', 'X'])))
                    rot[0,:] = data['Rigid Body', name, 'Rotation', 'X'].values
                    rot[1,:] = data['Rigid Body', name, 'Rotation', 'Y'].values
                    rot[2,:] = data['Rigid Body', name, 'Rotation', 'Z'].values
                    rot[3,:] = data['Rigid Body', name, 'Rotation', 'W'].values
                    
                    self.rigidBodies[name] = RigidBody(Position=pos, Rotation=rot, Name=name)
                    
            # Create Markers
            self.markers = {}
            if any('Marker' in s for s in data.columns.levels):
                Midx = [i for i, x in enumerate(data.columns) if x[0] == 'Marker']
                Mcount = int(np.round((len(Midx))/3))
                
                for i in range(0, Mcount):
                    # Find name
                    Mlabpos = Midx[i*3]
                    name = data.columns[Mlabpos][1]
                    
                    # Find position
                    pos = np.zeros((3, len(data['Marker', name, 'Position', 'X'])))
                    pos[0,:] = data['Marker', name, 'Position', 'X'].values
                    pos[1,:] = data['Marker', name, 'Position', 'Y'].values
                    pos[2,:] = data['Marker', name, 'Position', 'Z'].values
                    
                    self.markers[name] = Marker(Position = pos, Name = name)
            
            # Save time
            tidx = [i for i, x in enumerate(data.columns) if x[3] == 'Time (Seconds)']
            cols = data.keys()
            self.time = data[cols[tidx]].values
            
            # Then import recording info
            recordTab = pandas.read_csv(filename, nrows=1, header=None)
            self.recordingProperties = {
                    'frameRate' : recordTab[[i+1 for i, x in enumerate(recordTab.values[0]) if x == 'Capture Frame Rate']].values[0,0],
                    'startTime' : recordTab[[i+1 for i, x in enumerate(recordTab.values[0]) if x == 'Capture Start Time']].values[0,0],
                    'totalFrames' : recordTab[[i+1 for i, x in enumerate(recordTab.values[0]) if x == 'Total Exported Frames']].values[0,0],
                    'rotationType' : recordTab[[i+1 for i, x in enumerate(recordTab.values[0]) if x == 'Rotation Type']].values[0,0]
                    }
            if recordTab[[i+1 for i, x in enumerate(recordTab.values[0]) if x == 'Length Units']].values[0,0] == 'Millimeters':
                self.recordingProperties['lengthUnits'] = 'mm'
            elif recordTab[[i+1 for i, x in enumerate(recordTab.values[0]) if x == 'Length Units']].values[0,0] == 'Meters':
                self.recordingProperties['lengthUnits'] = 'm'
            else:
                raise Exception('Unrecognised Length Units')
                
    def listRB(self):
        """
        Function to list rigid bodies stored in optitrack object
        """
        rblist = list(self.rigidBodies.keys())
        return rblist
    
    def addRB(self, rigidBody):
        """
        Function to add a rigid body to the optitrack object
        """
        keys = self.listRB()
        if rigidBody.getName() in keys:
            raise Exception("Rigid Body with the same name already exists in OptiTrack object")
        else:
            if rigidBody.__position__.shape[1] == len(self.time):
                self.rigidBodies[rigidBody.getName()] = rigidBody
            else:
                raise Exception("Added rigid body position has different size to time save in OptiTrack object")
    
    def listMarkers(self):
        """
        Function to list all the markers in an optitrack object
        """
        markerlist = list(self.markers.keys())
        return markerlist
    
    def addMarker(self, marker):
        """
        Function to add a marker to the optitrack object  
        """
        keys = list(self.markers.keys())
        if marker.getName() in keys:
            raise Exception("Marker with same name already saved in OptiTrack object")
        else:
            if marker.__position__.shape[1] == len(self.time):
                self.markers[marker.getName()] = marker
            else:
                raise Exception("Added marker has different size to time as saved in OptiTrack object")
                
    def append(self, newOTobj):
        """
        Function to append two optitrack objects
        """
        
        # Check recording properties are compatible
        if self.recordingProperties['rotationType'] != newOTobj.recordingProperties['rotationType']:
            raise Exception("OptiTrack objects have different rotation types")
        
        # Get length of each object
        tlenSelf = len(self.time)
        tlenNew = len(newOTobj.time)
        
        if tlenSelf > 0:
            # Get marker list of self and new
            MLself = self.listMarkers()
            MLnew = newOTobj.listMarkers()
            # Create markers with NaNs in whichever object they are missing
            if sorted(MLself) != sorted(MLnew):
                # Print warning
                warnings.warn("Markers in original and appended optitrack objects do not match. Missing marker positions will be filled with NaNs")
                
                # Loop through each member of MLnew and save the names which are not in MLself
                missingMLself = []
                if len(MLself) != 0:
                    for i in MLnew:
                        for j in MLself:
                            if i == j:
                                break
                            elif j == len(MLself):
                                missingMLself.append(i)
                else:
                    missingMLself = MLnew
                missingMLnew = []
                if len(MLnew) != 0:
                    for i in MLself:
                        for j in MLnew:
                            if i == j:
                                break
                            elif j == len(MLnew):
                                missingMLnew.append(MLself[i])
                else:
                    missingMLnew = MLself
                
                # Create new markers
                for i in missingMLself:
                    self.addMarker(marker = Marker(Position=np.full((3, tlenSelf), np.nan), Name=i))
                for i in missingMLnew:
                    newOTobj.addMarker(marker = Marker(Position=np.full((3,tlenNew), np.nan), Name=i))
            
            # Append all markers
            MLself = self.listMarkers()
            for i in MLself:
                self.markers[i].appendPosition(newOTobj.markers[i].getPosition())
        
            # Get rigid body lists of self and new
            RBLself = self.listRB()
            RBLnew = newOTobj.listRB()
            if sorted(RBLself) != sorted(RBLnew):
                # Print warning
                warnings.warn("Rigid bodies in original and appended optitrack objects do not match. Missing rigid body positions and orientations will be filled with NaNs")
                
                # Loop through each member of RBLnew and save names which are not in RBLself
                missingRBLself = []
                for i in RBLnew:
                    for j in RBLself:
                        if RBLnew[i] == RBLself[j]:
                            break
                        elif j == len(MLself):
                            missingRBLself.append(RBLnew[i])
                missingRBLnew = []
                for i in RBLself:
                    for j in RBLnew:
                        if RBLself[i] == RBLnew[j]:
                            break
                        elif j == len(RBLnew):
                            missingRBLnew.append(RBLself[i])
                
                # Create new rigidbodies
                for i in missingRBLself:
                    self.addRB(rigidBody = RigidBody(Position = np.full((3, tlenSelf), np.nan), 
                                                     Rotation = np.full((3, tlenSelf), np.nan), Name = i))
                for i in missingRBLnew:
                    newOTobj.addRB(rigidBody = RigidBody(Position = np.full((3, tlenSelf), np.nan), 
                                                     Rotation = np.full((3, tlenSelf), np.nan), Name = i))
            # Append all rigid bodies
            RBself = self.listRB()
            for i in RBself:
                self.rigidBodies[i].appendPosition(newOTobj.rigidBodies[i].getPosition())
                self.rigidBodies[i].appendRotation(newOTobj.rigidBodies[i].getRotation())
                
            # Append Time
            self.time = np.append(self.time, newOTobj.time)
        else:
            self = copy.deepcopy(newOTobj)
        
        # Update recording properties
        self.recordingProperties['frameRate'] = np.nan
        self.recordingProperties['totalFrames'] = len(self.time)
        
    def filterOT(self, cutOffs, filtertype='lowpass'):
        """
        Filter optitrack data. Only for pre-recorded data - not a point-by-point,
        real-time filter.
        
        Parameters
        ----------
        cutOffs : double
            Cut off frequency/ies. One number of lowpass or highpass filter, 2 
            if bandpass or bandstop filter.
        filtertype : string
            Filter type. Options are: 'lowpass', 'highpass', 'bandpass' or 'bandstop'. 
            The default is 'lowpass'

        Returns
        -------
        None.

        """
        
        # Define filter
        b, a = butter(5, cutOffs, btype=filtertype, analog=False, fs = self.recordingProperties['frameRate'])
        
        # Rigid bodies
        RBself = self.listRB()
        for i in RBself:
            currentPos = self.rigidBodies[i].getPosition()
            currentRot = self.rigidBodies[i].getRotation()
            newPos = filtfilt(b, a, currentPos)
            newRot = filtfilt(b, a, currentRot)
            self.rigidBodies[i].setPosition(newPos)
            self.rigidBodies[i].setRotation(newRot)
        
        # Markers
        MLself = self.listMarkers()
        for i in MLself:
            currentPos = self.markers[i].getPosition()
            newPos = filtfilt(b, a, currentPos)
            self.markers[i].setPosition(newPos)
            
        
    def resample(self, tnew, method='linear'):
        """
        Resample optitrack data to fsample_new Hz

        Parameters
        ----------
        fsample_new : double
            new sampling frequency.
        method : string, optional
            Inerpolation method. Options match options for "kind" in 
            scipy.interpolate.interp1d. The default is 'linear'.

        Returns
        -------
        None.

        """

        # If new sampling freq is lower than original, filter data
        if 1/np.mean(np.diff(tnew)) < self.recordingProperties['frameRate']:
            warnings.warn('Optitrack data will automatically be filtered to half Nyquist limit of new sampling frequency')
            self.filterOT((1/np.mean(np.diff(tnew)))/4)
        
        # Get time arrays
        torig = np.squeeze(self.time)
        
        # Throw a warning message if tfinal and tinit are apart by more than 500 ms
        if (np.max(tnew) - np.max(torig)) > 0.5 or (np.min(tnew) - np.min(torig)) < -0.5:
            warnings.warn('Sampled times are not well matched; significant extrapolation will be required')
        
        # Rigid bodies
        for i in self.listRB():
            # Position
            currentPos = self.rigidBodies[i].getPosition()
            F = interp1d(torig, currentPos, kind=method, fill_value='extrapolate')
            newPos = F(tnew)
            self.rigidBodies[i].setPosition(newPos)
        
            # Rotation
            currentRot = self.rigidBodies[i].getRotation()
            F = interp1d(torig, currentRot, kind=method, fill_value='extrapolate')
            newRot = F(tnew)
            self.rigidBodies[i].setRotation(newRot)
        
        # Markers
        for i in self.listMarkers():
            # Position
            currentPos = self.markers[i].getPosition()
            F = interp1d(torig, currentPos, kind=method, fill_value='extrapolate')
            newPos = F(tnew)
            self.markers[i].setPosition(newPos)
        
        # Update time variable
        self.time = tnew
        
        # Update sampling frequency value
        self.recordingProperties['frameRate'] = 1/np.mean(np.diff(tnew))
        
        # Update total number of frames
        self.recordingProperties['totalFrames'] = len(tnew)
        
    def __str__(self):
        if len(self.time.shape) == 2:
            mint = self.time[0,0]
            maxt = self.time[-1,0]
        else:
            mint = self.time[0]
            maxt = self.time[-1]
        str = 'OptiTrack Data Info:\n' + ' Rigid bodies: ' + ', '.join(self.listRB()) + \
            '\n Markers: ' + ', '.join(self.listMarkers()) + '\n Time Range: ' + \
            '{:.2f} - {:.2f} s\n'.format(mint, maxt) + \
            ' Sampling Frequency: {} Hz\n Recording Start Time: {}\n Total Frames: {}\n Rotation Type: {}\n Length Units: {}'\
                .format(self.recordingProperties['frameRate'], self.recordingProperties['startTime'], \
                        self.recordingProperties['totalFrames'], self.recordingProperties['rotationType'], \
                            self.recordingProperties['lengthUnits'])
        return str
        
        
    def __removeTimePoints__(self, samplesToRemove):
        """
        Remove time points from optitrack data

        Parameters
        ----------
        samplesToRemove : integer, list of integers or slice
            Time points/samples to remove from optitrack data. 
            Recommend that they are all at either the start or end of the 
            recording, otherwise better to replace values with NaN to avoid 
            confusion with sampling frequency.

        Returns
        -------
        None.

        """
        
        # Create boolean mask
        torig = np.squeeze(self.time)
        mask = np.ones(len(torig), dtype=bool)
        mask[samplesToRemove] = False
        
        # Get time array
        tnew = torig[mask]
        
        # Rigid bodies
        for i in self.listRB():
            # Position
            currentPos = self.rigidBodies[i].getPosition()
            newPos = currentPos[:,mask]
            self.rigidBodies[i].setPosition(newPos)
        
            # Rotation
            currentRot = self.rigidBodies[i].getRotation()
            newRot = currentRot[:,mask]
            self.rigidBodies[i].setRotation(newRot)
        
        # Markers
        for i in self.listMarkers():
            # Position
            currentPos = self.markers[i].getPosition()
            newPos = currentPos[:,mask]
            self.markers[i].setPosition(newPos)
        
        # Update time variable
        self.time = tnew
        
        # Update total number of frames
        self.recordingProperties['totalFrames'] = len(tnew)
            
        
# Functions involving optitrack

# Synchronise optitrack and opm data
def sync_optitrack_and_opm(opti_data, opm_data, LengthsAlreadyMatch=False, TriggerChannelName='OptiTrig', \
                           Trigger = None, ResampleOPMdata=False):
    """
    Synchronise OPM and Optitrack recordings
    Trims and resamples each dataset. Should occur before epoching.
    
    The trigger is asssumed to be in the same time as the OPM data and by
    default the OptiTrack data is resampled to match the OPMs, rather than 
    the other way around, although this can be adjusted using the
    "ResampleOPMdata" input option.

    INTPUT:
    - opti_data: the OptiTrack recordings
    - opm_data options: MNE python meeg object containing OPM data.
    - Additional Options:
        - LengthsAlreadyMatch (default: false): Boolean to determine
            whether or not to trim the OptiTrack data to match a given
            trigger. When set to true, it is assumed that the start and end
            of the OptiTrack and OPM recordings are already syncronised so
            there so any trimming of the data will be skipped.
        - TriggerChannelName (default: 'OptiTrig'): If trimming of the data
            is required and the OPM data is a meeg or FieldTrip object, you
            can specify the name of the trigger channel here. This trigger
            must be high when the OptiTrack is recording and low when it is 
            off.
        - Trigger (default: None) Double array of zeros and ones for when the Optitrack is
            recording (1) or off (0). Particularly useful way to enter the
            Optitrack trigger if you either put the OPM data in as a matrix
            or recorded triggers in such a way that it wasn't a
            straightforward high-gated (optitrack on when high) trigger.
            N.B. the assumption is made that if only one change in the
            trigger is found, then this is the Optitrack turning on and the
            moment it turns off was missed by the OPMs, so the Movement
            data is also trimmed. A warning will be produced if this is the
            case.
        - ResampleOPMdata (default: false): Boolean to choose to resample
            the OPM data to match the OptiTrack data, rather than the other
            way around. All of the trimming will still match the previously
            given pattern (i.e. trim the OPM data unless only one trigger
            point is available).

    OUTPUT:
    - opti_data: Resampled (or untouched if ResampleOPMdata = true) 
        OptiTrack recordings, output in the same format as they are input
    - opm_data: Trimmed and/or resampled OPM data, output in the same
        format as it is input
                         
    """
    
    # Copy opm_data
    opm_data = opm_data.copy()
    
    # Trim data
    if not LengthsAlreadyMatch:
        # Find optitrack trigger
        if Trigger is None:
            Trigger = opm_data.get_data(picks=[TriggerChannelName])
        
        # Find steps
        Trigger = np.squeeze(Trigger)
        trigevs, props = find_peaks(np.abs(np.diff(Trigger)), height=0.5*np.max(np.abs(np.diff(Trigger))))
        trigt = opm_data.times[trigevs]
        
        # Check length of trigevs is 2. If only one step recorded, print a warning
        if len(trigevs) > 2:
            raise ValueError('Too many steps in the trigger were found. \nCheck that \
                             the optitrack trigger only steps up and down')
        elif len(trigevs) == 0:
            raise ValueError('No trigger steps found. Check given trigger channel')
        elif len(trigevs) == 1:
            warnings.warn("Only one step in the trigger was found. \n\
                          It will be assumed that this is from the start of the Optitrack recording")
            trigt = np.append(trigt, np.max(opm_data.times))
        elif len(trigevs) == 2:
            print('Start-End of OPM Trigger: {:.2f} s'\
                  .format(opm_data.times[trigevs[1]]-opm_data.times[trigevs[0]]))
        print('Length of Optitrack data: {:.2f} s'.format(np.ptp(opti_data.time)))
        
        # Trim OPM data
        opm_data.crop(tmin=trigt[0], tmax=trigt[1], include_tmax=True)
        
        # Trim movement data (if needed)
        if len(trigevs) == 1:
            topm = opm_data.times
            tmov = opti_data.time
            trim_mask = tmov > np.max(topm)
            opti_data.__removeTimePoints__(np.asarray(trim_mask).nonzero())
    else:
        trigevs = 0
        
    # Resample data
    if not ResampleOPMdata:
        opti_data.resample(opm_data.times)
    else:
        opm_data.resample(opti_data.recordingProperties['frameRate'])
        
    print('OPM and Optitrack Data should now be synchronised')
        
    return opti_data, opm_data, trigt

# Append rigid body position and orientation to opm data (useful for plotting and noise cancallation)
def append_RB_to_opm_data(opm_data, opti_data, RBname=None, Quaternion=False):
    """
    Append position and orientation of a given rigid body to the OPM data

    Parameters
    ----------
    opm_data : mne.io.Raw
        MNE python meeg object containing opm data.
    opti_data : OTobj
        The optitrack recordings. Must contain at least one rigid body.
    RBname : string, optional
        Name of the rigid body to append to the opm data. The default is the 
        first rigid body in the data (sensible only if there is only 1 rigid 
        body in the data).
    Quaternion : bool, optional
        If True, rigid body rotation is output in Quaternions. Else, it is 
        output in Euler coordinates. The default is False

    Returns
    -------
    opm_data_out : mne.io.Raw
        opm_data with rigid body data appended as misc channel.

    """
    # Copy opm_data
    opm_data_out = opm_data.copy()
    
    # Get rigid body name
    if RBname == None:
        RBname = opti_data.listRB()[0]
    
    # Get rigid body data into the right format
    rotation = opti_data.rigidBodies[RBname].getRotation()
    if Quaternion:
        chanlist = ['RB_' + RBname + '_' + s for s in ['posX', 'posY', 'posZ', 'rotX', 'rotY', 'rotZ', 'rotW']] 
        if opti_data.recordingProperties['rotationType'] != 'Quaternion':
            r = R.from_euler('XYZ', rotation.T)
            rotation = r.as_quat().T
    else:
        chanlist = ['RB_' + RBname + '_' + s for s in ['posX', 'posY', 'posZ', 'rotX', 'rotY', 'rotZ']]
        if opti_data.recordingProperties['rotationType'] == 'Quaternion':
            r = R.from_quat(rotation.T)
            rotation = r.as_euler('XYZ').T
    data = np.vstack([opti_data.rigidBodies[RBname].getPosition(), rotation])
    
    # Create mne info object
    info = mne.create_info(ch_names=chanlist,
                           ch_types='misc',
                           sfreq=opti_data.recordingProperties['frameRate'])

    # Create a new mne.io.Raw object to save the rigid body data in
    opti_data_raw = mne.io.RawArray(data, info)
    
    # Append opti_data_raw to opm_data_out
    opm_data_out.add_channels([opti_data_raw], force_update_info=True)
    
    return opm_data_out
