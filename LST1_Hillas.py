#!/usr/bin/env python3

"""
Example of Hillas parameters calculation from simtelarray file, and writing out ntuple as and astropy table
"""

import matplotlib.pylab as plt
import astropy.units as u
import numpy as np
import ctapipe
import os
from ctapipe.core import Container, Field, Map
from ctapipe.instrument import CameraGeometry
from ctapipe.visualization import CameraDisplay
from ctapipe.image import hillas_parameters, hillas_parameters_2, tailcuts_clean
from ctapipe.io.eventsourcefactory import EventSourceFactory
from ctapipe.image.charge_extractors import LocalPeakIntegrator
from astropy.visualization import quantity_support
from astropy.table import Table
from astropy.io import fits
import copy

class EventContainer(Container):
    event = Field(Map(),"Event")


if __name__ == '__main__':

    level1 = {'LSTCam' : 6.}
    level2 = level1.copy()
    # We use as second cleaning level just half of the first cleaning level
    for key in level2:
        level2[key] *= 0.5
    print(level2)


    source01 = EventSourceFactory.produce(input_url="/scratch/bernardos/LST1/gamma_20deg_180deg_run2200___cta-prod3-demo-2147m-LaPalma-demo3-FC_cone8.simtel.gz",allowed_tels={1});
#    source02 = EventSourceFactory.produce(input_url="/scratch/bernardos/LST1/gamma_20deg_180deg_run2201___cta-prod3-demo-2147m-LaPalma-demo3-FC_cone8.simtel.gz",allowed_tels={1});
    camtype = []   # one entry per image
    width = np.array([])
    length = np.array([])
    size = np.array([])
    
    #Event Parameters
    ObsID = np.array([])
    EvID = np.array([])
    
    #MC Parameters:
    mcEnergy = np.array([])
    mcAlt  = np.array([])
    mcAz = np.array([])
    mcCore_x = np.array([])
    mcCore_y = np.array([])
    mcHfirst = np.array([])
    mcType = np.array([])

    fitsdata = np.array([])

    log10pixelHGsignal = {}
    survived = {}

    ev = EventContainer()

    n=0;
    for event in source01:
        ev.event[n] = ctapipe.io.containers.DataContainer()
        ev.event[n] = copy.deepcopy(event)
        n=n+1
    '''''
    for event in source02:
        ev.event[n] = ctapipe.io.containers.DataContainer()
        ev.event[n] = copy.deepcopy(event)
        n=n+1
    '''''
    for key in level1:

        log10pixelHGsignal[key] = []
        survived[key] = []
    i=0
    for i in range(0,n-1):
        event = ev.event[i]
        
        if i%100==0:
            print("EVENT_ID: ", event.r0.event_id, "TELS: ",
                  event.r0.tels_with_data,
                  "MC Energy:", event.mc.energy )

        ntels = len(event.r0.tels_with_data)

        '''
        if i > 100:   # for quick tests
            break
        '''
        for ii, tel_id in enumerate(event.r0.tels_with_data):
            
            geom = event.inst.subarray.tel[tel_id].camera

            if str(geom) == 'FlashCam':  # FlashCam now skipped because of no proper pulse extraction
                continue
            elif str(geom) == 'ASTRI':   # excluded because no proper simulation in Prod3
                continue
            
            data = event.r0.tel[tel_id].waveform
            ped = event.mc.tel[tel_id].pedestal
            # the pedestal is the average (for pedestal events) of the *sum* of all samples, from sim_telarray

            nsamples = data.shape[2]  # total number of samples
            pedcorrectedsamples = data - np.atleast_3d(ped)/nsamples    # Subtract pedestal baseline. atleast_3d converts 2D to 3D matrix

            integrator = LocalPeakIntegrator(None, None)
            integration, peakpos, window = integrator.extract_charge(pedcorrectedsamples) # these are 2D matrices num_gains * num_pixels

            chan = 0  # high gain used for now...
            signals = integration[chan].astype(float)

            dc2pe = event.mc.tel[tel_id].dc_to_pe   # numgains * numpixels
            signals *= dc2pe[chan]

            # Add all individual pixel signals to the numpy array of the corresponding camera inside the log10pixelsignal dictionary
            log10pixelHGsignal[str(geom)].extend(np.log10(signals))  # This seems to be faster like this, with normal python lists

            # Apply image cleaning
            cleanmask = tailcuts_clean(geom, signals, picture_thresh=level1[str(geom)],
                                       boundary_thresh=level2[str(geom)], keep_isolated_pixels=False, min_number_picture_neighbors=1)
            survived[str(geom)].extend(cleanmask)  # This seems to be faster like this, with normal python lists

            clean = signals.copy()
            clean[~cleanmask] = 0.0   # set to 0 pixels which did not survive cleaning
            if np.max(clean) < 1.e-6: # skip images with no pixels
                continue



            # Calculate image parameters
            hillas = hillas_parameters(geom, clean) # this one gives some warnings invalid value in sqrt            
            foclen = event.inst.subarray.tel[tel_id].optics.equivalent_focal_length

            w = np.rad2deg(np.arctan2(hillas.width,foclen));
            l = np.rad2deg(np.arctan2(hillas.length,foclen));

            camtype.append(str(geom))
            width = np.append(width, w.value)
            length = np.append(length, l.value)
            size = np.append(size, hillas.size)


            #Store parameters from event and MC:
            ObsID = np.append(ObsID,event.r0.obs_id)
            EvID = np.append(EvID,event.r0.event_id)
            
            mcEnergy = np.append(mcEnergy,event.mc.energy)
            mcAlt = np.append(mcAlt,event.mc.alt)
            mcAz = np.append(mcAz,event.mc.az)
            mcCore_x = np.append(mcCore_x,event.mc.core_x)
            mcCore_y = np.append(mcCore_y,event.mc.core_y)
            mcHfirst = np.append(mcHfirst,event.mc.h_first_int)
            mcType = np.append(mcType,event.mc.shower_primary_id)
            

    #print(np.shape(camtype),np.shape(ObsID),np.shape(EvID),np.shape(mcEnergy),np.shape(mcAlt),np.shape(mcAz),np.shape(mcCore_x),np.shape(mcCore_y),np.shape(mcHfirst),np.shape(mcType),np.shape(width),np.shape(length),np.shape(size))

    output = {'camtype':camtype,'ObsID':ObsID,'EvID':EvID,'mcEnergy':mcEnergy,'mcAlt':mcAlt,'mcAz':mcAz, 'mcCore_x':mcCore_x,'mcCore_y':mcCore_y,'mcHfirst':mcHfirst,'mcType':mcType, 'width':width, 'length':length, 'size':size}
    ntuple = Table(output)

    if os.path.isfile('events.txt'):
        with open('events.txt',mode='a') as f:
            f.seek(0,os.SEEK_END)
            ntuple.write(f,format='ascii.no_header')

    else:
        ntuple.write('events.txt',format='ascii')

        
    table = Table.read('events.txt',format='ascii')
    table.write('events.fits',overwrite=True)

    '''
    hdu = fits.BinTableHDU.from_columns(
        [fits.Column(name='camtype',format = '20A', array=camtype),
        fits.Column(name='ObsID',format='E',array=ObsID),
        fits.Column(name='EvID',format='E',array=EvID),
        fits.Column(name='mcEnergy',format='E',array=mcEnergy),
        fits.Column(name='mcAlt',format='E',array=mcAlt),
        fits.Column(name='mcAz',format='E',array=mcAz),
        fits.Column(name='mcCore_x',format='E',array=mcCore_x),
        fits.Column(name='mcCore_y',format='E',array=mcCore_y),
        fits.Column(name='mcHfirst',format='E',array=mcHfirst),
        fits.Column(name='mcType',format='E',array=mcType),
        fits.Column(name='width',format='E',array=width),
        fits.Column(name='length',format='E',array=length),
        fits.Column(name='size',format='E',array=size)])

    hdu.writeto("out.fits",overwrite=True)
    ntuple.write("out.fits",format='fits',overwrite=True)
    
    fitsfile = 'out.fits'
    hdul = fits.open(fitsfile)
    hdr = hdul[1].header
    hdu.append(hdu)
    hdu.writeto("new.fits",overwrite=True)
    '''
    '''
    for key in level1:
        if level1[key] > 0:
            Table({'log10HGsignal': np.array(log10pixelHGsignal[key]), 'survivedcleaning' : np.array(survived[key])}).write("SPS_%s.fits" % str(key), overwrite=True)
    '''
