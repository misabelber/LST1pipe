import matplotlib.pylab as plt
import astropy.units as u
import numpy as np
import ctapipe
from ctapipe.core import Container, Field, Map
from ctapipe.instrument import CameraGeometry
from ctapipe.visualization import CameraDisplay
from ctapipe.image import hillas_parameters, hillas_parameters_2, tailcuts_clean
from ctapipe.io.eventsourcefactory import EventSourceFactory
from ctapipe.image.charge_extractors import LocalPeakIntegrator
from astropy.visualization import quantity_support
from astropy.table import Table
import copy

class EventContainer(Container):
    event = Field(Map(),"Event")

#source01 = EventSourceFactory.produce(input_url="/scratch/bernardos/LST1/gamma_20deg_0deg_run5805___cta-prod3-lapalma3-2147m-LaPalma_cone10.simtel.gz",allowed_tels={1,2,3,4});

source02 = EventSourceFactory.produce(input_url="/scratch/bernardos/LST1/gamma_20deg_180deg_run805___cta-prod3-merged_desert-2150m-Paranal-subarray-1-nosct.simtel.gz",allowed_tels={1,2,3,4});

camtype = []   # one entry per image                                                           
width = np.array([])                                                                           
length = np.array([])                                                                          
size = np.array([])                                                                            

log10pixelHGsignal = {}                                                                        
survived = {}


ev = EventContainer()

n=0;
#for event in source01:
    
#    ev.event[n] = ctapipe.io.containers.DataContainer()
#    ev.event[n] = copy.deepcopy(event)
#    n=n+1
#print(event.count)

for event in source02:
    ev.event[n] = ctapipe.io.containers.DataContainer()
    ev.event[n] = copy.deepcopy(event)
    n=n+1
i=0;
for i in range(0,n-1):
    event = ev.event[i]
    print("EVENT_ID: ", event.r0.event_id, "TELS: ",
          event.r0.tels_with_data,
          "MC Energy:", event.mc.energy )

    ntels = len(event.r0.tels_with_data)
    for ii, tel_id in enumerate(event.r0.tels_with_data):
        # print("\t cam {}...".format(tel_id))
        geom = event.inst.subarray.tel[tel_id].camera
        # print(event.inst.subarray.tel[tel_id])

        if str(geom) == 'FlashCam':  # FlashCam now skipped because of no proper pulse extraction
            continue
        elif str(geom) == 'ASTRI':   # excluded because no proper simulation in Prod3
            continue
                    
        # print(geom)

        # print(event.r0.tel[tel_id])
        data = event.r0.tel[tel_id].waveform
        # this is a 3D matrix num_gains * num_pixels * num_samples
    
        # print(data)
        # print(data.shape)
        
        # print(event.mc.tel[tel_id])
        ped = event.mc.tel[tel_id].pedestal
        # This one instead, is only numgains * num_pixels
        # the pedestal is the average (for pedestal events) of the *sum* of all samples, from sim_telarray

        # print(ped.shape)

        nsamples = data.shape[2]  # total number of samples
        pedcorrectedsamples = data - np.atleast_3d(ped)/nsamples    # Subtract pedestal baseline. atleast_3d converts 2D to 3D matrix

        integrator = LocalPeakIntegrator(None, None)
        integration, peakpos, window = integrator.extract_charge(pedcorrectedsamples) # these are 2D matrices num_gains * num_pixels
        # print(integration)
        # print(peakpos)
        # print(window)
        
        chan = 0  # high gain used for now...
        signals = integration[chan].astype(float)
        
        dc2pe = event.mc.tel[tel_id].dc_to_pe   # numgains * numpixels
        signals *= dc2pe[chan]
        # print(signals)

        # Add all individual pixel signals to the numpy array of the corresponding camera inside the log10pixelsignal dictionary
                # log10pixelHGsignal[str(geom)] = np.append(log10pixelHGsignal[str(geom)], np.log10(signals))
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
        # hillas = hillas_parameters_2(geom clean) # this one also gives some warnings invalid value in sqrt
        # print(hillas)


        # if ii == 0 :   # len(event.r0.tels_with_data)-1 :
        #     disp = CameraDisplay(geom)
        #     set_limits_minmax(0, np.amax(clean))
        #     disp.add_colorbar()
        #     # disp.image = signals
        #     disp.image = clean
        #     disp.highlight_pixels(cleanmask, color='black')
        #     plt.show()

        foclen = event.inst.subarray.tel[tel_id].optics.equivalent_focal_length

        w = np.rad2deg(np.arctan2(hillas.width,foclen));
        l = np.rad2deg(np.arctan2(hillas.length,foclen));

        camtype.append(str(geom))
        width = np.append(width, w.value)
        length = np.append(length, l.value)
        size = np.append(size, hillas.size)

