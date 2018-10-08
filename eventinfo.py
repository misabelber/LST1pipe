#!/usr/bin/env python3


import matplotlib.pylab as plt
import astropy.units as u
import numpy as np
import ctapipe
import os
import copy
import sys
from ctapipe.core import Container, Field, Map
from ctapipe.instrument import CameraGeometry
from ctapipe.visualization import CameraDisplay
from ctapipe.image import hillas_parameters, hillas_parameters_2, tailcuts_clean
from ctapipe.io.eventsourcefactory import EventSourceFactory
from ctapipe.image.charge_extractors import LocalPeakIntegrator
from astropy.visualization import quantity_support
from astropy.table import Table
from astropy.io import fits

class EventContainer(Container):
    event = Field(Map(),"Event")

DATA_PATH="/scratch/bernardos/LST1/Gamma/"

source = EventSourceFactory.produce(input_url=DATA_PATH+"gamma_20deg_180deg_run2200___cta-prod3-demo-2147m-LaPalma-demo3-FC_cone8.simtel.gz",allowed_tels={1});

ev = EventContainer()

n=0;
for event in source:
    ev.event[n] = ctapipe.io.containers.DataContainer()
    ev.event[n] = copy.deepcopy(event)
    n=n+1
