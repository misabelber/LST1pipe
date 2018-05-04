from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np

from ctapipe.io import event_source
from ctapipe.visualization import CameraDisplay
from ctapipe.instrument.camera import CameraGeometry
from ctapipe.calib import CameraCalibrator
from ctapipe.reco.HillasReconstructor import HillasReconstructor
from ctapipe.image.hillas import hillas_parameters
from ctapipe.image.cleaning import tailcuts_clean
from ctapipe.utils import linalg
from ctapipe.utils import datasets

filename = "/scratch/bernardos/LST1/gamma_20deg_180deg_run805___cta-prod3-merged_desert-2150m-Paranal-subarray-1-nosct.simtel.gz"
#filename = datasets.get_dataset_path("gamma_test_large.simtel.gz")
source = event_source(filename, allowed_tels={1, 2, 3, 4})

# for each event
off_angles = []      

reco = HillasReconstructor()
calib = CameraCalibrator()

for event in source:

    # shower direction
    # converting MC shower angular parameters to 3D spatial vector
    shower_azimuth = event.mc.az                        # same as in MC file
    shower_altitude = np.pi * u.rad / 2 - event.mc.alt  # altitude = 90 - mc.alt
    shower_direction = linalg.set_phi_theta(shower_azimuth, shower_altitude)

    # calibrating the event
    calib.calibrate(event)

    # for each telescope and event
    point_azimuth = {}
    point_altitude = {}
    hillas_params = {}

    # get hillas params for each event in different telescopes
    for tel_id in event.dl0.tels_with_data:

        # telescope pointing direction extracted from MC dataset
        point_azimuth[tel_id] = event.mc.tel[tel_id].azimuth_raw * u.rad
        point_altitude[tel_id] = (
            np.pi / 2 - event.mc.tel[tel_id].altitude_raw) * u.rad

        # camera geometry required for hillas parametrization
        camgeom = event.inst.subarray.tel[tel_id].camera

        # note that [0] is for channel 0 which is high-gain channel
        image = event.dl1.tel[tel_id].image[0]

        # cleaning  the image
        cleaned_image = image
        # create a clean mask of pixels above the threshold
        cleanmask = tailcuts_clean(
            camgeom, image, picture_thresh=10, boundary_thresh=5)
        # set all rejected pixels to zero
        cleaned_image[~cleanmask] = 0

        # calulate hillas parameters
        # it fails for empty pixels
        try:
            hillas_params[tel_id] = hillas_parameters(camgeom, cleaned_image)
        except:
            pass
    if len(hillas_params) < 2:
        continue
    # fit stereo event direction
    reco.get_great_circles(
        hillas_params, event.inst.subarray, point_azimuth, point_altitude)

    # reconstruct direction (3 components) with errors on the values
    reco_direction, reco_dir_err = reco.fit_origin_crosses()
    # in case fit reconstruct fails to get any real value -> next event
    if np.isnan(reco_direction).any():
        continue
    # get angular offset between reconstructed event direction and MC shower direction
    off_angle = linalg.angle(reco_direction, shower_direction)

    # appending all estimated off angles of each event
    off_angles.append(off_angle.to(u.deg).value)

# square the angles
thetasq = []
for i in off_angles:
    thetasq.append(i**2)

# plot 
plt.hist(thetasq, bins=np.linspace(0, 0.5, 300))
plt.xlabel(r'$\theta^2$ (deg)')
plt.ylabel("# of events")
plt.show()

