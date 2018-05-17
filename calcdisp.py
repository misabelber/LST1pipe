from astropy.io import fits
import matplotlib.pylab as plt
import numpy as np
import math as m
from ctapipe.instrument import CameraGeometry,OpticsDescription
from ctapipe.visualization import CameraDisplay
import astropy.units as u
import ctapipe.coordinates as c

hdu_list = fits.open("events.fits") #File with events
hdu_list[1].data

tel = OpticsDescription.from_name('LST') #Telescope description
focal_length = tel.equivalent_focal_length.value #Telescope focal length
geom = CameraGeometry.from_name("LSTCam") #Camera geometry

nevents = hdu_list[1].data.field(0).size #Total number of events
disp = np.array([]) #Disp quantity

for i in range(0,nevents):
    size = hdu_list[1].data.field(18)[i]
    if size < 150:
        continue

    #Calculate source position    
    mcAlt = hdu_list[1].data.field(4)[i] 
    mcAz = hdu_list[1].data.field(5)[i]
    mcAlttel = hdu_list[1].data.field(19)[i]
    mcAztel = hdu_list[1].data.field(20)[i]
    
    mcAlt = c.alt_to_theta(mcAlt*u.rad).value
    mcAz = c.az_to_phi(mcAz*u.rad).value
    mcAlttel = c.alt_to_theta(mcAlttel*u.rad).value
    mcAztel = c.az_to_phi(mcAztel*u.rad).value

    #Sines and cosines of direction angles
    cp = np.cos(mcAz)
    sp = np.sin(mcAz)
    ct = np.cos(mcAlt)
    st = np.sin(mcAlt)


    #Shower direction coordinates

    source = np.zeros(3)
    source[0] = st*cp
    source[1] = st*sp
    source[2] = ct
    
    

    #Rotation matrices towars the camera frame
    rot_Y = np.matrix([np.zeros(3),np.zeros(3),np.zeros(3)])
    rot_Z = np.matrix([np.zeros(3),np.zeros(3),np.zeros(3)])
 

    rot_Y[0,0] = np.cos(mcAlttel)
    rot_Y[0,1] = 0
    rot_Y[0,2] = np.sin(mcAlttel) 
    rot_Y[1,0] = 0
    rot_Y[1,1] = 1
    rot_Y[1,2] = 0
    rot_Y[2,0] = -np.sin(mcAlttel)
    rot_Y[2,1] = 0
    rot_Y[2,2] = np.cos(mcAlttel)
    
    rot_Z[0,0] = np.cos(mcAztel) 
    rot_Z[0,1] = -np.sin(mcAztel)
    rot_Z[0,2] = 0
    rot_Z[1,0] = np.sin(mcAztel)
    rot_Z[1,1] = np.cos(mcAztel)
    rot_Z[1,2] = 0
    rot_Z[2,0] = 0
    rot_Z[2,1] = 0
    rot_Z[2,2] = 1

    tmp = np.dot(rot_Y,source)
    tmp=tmp.getT()
    res = np.squeeze(np.asarray(np.dot(rot_Z,tmp)))

    Source_X = -focal_length*res[0]/res[2]
    Source_Y = -focal_length*res[1]/res[2]
    
    cen_x = hdu_list[1].data.field(16)[i]
    cen_y = hdu_list[1].data.field(17)[i]

    disp = np.append(disp,m.sqrt((Source_X-cen_x)**2+(Source_Y-cen_y)**2))


    display = CameraDisplay(geom)
    display.add_colorbar()
    
    image = hdu_list[2].data[i]

    display.image = image
    display.cmap = 'CMRmap'
    
    
    plt.plot([Source_X],[Source_Y],marker='o',markersize=10,color="green")
    plt.plot([cen_x],[cen_y],marker='x',markersize=10,color="blue")
    plt.plot([Source_X,cen_x],[Source_Y,cen_y],'-',color="red")
    print(size,mcAlttel,mcAztel,mcAlt,mcAz)
    plt.show()

plt.plot(disp,hdu_list[1].data.field(11)/hdu_list[1].data.field(12),'.')
plt.yscale('log')
plt.show()  
