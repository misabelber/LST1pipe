from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import math as m
hdu_list = fits.open("/scratch/bernardos/LST1/Events/events.fits")
hdu_list[1].data

nevents = hdu_list[1].data.field(0).size
disp = np.array([])

for i in range(0,nevents):
    #Calculate source position    
    mcAz = hdu_list[1].data.field(4)[i]
    mcAlt = hdu_list[1].data.field(5)[i]
    mcPhitel = hdu_list[1].data.field(19)[i]
    mcThetatel = hdu_list[1].data.field(20)[i]
    
    #Sines and cosines of direction angles
    cp = np.cos(mcAlt)
    sp = np.sin(mcAlt)
    ct = np.cos(mcAz)
    st = np.sin(mcAz)

    #Shower direction coordinates

    source = np.zeros(3)
    source[0] = st*cp
    source[1] = st*sp
    source[2] = ct
    
    

    #Rotation matrices towars the camera frame
    rot_Y = np.matrix([np.zeros(3),np.zeros(3),np.zeros(3)])
    rot_Z = np.matrix([np.zeros(3),np.zeros(3),np.zeros(3)])
 

    rot_Y[0,0] = np.cos(mcThetatel)
    rot_Y[0,1] = 0
    rot_Y[0,2] = np.sin(mcThetatel) 
    rot_Y[1,0] = 0
    rot_Y[1,1] = 1
    rot_Y[1,2] = 0
    rot_Y[2,0] = -np.sin(mcThetatel)
    rot_Y[2,1] = 0
    rot_Y[2,2] = np.cos(mcThetatel)
    
    rot_Z[0,0] = np.cos(mcPhitel) 
    rot_Z[0,1] = -np.sin(mcPhitel)
    rot_Z[0,2] = 0
    rot_Z[1,0] = np.sin(mcPhitel)
    rot_Z[1,1] = np.cos(mcPhitel)
    rot_Z[1,2] = 0
    rot_Z[2,0] = 0
    rot_Z[2,1] = 0
    rot_Z[2,2] = 1

    tmp = np.dot(rot_Y,source)
    tmp=tmp.getT()
    res = np.squeeze(np.asarray(np.dot(rot_Z,tmp)))

    Source_X = -28000*res[1]/res[2]
    Source_Y = 28000*res[0]/res[2]
    
    cen_x = hdu_list[1].data.field(16)[i]*1000
    cen_y = hdu_list[1].data.field(17)[i]*1000

    disp = np.append(disp,m.sqrt((Source_X-cen_x)**2+(Source_Y-cen_y)**2))

plt.plot(disp,hdu_list[1].data.field(11)/hdu_list[1].data.field(12),'.')
#plt.yscale('log')
plt.show()  
