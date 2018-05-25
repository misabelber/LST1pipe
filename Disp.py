import numpy as np
import ctapipe.coordinates as c
import astropy.units as u

def calc_CamSourcePos(mcAlt,mcAz,mcAlttel,mcAztel,focal_length):
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

    sourcex = st*cp
    sourcey = st*sp
    sourcez = ct

    source = np.array([sourcex,sourcey,sourcez])
    
    

    #Rotation matrices towars the camera frame
    rot_Y = np.matrix([np.zeros(3),np.zeros(3),np.zeros(3)])
    rot_Z = np.matrix([np.zeros(3),np.zeros(3),np.zeros(3)])
    
    rot_Y = np.array([[np.cos(mcAlttel),0,np.sin(mcAlttel)],
                      [0,1,0],
                      [-np.sin(mcAlttel),0,np.cos(mcAlttel)]])

    rot_Z = np.array([[np.cos(mcAztel),-np.sin(mcAztel),0],
                      [np.sin(mcAztel),np.cos(mcAztel),0],
                      [0,0,1]])

    '''
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
    '''

    tmp = np.dot(rot_Y,source)
    tmp=tmp.getT()
    res = np.squeeze(np.asarray(np.dot(rot_Z,tmp)))

    Source_X = -focal_length*res[0]/res[2]
    Source_Y = -focal_length*res[1]/res[2]
    return Source_X, Source_Y

def calc_DISP(Source_X,Source_Y,cen_x,cen_y):
    disp = np.sqrt((Source_X-cen_x)**2+(Source_Y-cen_y)**2)
    return disp
