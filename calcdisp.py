import sys
from astropy.io import fits
import matplotlib.pylab as plt
import numpy as np
import math as m
from ctapipe.instrument import CameraGeometry,OpticsDescription
from ctapipe.visualization import CameraDisplay
import astropy.units as u
import ctapipe.coordinates as c
from scipy.optimize import minimize, newton
import Disp

def function(pars):
    return sum((disp - (pars[0]+pars[1]*(width/length)+pars[2]*size))**2)

hdu_list = fits.open("/scratch/bernardos/LST1/Events/events.fits") #File with events
hdu_list[1].data

tel = OpticsDescription.from_name('LST') #Telescope description
focal_length = tel.equivalent_focal_length.value #Telescope focal length
geom = CameraGeometry.from_name("LSTCam") #Camera geometry

nevents = hdu_list[1].data.field(0).size #Total number of events
disp = np.array([]) #Disp quantity

width = np.array([])
length = np.array([])
size = np.array([])
energy = np.array([])
ntrain=0;
for i in range(0,nevents):
    if i%2==0:
        continue
    ntrain=ntrain+1
    this_size = hdu_list[1].data.field(18)[i]
    if this_size < 180:
        continue
    width = np.append(width,hdu_list[1].data.field(11)[i])
    length = np.append(length,hdu_list[1].data.field(12)[i])
    size = np.append(size,hdu_list[1].data.field(18)[i])
    energy = np.append(energy,hdu_list[1].data.field(3)[i])

    #Calculate source position    
    mcAlt = hdu_list[1].data.field(4)[i] 
    mcAz = hdu_list[1].data.field(5)[i]
    mcAlttel = hdu_list[1].data.field(19)[i]
    mcAztel = hdu_list[1].data.field(20)[i]

    srcpos = Disp.calc_CamSourcePos([mcAlt],[mcAz],[mcAlttel],[mcAztel],focal_length)
    
    Source_X = srcpos[0]
    Source_Y = srcpos[1]
    
    cen_x = hdu_list[1].data.field(16)[i]
    cen_y = hdu_list[1].data.field(17)[i]

    disp = Disp.calc_DISP(Source_X,Source_Y,cen_x,cen_y)
        
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
    
#plt.plot(width/length,disp,'.')
#plt.yscale('log')
#plt.show()  



print(ntrain)
 
pars = [0,1,1]

result = minimize(function,pars,method='nelder-mead')
print(function(result.x))

print(result.x[0],result.x[1],result.x[2])

disp_ = np.array([]) #Disp quantity
estimate = np.array([])

width_ = np.array([])
length_ = np.array([])
size_ = np.array([])
energy_ = np.array([])
ntest=0;


for i in range(0,nevents):
    if i%2!=0:
        continue
    ntest=ntest+1
    this_size_ = hdu_list[1].data.field(18)[i]
    if this_size_ < 180:
        continue
    width_ = np.append(width,hdu_list[1].data.field(11)[i])
    length_ = np.append(length,hdu_list[1].data.field(12)[i])
    size_ = np.append(size,hdu_list[1].data.field(18)[i])
    energy_ = np.append(energy,hdu_list[1].data.field(3)[i])

    #Calculate source position    
    mcAlt = hdu_list[1].data.field(4)[i] 
    mcAz = hdu_list[1].data.field(5)[i]
    mcAlttel = hdu_list[1].data.field(19)[i]
    mcAztel = hdu_list[1].data.field(20)[i]

    srcpos = Disp.calc_CamSourcePos([mcAlt],[mcAz],[mcAlttel],[mcAztel],focal_length)
    
    Source_X = srcpos[0]
    Source_Y = srcpos[1]

    cen_x = hdu_list[1].data.field(16)[i]
    cen_y = hdu_list[1].data.field(17)[i]

    disp_ = Disp.calc_DISP(Source_X,Source_Y,cen_x,cen_y)
    
    
estimate = result.x[0]+result.x[1]*width_/length_+result.x[2]*size_



    
