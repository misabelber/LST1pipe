from astropy.io import fits 
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.multioutput import MultiOutputRegressor
import math as m
from ctapipe.instrument import CameraGeometry,OpticsDescription
from ctapipe.visualization import CameraDisplay
import astropy.units as u
import ctapipe.coordinates as c
import matplotlib as mpl
import Disp


hdu_list = fits.open("/home/queenmab/DATA/LST1/Events/events.fits") #File with events
hdu_list[1].data

tel = OpticsDescription.from_name('LST') #Telescope description
focal_length = tel.equivalent_focal_length.value #Telescope focal length
geom = CameraGeometry.from_name("LSTCam") #Camera geometry

nevents = hdu_list[1].data.field(0).size #Total number of events
disp = np.array([]) #Disp quantity

data = hdu_list[1].data

width = data.field('width')
length = data.field('length')
size = data.field('size')
phi = data.field('phi')
energy = np.log10(data.field('mcEnergy')*1e3) #Log of energy in GeV
cen_x = data.field('cen_x')
cen_y = data.field('cen_y')
psi = data.field('psi')
dist = np.sqrt(cen_x*cen_x + cen_y*cen_y) 

mcAlt = data.field('mcAlt')
mcAz = data.field('mcAz')
mcAlttel = data.field('mcAlttel')
mcAztel = data.field('mcAztel')

sourcepos = Disp.calc_CamSourcePos(mcAlt,mcAz,mcAlttel,mcAztel,focal_length)
disp = Disp.calc_DISP(sourcepos[0],sourcepos[1],cen_x,cen_y)

X_e = np.array([np.log10(size),dist,width,length,phi]).T
#X = np.array([size,width,length,phi]).T
X_etrain, X_etest, E_train, E_test = train_test_split(X_e, energy,
                                                    train_size=int(2*nevents/3),
                                                    random_state=4)

max_depth = 50
regr_rf = RandomForestRegressor(max_depth=max_depth, random_state=2)                                                           
regr_rf.fit(X_etrain, E_train)
erec = regr_rf.predict(X_etest)

difE = ((E_test-erec))
print(difE.mean(),difE.std()*difE.std())

figE, ax = plt.subplots()
hE = ax.hist2d(E_test,erec,bins=50)
plt.colorbar(hE[3],ax=ax)
figE.show()

X_disp = np.array([width/length,size,phi]).T
X_dtrain, X_dtest, Disp_train, Disp_test = train_test_split(X_disp, disp,
                                                    train_size=int(2*nevents/3),
                                                    random_state=4)

regr_rf.fit(X_dtrain, Disp_train)
Disprec = regr_rf.predict(X_dtest)

difD = ((Disp_test-Disprec))
print(difD.mean(),difD.std()*difD.std())

figD, aax = plt.subplots()
hD = aax.hist2d(Disp_test,Disprec,bins=50)
plt.colorbar(hD[3],ax=aax)
figD.show()



'''
rng = np.random.RandomState(1)
X = np.array([(np.sort(200 * rng.rand(600, 1) - 100,axis=0)).ravel(),(np.sort(200 * rng.rand(600, 1) - 100,axis=0)).ravel()])
y = np.array(2*X[0]+5*X[1])
y += (0.5 - rng.rand(*y.shape))

X_train, X_test, y_train, y_test = train_test_split(X.T, y,
                                                    train_size=400,
                                                    random_state=4)

max_depth = 30
#regr_multirf = MultiOutputRegressor(RandomForestRegressor(max_depth=max_depth,
#                                                          random_state=0))
#regr_multirf.fit(X_train, y_train)
regr_rf = RandomForestRegressor(max_depth=max_depth, random_state=2)
regr_rf.fit(X_train, y_train)

#y_multirf = regr_multirf.predict(X_test)
y_rf = regr_rf.predict(X_test)
'''
