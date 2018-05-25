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


hdu_list = fits.open("events.fits") #File with events
hdu_list[1].data

tel = OpticsDescription.from_name('LST') #Telescope description
focal_length = tel.equivalent_focal_length.value #Telescope focal length
geom = CameraGeometry.from_name("LSTCam") #Camera geometry

nevents = hdu_list[1].data.field(0).size #Total number of events
disp = np.array([]) #Disp quantity

width = hdu_list[1].data.field(11)
length = hdu_list[1].data.field(12)
size = hdu_list[1].data.field(18)
phi = hdu_list[1].data.field(13)
energy = hdu_list[1].data.field(3)
cen_x = hdu_list[1].data.field(16)
cen_y = hdu_list[1].data.field(17)
psi = hdu_list[1].data.field(14)

X = np.array([width/length,size,phi,cen_x,cen_y,psi]).T
#X = np.array([size,width,length,phi]).T
X_train, X_test, y_train, y_test = train_test_split(X, energy,
                                                    train_size=int(9*nevents/10),
                                                    random_state=4)

max_depth = 50
regr_rf = RandomForestRegressor(max_depth=max_depth, random_state=2)                                                           
regr_rf.fit(X_train, y_train)
erec = regr_rf.predict(X_test)

h = (y_test-erec)
print(h.mean()*np.log(10),h.std()*h.std()*np.log(10))
plt.hist(h*np.log(10),bins=100,range=[-1,1])

plt.show()




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
