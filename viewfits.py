from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
hdu_list = fits.open("SPS_LSTCam.fits")
hdu_list.info()
hdu_list[1].data

plt.hist(hdu_list[1].data.field(0),bins=100,range=(0.,2.))
plt.yscale('log')
plt.show()  
