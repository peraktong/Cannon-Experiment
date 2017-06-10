import matplotlib.pyplot as plt

import pickle
import math
import numpy as np
import os
from astropy.io import fits
import pickle
from astropy.table import Table
import AnniesLasso_2 as tc

import plot_class_4_star_v2_1107

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()


def log10(x):
    return math.log10(x)

# sinc interpolation

# Distance between s is equal and they are your raw data
# u is what you want. Then length of u is not equal to s.

def sinc_interp(x, s, u):
    """
    Interpolates x, sampled at "s" instants
    Output y is sampled at "u" instants ("u" for "upsampled")
    """

    # Your x is the raw flux
    # Your s is the wave length of the raw flux. s should have equal distance

    # Your u is the log wl, which has equal distance between the neighborhoods.
    # The length of u is 8575 and you can use log wl

    if len(x) != len(s):
        print("len(x) should be equal to len(s")

    # Find the period
    # T_r = s[1] - s[0]

    # I don't think these two methods have a big different.

    # Can we use this?
    N = len(s)
    T = (s[N-1]-s[0])/N

    sincM = np.tile(u, (len(s), 1)) - np.tile(s[:, np.newaxis], (1, len(u)))
    y = np.dot(x, np.sinc(sincM / T))
    return y



# import data
image_path = "apVisit-r3-5094-55874-010.fits"
image = fits.open(image_path,ignore_missing_end=True)

dat = Table.read(image_path)


i=0

flux =image[1].data[i]

# HDU 4 is the wave length
print(flux.shape)
print(image[4].data.shape)
#print(wl[8574],wl[0])


"""
plt.plot(image[4].data[0],"ro")



axes = plt.gca()
# axes.set_xlim([15660,15780])
axes.set_xlim([1,30])
axes.set_ylim([16900,17000])
#axes.set_yticks(np.arange(0.8, 1.21, 0.1))

plt.show()
"""



wl_raw = image[4].data[0]


wl_raw = list(map(log10,wl_raw))

N = len(wl_raw)
print((wl_raw[N-1]-wl_raw[0])/(N-1),wl_raw[3]-wl_raw[2])

