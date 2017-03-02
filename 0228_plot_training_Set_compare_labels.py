import os
import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
import pickle

import astropy.io.fits as ts
import AnniesLasso_2 as tc

import os
import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
import pickle

import astropy.io.fits as ts
import AnniesLasso_2 as tc


import math

def get_pixmask(flux, err):
    bad_flux = ~np.isfinite(flux)
    bad_err = (~np.isfinite(err)) | (err <= 0)
    bad_pixels = bad_err | bad_flux
    return bad_pixels




# train the cannon


## training the cannon

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()

training_set_path = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/reference_labels.csv"
training_set_spectrum_dir = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/Data/"

training_set = Table.read("reference_labels.csv")

N = len(training_set)
keep = np.ones(N, dtype=bool)

training_set_flux = []
training_set_ivar = []
training_set_error = []

for i, row in enumerate(training_set):
    image_path = os.path.join(training_set_spectrum_dir, row["ID"])
    if not os.path.exists(image_path):
        print("{}/{} could not be found: {}".format(i + 1, N, image_path))
        keep[i] = False
        continue
    print("{}/{}: {}".format(i + 1, N, image_path))
    image = fits.open(image_path)
    flux = image[1].data
    flux_err = image[2].data
    badpix = get_pixmask(flux, flux_err)
    ivar = 1.0 / flux_err ** 2
    error = flux_err
    # badpix is a array and the length is 8575
    flux[badpix] = 1.0
    ivar[badpix] = 0.0
    training_set_flux.append(flux)
    training_set_ivar.append(ivar)
    training_set_error.append(error)

training_set_flux = np.array(training_set_flux)
training_set_ivar = np.array(training_set_ivar)
training_set_error = np.array(training_set_error)






pkl_file = open('training_c20.pkl', 'rb')
c20 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('training_inferred_labels.pkl', 'rb')
inf = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('training_reference_labels.pkl', 'rb')
ref = pickle.load(pkl_file)
pkl_file.close()

# flux;
pkl_file = open('training_inf_flux.pkl', 'rb')
inf_flux = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('training_inf_flux_sim.pkl', 'rb')
inf_flux_sim = pickle.load(pkl_file)
pkl_file.close()

inf = c20[:,0:3]


##1
plt.suptitle("Inferred labels from fitting simultaneously ",fontsize =25)
plt.legend()

ax1 = plt.subplot(3,1,1)
ax1.plot(ref[:,0],inf[:,0],"go",alpha =0.3)
ax1.plot(ref[:,0],ref[:,0],"k",linewidth=0.5)

ax1.set_xlabel("Reference Teff",fontsize =15)
ax1.set_ylabel("Inferred Teff",fontsize =15)
pos1 = ax1.get_position()

pos2 = [pos1.x0, pos1.y0 ,  pos1.width, pos1.height/1.2]

ax1.set_position(pos2)



## 2

ax2 = plt.subplot(3,1,2)
ax2.plot(ref[:,1],inf[:,1],"bo",alpha =0.3)
ax2.plot(ref[:,1],ref[:,1],"k",linewidth=0.5)

ax2.set_xlabel("Reference Logg",fontsize =15)
ax2.set_ylabel("Inferred logg",fontsize =15)
pos1 = ax2.get_position()

pos2 = [pos1.x0, pos1.y0 ,  pos1.width, pos1.height/1.2]

ax2.set_position(pos2)

## 3

ax3 = plt.subplot(3,1,3)
ax3.plot(ref[:,2],inf[:,2],"go",alpha =0.3)
ax3.plot(ref[:,2],ref[:,2],"k",linewidth=0.5)

ax3.set_xlabel("Reference Fe/H",fontsize =15)
ax3.set_ylabel("Inferred Fe/H",fontsize =15)
pos1 = ax3.get_position()

pos2 = [pos1.x0, pos1.y0 ,pos1.width, pos1.height/1.2]

ax3.set_position(pos2)

plt.show()

## 4

"""
plt.subplot(2,2,4)

plt.suptitle("The chi-squared of the inferred spectrum")
print("calculating chi-squared -old")
plt.plot(np.sum(training_set_ivar*(inf_flux-training_set_flux)**2,axis=1),"ko",label="chi-squared of original inferred flux")

print("calculating chi-squared -new")
plt.plot(np.sum(training_set_ivar*(inf_flux_sim-training_set_flux)**2,axis=1),"bo",label="chi-squared of inferred flux from simultaneous fitting")

plt.xlabel("Star")
plt.ylabel("Chi-squared")


"""



