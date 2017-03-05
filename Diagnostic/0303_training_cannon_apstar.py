
import astropy.io.fits as ts
from TheCannon_2 import dataset,apogee
from astropy.table import Table
import numpy as np
import os
from astropy.io import fits
import pickle
import AnniesLasso_2 as tc
import time
import matplotlib.pyplot as plt

## training the cannon

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()



training_set_path = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/reference_labels.csv"
training_set_spectrum_dir = "/Users/caojunzhi/Desktop/Data/training_set/"


training_set = Table.read("reference_labels.csv")



def get_pixmask(flux, err):
    bad_flux = ~np.isfinite(flux)
    bad_err = (~np.isfinite(err)) | (err <= 0)
    bad_pixels = bad_err | bad_flux
    return bad_pixels


N = len(training_set)
keep = np.ones(N, dtype=bool)

training_set_flux = []
training_set_ivar = []
training_set_error = []
ref_labels = []

for i, row in enumerate(training_set):

    name_i = row["ID"]
    name_i = name_i.replace("aspcapStar-v304-", "")
    name_i = "apStar-s3-"+name_i

    image_path = os.path.join(training_set_spectrum_dir, name_i)
    if not os.path.exists(image_path):
        print("{}/{} could not be found: {}".format(i + 1, N, image_path))
        keep[i] = False
        continue
    print("{}/{}: {}".format(i + 1, N, image_path))
    image = fits.open(image_path)
    flux = image[0].data
    flux_err = (image[1].data)**(-0.5)

    badpix = get_pixmask(flux, flux_err)
    ivar = 1.0/flux_err**2
    error = flux_err

    # badpix is a array and the length is 8575
    flux[badpix] = 1.0
    ivar[badpix] = 0.0
    #labels
    ref_labels.append([row["Teff_{corr}"], row["logg_{corr}"], row["[M/H]_{corr}"]])

    flux = flux.ravel()
    ivar = ivar.ravel()

    training_set_flux.append(flux)
    training_set_ivar.append(ivar)


training_set_flux = np.array(training_set_flux)
training_set_ivar = np.array(training_set_ivar)
ref_labels = np.array(ref_labels)
training_set = training_set[keep]

print(training_set_flux.shape,training_set_ivar.shape,np.array(training_set).shape,ref_labels.shape)
#assert all(keep)

# Construct model.
model = tc.L1RegularizedCannonModel(
    training_set, training_set_flux, training_set_ivar,threads=8)
model.s2 = 0
model.regularization = 0
model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(training_set,
    tc.vectorizer.polynomial.terminator(("Teff_{corr}", "logg_{corr}", "[M/H]_{corr}"), 2))

model.train()

### test:
inf_labels =model.fit(training_set_flux,training_set_ivar)
v = model.vectorizer.get_label_vector(inf_labels)
inf_flux = np.dot(v, model.theta.T)

chi = (training_set_flux-inf_flux)**2*training_set_ivar

chi = np.sum(chi,axis=1)

print("chi = %.2f"%(np.mean(chi)))

"""
plt.plot(ref_labels[:,0],inf_labels[:,0],"ro")
plt.plot(ref_labels[:,0],ref_labels[:,0],"k")
plt.show()

"""

