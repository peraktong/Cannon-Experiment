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

assert all(keep)

# Construct model.
model = tc.L1RegularizedCannonModel(
    training_set, training_set_flux, training_set_ivar, threads=8)
model.s2 = 0
model.regularization = 0
model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(training_set,
                                                                tc.vectorizer.polynomial.terminator((
                                                                    "Teff_{corr}",
                                                                    "logg_{corr}",
                                                                    "[M/H]_{corr}"),
                                                                    2))

model.train()




pkl_file = open('n_900_path_fits.pkl', 'rb')
path_fits = pickle.load(pkl_file)
pkl_file.close()


pkl_file = open('n_900_path_flux.pkl', 'rb')
path_flux = pickle.load(pkl_file)
pkl_file.close()


# mean_ivar

pkl_file = open('n_900_mean_ivar.pkl', 'rb')
mi = pickle.load(pkl_file)
pkl_file.close()

N = len(path_fits)

# obtain theta:

theta = model.theta

# Set the boundary to be 0
one = 0 * np.ones(len(theta[0, :]))

row = len(theta[:, 0])

# x
theta_x = np.vstack((theta, one))
theta_x = theta_x[1:row + 1, :]

# y
theta_y = theta

# z

theta_z = np.vstack((one, theta))
theta_z = theta_z[0, :row]

for i in range(0,N):

    star = fits.open(path_fits[i])

    flux = star[0].data
    ivar = star[1].data
    parameters = star[4].data[:,0:3]

    inf_labels = star[9].data

    #initial data
    initial_labels = np.hstack((inf_labels, parameters))

    n_i = len(flux[:,0])

    result = []

    inf_flux_sim = []


    for j in range(0,n_i):

        flux_j = np.atleast_2d(flux[j,:])
        ivar_j = np.atleast_2d(ivar[j,:])
        initial_j = np.atleast_2d(initial_labels[j,:])

        label_6,un = model.fit_opt(flux_j, ivar_j, initial_labels=initial_j)


        # new inferred flux:

        label_6 = np.ravel(label_6)


        a = label_6[3]
        b = label_6[4]
        c = label_6[5]

        inf_labels_sim = label_6[0:3]
        v_sim = model.vectorizer.get_label_vector(inf_labels_sim)
        inf_flux = a*np.dot(v_sim, theta_x.T)+b*np.dot(v_sim, theta_y.T)+c*np.dot(v_sim, theta_z.T)

        result.append(label_6)
        inf_flux_sim.append(inf_flux)



    result = np.array(result)
    inf_flux_sim = np.array(inf_flux_sim)

    print("saving %d"%(i+1))
    ts.append(path_fits[i], result)
    ts.append(path_fits[i], inf_flux_sim)



