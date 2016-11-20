import os
import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle
import matplotlib.pyplot as plt
import matplotlib

import AnniesLasso_2 as tc
import time

# Get training set data

# get the time of the program
start_time = time.time()


training_set_path = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/reference_labels.csv"
training_set_spectrum_dir = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/Data/"


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
    ivar = 1.0/flux_err**2
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
    training_set, training_set_flux, training_set_ivar,threads=8)
model.s2 = 0
model.regularization = 0
model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(training_set,
    tc.vectorizer.polynomial.terminator(("Teff_{corr}", "logg_{corr}", "[M/H]_{corr}"), 2))

#train it
model.train()

# infer training set flux:
inf_label_tr = model.fit(training_set_flux,training_set_ivar)

v_tr = model.vectorizer.get_label_vector(inf_label_tr)
# v_500 = model.vectorizer.get_label_vector(test_label)

inf_flux_tr = np.dot(v_tr,model.theta.T)

print(v_tr.shape,model.theta.shape)

## obtain parameters for the training set

opt_flux_tr , parameters_tr = model.fitting_spectrum_parameters_single(normalized_flux=training_set_flux,
                                                                       normalized_ivar=training_set_ivar,inf_flux=inf_flux_tr)


## save them:

output = open('inf_labels_tr.pkl', 'wb')
pickle.dump(inf_label_tr, output)
output.close()

output = open('inf_flux_tr.pkl', 'wb')
pickle.dump(inf_flux_tr, output)
output.close()

output = open('opt_flux_tr.pkl', 'wb')
pickle.dump(opt_flux_tr, output)
output.close()

output = open('parameters_tr.pkl', 'wb')
pickle.dump(parameters_tr, output)
output.close()



np.savetxt("parameters_tr.csv", parameters_tr, delimiter=",")


## Let's retrain with abc

print("Retrain the model with a b and c")


training_set = Table.read("reference_labels_abc.csv")

assert all(keep)

# Construct model.
model = tc.L1RegularizedCannonModel(
    training_set, training_set_flux, training_set_ivar,threads=8)
model.s2 = 0
model.regularization = 0
model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(training_set,
    tc.vectorizer.polynomial.terminator(("Teff_{corr}", "logg_{corr}", "[M/H]_{corr}","a","b","c"), 2))

#train it
model.train()

# infer training set flux:
inf_label_tr = model.fit(training_set_flux,training_set_ivar)

v_tr = model.vectorizer.get_label_vector(inf_label_tr)
# v_500 = model.vectorizer.get_label_vector(test_label)

inf_flux_tr = np.dot(v_tr,model.theta.T)

print(v_tr.shape,model.theta.shape)

print("inf_label check",inf_label_tr.shape)

# save inf_abc


output = open('inf_labels_tr_opt.pkl', 'wb')
pickle.dump(inf_label_tr[:,0:3], output)
output.close()

output = open('parameters_tr_inf.pkl', 'wb')
pickle.dump(inf_label_tr[:,3:6], output)
output.close()

print(inf_label_tr[:,0:3].shape,inf_label_tr[:,3:6].shape)



np.savetxt("parameters_tr_inf.csv", inf_label_tr[:,3:6], delimiter=",")

# You'd better choose another fluxes and check:

# load data
# load testing set

pkl_file = open('n_testing_labels_900.pkl', 'rb')
test_label = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n_testing_flux_900.pkl', 'rb')
test_flux = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n_testing_error_900.pkl', 'rb')
all_set_error = pickle.load(pkl_file)
pkl_file.close()



pkl_file = open('n_testing_ivar_900.pkl', 'rb')
test_ivar = pickle.load(pkl_file)
pkl_file.close()




inf_labels_opt = model.fit(test_flux,test_ivar)


# save:
output = open('n_testing_labels_900.pkl', 'wb')
pickle.dump(test_label, output)
output.close()


output = open('n_inf_labels_900_opt.pkl', 'wb')
pickle.dump(inf_labels_opt, output)
output.close()




