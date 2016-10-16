
import os
import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle

import AnniesLasso_2 as tc


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

# save training_set_flux, ivar and error
"""
output = open('training_set_flux.pkl', 'wb')
pickle.dump(training_set_flux, output)
output.close()

output = open('training_set_ivar.pkl', 'wb')
pickle.dump(training_set_ivar, output)
output.close()

output = open('training_set_error.pkl', 'wb')
pickle.dump(training_set_error, output)
output.close()

"""

assert all(keep)

# Construct model.
model = tc.L1RegularizedCannonModel(
    training_set, training_set_flux, training_set_ivar,threads=8)

model.s2 = 0
model.regularization = 0
model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(training_set,
    tc.vectorizer.polynomial.terminator(("Teff_{corr}", "logg_{corr}", "[M/H]_{corr}"), 2))



model.train()

# save model
# model.save("model")

# compare_v1


inferred_label =model.fit_labelled_set()
inferred_label_opt = model.fit_labelled_set_opt()

output = open('inferred_labels.pkl', 'wb')
pickle.dump(inferred_label, output)
output.close()

output = open('inferred_labels_opt.pkl', 'wb')
pickle.dump(inferred_label_opt, output)
output.close()



# optimized spectrum
opt_flux,theta_opt = model.fitting_spectrum_parameters(training_set_flux,training_set_ivar)

inferred_labels = model.fit_labelled_set()
v = model.vectorizer.get_label_vector(inferred_labels)

inf_flux = np.dot(v,model.theta.T)

output = open('inf_flux.pkl', 'wb')
pickle.dump(inf_flux, output)
output.close()

output = open('opt_flux.pkl', 'wb')
pickle.dump(opt_flux, output)
output.close()




# compare_v2 Re-train model

"""
inf_opt,theta_opt = model.fitting_spectrum_parameters(normalized_flux=training_set_flux
                                                      ,normalized_ivar=training_set_ivar)
model.theta = theta_opt
#re-train
print("retrain model")
model.train()
inferred_label_opt_v2 = model.fit_labelled_set()

output = open('inferred_labels_opt_v2.pkl', 'wb')
pickle.dump(inferred_label_opt_v2, output)
output.close()

"""




