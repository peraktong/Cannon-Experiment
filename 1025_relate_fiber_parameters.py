import os
import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle
import matplotlib.pyplot as plt
import matplotlib
from astropy.io.fits import getdata
from numpy.linalg import inv
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


assert all(keep)

# Construct model.
model = tc.L1RegularizedCannonModel(
    training_set, training_set_flux, training_set_ivar,threads=8)

"""
model.s2 = 0
model.regularization = 0
model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(training_set,
    tc.vectorizer.polynomial.terminator(("Teff_{corr}", "logg_{corr}", "[M/H]_{corr}"), 2))

model.train()

"""


# load parameters
pkl_file = open('parameters_500.pkl', 'rb')
parameters_500 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('FiberID.pkl', 'rb')
FiberID = pickle.load(pkl_file)
pkl_file.close()

# Let's rock

para_p, inf_fiber = model.fitting_fiber_number(parameters_500,FiberID)

print(para_p)

"""
# plot:
font = {'weight': 'bold','size': 30}
matplotlib.rc('font', **font)
fig = plt.figure()

plt.plot(FiberID,inf_fiber, "bo", label="a", linewidth=3.0)
plt.plot(parameters_500[:,1],parameters_500[:,1], "k", linewidth=3.0)


plt.legend(loc="best")

fig.suptitle('Fiber vs inf_fiber', fontsize=40)
plt.xlabel('fiber_number', fontsize=38)
plt.ylabel('inf_fiber_number', fontsize=36)
plt.show()

"""
