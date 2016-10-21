import os
import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle
import matplotlib.pyplot as plt
import matplotlib
from astropy.io.fits import getdata

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
model.s2 = 0
model.regularization = 0
model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(training_set,
    tc.vectorizer.polynomial.terminator(("Teff_{corr}", "logg_{corr}", "[M/H]_{corr}"), 2))

model.train()


# load testing set

pkl_file = open('testing_labels_500.pkl', 'rb')
test_label = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('testing_flux_500.pkl', 'rb')
test_flux = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('testing_error_500.pkl', 'rb')
all_set_error = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('nor_testing_flux_500.pkl', 'rb')
nor_test_flux = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('testing_ivar_500.pkl', 'rb')
test_ivar = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('nor_testing_ivar_500.pkl', 'rb')
nor_test_ivar = pickle.load(pkl_file)
pkl_file.close()


pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('FiberID.pkl', 'rb')
FiberID = pickle.load(pkl_file)
pkl_file.close()




# optimize
# Try to use inferred labels instead of the true labels

inf_label_500 = model.fit(nor_test_flux,nor_test_ivar)

v_500 = model.vectorizer.get_label_vector(inf_label_500)
# v_500 = model.vectorizer.get_label_vector(test_label)

inf_flux_500 = np.dot(v_500,model.theta.T)

opt_flux_500,parameters_500 = model.fitting_spectrum_parameters_single\
    (test_flux,test_flux,inf_flux_500)



# save

output = open('inf_flux_500.pkl', 'wb')
pickle.dump(inf_flux_500, output)
output.close()

output = open('opt_flux_500.pkl', 'wb')
pickle.dump(opt_flux_500, output)
output.close()

output = open('inf_label_500.pkl', 'wb')
pickle.dump(inf_label_500, output)
output.close()





output = open('parameters_500.pkl', 'wb')
pickle.dump(parameters_500, output)
output.close()


# do some plot

#Let's plot them

p = np.random.randint(0,len(nor_test_flux[:,0])-1)

font = {'weight': 'bold','size': 13}
matplotlib.rc('font', **font)
fig = plt.figure()


plt.plot(wl,nor_test_flux[p,:],"k",label="data")
plt.plot(wl,inf_flux_500[p,:],"r",label = "inf")
plt.plot(wl,opt_flux_500[p,:],"g",label="opt")


axes = plt.gca()
#axes.set_xlim([15660,15780])
axes.set_xlim([16160,16280])
axes.set_ylim([0.8,1.21])
axes.set_yticks(np.arange(0.8,1.21,0.1))

plt.show()

