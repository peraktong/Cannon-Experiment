import TheCannon_2
from TheCannon_2 import dataset,apogee
from astropy.table import Table
import numpy as np
import os
from astropy.io import fits
import pickle
import matplotlib.pyplot as plt
import matplotlib
import AnniesLasso_2 as tc


# function
def get_pixmask(flux, err):
    bad_flux = ~np.isfinite(flux)
    bad_err = (~np.isfinite(err)) | (err <= 0)
    bad_pixels = bad_err | bad_flux
    return bad_pixels

def get_error_mask(err):
    bad_err = (~np.isfinite(err))
    bad_pixels = bad_err
    return bad_pixels

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()

# read data

# open the choose file
pkl_file = open('n_tr_ID_900.pkl', 'rb')
tr_ID_900 = pickle.load(pkl_file)
pkl_file.close()




# train the model:

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


# read fluxes for all visits

n_star = len(tr_ID_900)
inf_labels_deviation_teff = []
inf_labels_deviation_logg = []
inf_labels_deviation_fe = []

for i in range(0,n_star):
    image_path = tr_ID_900[i]
    image = fits.open(image_path, ignore_missing_end=True)
    dat = Table.read(image_path)
    flux = image[1].data
    flux_err = image[2].data
    print("Working on %.2f and %d"%(i/n_star,i))

    #put a mask and normalize

    badpix = get_pixmask(flux, flux_err)
    ivar = 1.0 / flux_err ** 2
    error = flux_err
    # badpix is a array and the length is 8575
    flux = np.array(flux, dtype=np.float64)
    ivar = np.array(ivar, dtype=np.float64)

    flux[badpix] = 1.0
    ivar[badpix] = 0.0

    # value
    tr_flux = flux
    tr_ivar = ivar

    test_labels_all_i = np.array((dat[0]["TEFF"], dat[0]["LOGG"], dat[0]["FEH"]))

    tr_ID = tr_ID_900[i]
    test_flux = tr_flux
    test_ivar = tr_ivar

    # normalize them:

    ds = dataset.Dataset(wl, tr_ID, flux, ivar,
                         test_labels_all_i, tr_ID, flux, ivar)

    ds.ranges = [[371, 3192], [3697, 5997], [6461, 8255]]

    # set sudo-continuous spectrum
    pseudo_tr_flux, pseudo_tr_ivar = ds.continuum_normalize_training_q \
        (q=0.90, delta_lambda=50)

    # set mask
    contmask = ds.make_contmask(pseudo_tr_flux, pseudo_tr_ivar, frac=0.07)

    # get continuous mask

    ds.set_continuum(contmask)

    # fit the normalized-spectrum in the continuous region

    cont = ds.fit_continuum(3, "sinusoid")

    # Obtain the normalized flux
    norm_tr_flux, norm_tr_ivar, norm_test_flux, norm_test_ivar = \
        ds.continuum_normalize(cont)

    # get inf_labels for teff
    inf_labels_i = model.fit(norm_tr_flux, norm_tr_ivar)
    teff = inf_labels_i[:,0]
    logg = inf_labels_i[:,1]
    fe = inf_labels_i[:,2]


    one = np.ones(len(teff))

    inf_labels_deviation_teff.extend(teff-one*teff[0])
    print(len(inf_labels_deviation_teff))

    one = np.ones(len(logg))
    inf_labels_deviation_logg.extend(logg-one*logg[0])

    one = np.ones(len(fe))
    inf_labels_deviation_fe.extend(fe-one*fe[0])


inf_labels_deviation_teff = np.array(inf_labels_deviation_teff)
inf_labels_deviation_logg = np.array(inf_labels_deviation_logg)
inf_labels_deviation_fe = np.array(inf_labels_deviation_fe)

# save them and check

print(inf_labels_deviation_teff.shape,inf_labels_deviation_logg.shape,inf_labels_deviation_fe.shape)

output = open('inf_labels_deviation_teff_900.pkl', 'wb')
pickle.dump(inf_labels_deviation_teff, output)
output.close()

output = open('inf_labels_deviation_logg_900.pkl', 'wb')
pickle.dump(inf_labels_deviation_logg, output)
output.close()

output = open('inf_labels_deviation_fe_900.pkl', 'wb')
pickle.dump(inf_labels_deviation_fe, output)
output.close()




# plot histogram

font = {'weight': 'bold','size': 30}
matplotlib.rc('font', **font)
fig = plt.figure()



colors = ['b', 'g', 'r']
name = ["Delta Teff","Delta logg","Delta Fe/H"]

# Plot histogram
fig = plt.figure()
plt.hist(inf_labels_deviation_teff,bins=15,color=colors[0],label=name[0])
plt.legend(prop={'size': 30})
fig.suptitle('Histogram for the deviations of the inferred Teff',fontsize =30)
plt.xlabel('values of Delta Teff', fontsize=30)
plt.ylabel('Number of stars', fontsize=30)
plt.show()


colors = ['b', 'g', 'r']
name = ["Delta Teff","Delta logg","Delta Fe/H"]

# Plot histogram
fig = plt.figure()
plt.hist(inf_labels_deviation_logg,bins=15,color=colors[1],label=name[1])
plt.legend(prop={'size': 30})
fig.suptitle('Histogram for the deviations of the inferred logg',fontsize =30)
plt.xlabel('values of Delta logg', fontsize=30)
plt.ylabel('Number of stars', fontsize=30)
plt.show()


colors = ['b', 'g', 'r']
name = ["Delta Teff","Delta logg","Delta Fe/H"]

# Plot histogram
fig = plt.figure()
plt.hist(inf_labels_deviation_fe,bins=15,color=colors[2],label=name[2])
plt.legend(prop={'size': 30})
fig.suptitle('Histogram for the deviations of the inferred Fe/H',fontsize =30)
plt.xlabel('values of Delta Fe/H', fontsize=30)
plt.ylabel('Number of stars', fontsize=30)
plt.show()

