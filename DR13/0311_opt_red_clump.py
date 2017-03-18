import astropy.io.fits as ts
from TheCannon_2 import dataset,apogee
from astropy.table import Table
import numpy as np
import os
from astropy.io import fits
import pickle
import AnniesLasso_2 as tc
import time
from os import listdir
from os.path import isfile, join

## training the cannon

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()



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
    flux[badpix] = np.median(flux)
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



# function
def get_pixmask(flux, err):
    bad_flux = ~np.isfinite(flux)
    bad_err = (~np.isfinite(err)) | (err <= 0)
    bad_pixels = bad_err | bad_flux
    return bad_pixels


# find mask


def find_mask(flux):

    # use 10 pixel:
    N_pixel = len(flux[0,:])
    N_star = len(flux[:,0])

    mask = []


    for i in range(0,N_star):

        print("Dealing with star %d"%(i+1))

        # one pixel is approximately 0.216\AA

        mask_i = np.zeros(N_pixel)

        width = 5
        limit = 0.95

        for j in range(10,N_pixel-10):
            flux_j = flux[i,j]
            min_j = np.min(flux[i,j-width:j+width])

            if flux_j <= min_j and flux_j < 1:
                mask_i[j-width:j+width] = 1
            else:
                nothing=1

        for j in range(0,N_pixel):
            flux_j = flux[i, j]
            if flux_j > limit:
                mask_i[j] = 0
            else:
                nothing=1

        mask.append(mask_i)
    mask = np.array(mask)
    print("mask shape")


    print(mask.shape)
    return mask


# theta:


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
theta_z = theta_z[0:row,:]


##########
#Run them on our stars and count time
##########


def fitting_ve(name):

    image_path = name
    if not os.path.exists(image_path):
        print("{}/{} could not be found: {}".format(i + 1, N, image_path))
        keep[i] = False

    # We only store flux,ivar,inf_flux,parameters,parameters_new,parameters_sim,ve(n*3)(include ve, ve_new,ve_sim)

    image = fits.open(image_path, ignore_missing_end=True)
    dat = Table.read(image_path)

    flux = image[1].data
    flux_err = image[2].data

    n_i = len(flux[:, 0])

    # You need to change a little bit if you need combined spectra

    """

    flux = flux[2:n_i, :]
    flux_err = flux_err[2:n_i, :]


    """

    badpix = get_pixmask(flux, flux_err)
    ivar = 1.0 / flux_err ** 2
    error = flux_err
    # badpix is a array and the length is 8575
    flux = np.array(flux, dtype=np.float64)
    ivar = np.array(ivar, dtype=np.float64)

    flux[badpix] = np.median(flux)
    ivar[badpix] = 0.0

    flux = np.array(flux)
    ivar = np.array(ivar)

    # normalize flux:
    # value

    tr_ID = image_path

    test_labels_all_i = np.array([5000, 1, 1])

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

    # infer labels


    #inf_labels = model.fit(norm_tr_flux, norm_tr_ivar)


    # Use inferred labels from the combined spectra:


    inf_labels = model.fit(norm_tr_flux, norm_tr_ivar)
    # only use the inf labels from the combined spectra

    com = len(inf_labels[:,0])

    inf_labels_com = inf_labels[0, :]

    inf_labels = []
    for z in range(0,com):
        inf_labels.append(inf_labels_com)


    inf_labels = np.array(inf_labels)





    v = model.vectorizer.get_label_vector(inf_labels)
    inf_flux = np.dot(v, model.theta.T)
    opt_flux, parameters = model.fitting_spectrum_parameters_single \
        (norm_tr_flux, norm_tr_ivar, inf_flux)

    ve = (parameters[:, 2] - parameters[:, 0]) / (parameters[:, 0] + parameters[:, 1] + parameters[:, 2]) * 4144.68

    ve_un = model.uncertainty

    # old
    a0 = parameters
    a1 = ve
    a2 = ve_un

    # spectra

    a3 = norm_tr_flux
    a4 = norm_tr_ivar
    a5 = inf_flux
    a6 = opt_flux

    a7 = inf_labels

    # save them

    # pay attention to the fits file saving

    path_fits_i = image_path.replace("apStar-r6-", "")


    print("saving files" + path_fits_i)

    hdu = fits.PrimaryHDU(data=a0)
    hdu.header[
        'COMMENT'] = "Inf labels are from the combined spectra"

    hdu.writeto(path_fits_i, clobber=True)

    ts.append(path_fits_i, a1)
    ts.append(path_fits_i, a2)
    ts.append(path_fits_i, a3)
    ts.append(path_fits_i, a4)
    ts.append(path_fits_i, a5)
    ts.append(path_fits_i, a6)
    ts.append(path_fits_i, a7)


# Let's rock!


start_time = time.time()

mypath = "/Users/caojunzhi/Desktop/Data/suspect_14"

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

N = len(onlyfiles)
onlyfiles = onlyfiles[1:N]
print(onlyfiles)

for i in range(0,len(onlyfiles)):

    name_i = mypath+"/"+onlyfiles[i]
    fitting_ve(name=name_i)

onlyfiles = np.array(onlyfiles)
output = open('suspect_14.pkl', 'wb')
pickle.dump(onlyfiles, output)
output.close()



stop_time = time.time()

print("The time we use %.2f s"%(stop_time-start_time))















