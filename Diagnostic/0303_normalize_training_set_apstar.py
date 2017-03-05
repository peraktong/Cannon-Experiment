import astropy.io.fits as ts
from TheCannon_2 import dataset,apogee
from astropy.table import Table
import numpy as np
import os
from astropy.io import fits
import pickle
import AnniesLasso_2 as tc
import time

## training the cannon

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()


training_set_path = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/reference_labels.csv"
training_set_spectrum_dir = "/Volumes/Data_2TB/Data/APOGEE_DR10_Apstar/"


training_set = Table.read("reference_labels.csv")



def get_pixmask(flux, err):
    bad_flux = ~np.isfinite(flux)
    bad_err = (~np.isfinite(err)) | (err <= 0)
    bad_pixels = bad_err | bad_flux
    return bad_pixels


N = len(training_set)
keep = np.ones(N, dtype=bool)



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
    flux = np.atleast_2d(image[1].data)[0,:]
    flux_err = np.atleast_2d(image[2].data)[0,:]

    badpix = get_pixmask(flux, flux_err)
    ivar = 1.0 / flux_err ** 2
    error = flux_err
    # badpix is a array and the length is 8575
    flux[badpix] = 1.0
    ivar[badpix] = 0.0

    ####
    # value

    test_labels_all_i = [row["Teff_{corr}"], row["logg_{corr}"], row["[M/H]_{corr}"]]
    tr_ID = image_path

    # dataset
    flux = np.atleast_2d(flux)
    ivar = np.atleast_2d(ivar)


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

    path_fits_i = "/Users/caojunzhi/Desktop/Data/training_set/"+name_i
    path_flux_i = image_path



    # save as separate fits
    print("doing %d"%(i+1))
    print("saving files " + path_fits_i)

    hdu = fits.PrimaryHDU(data=norm_tr_flux)
    hdu.header[
        'COMMENT'] = "normalized_flux and normalized ivar"

    hdu.writeto(path_fits_i, clobber=True)

    ts.append(path_fits_i, norm_tr_ivar)


