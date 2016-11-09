from numpy import genfromtxt
import matplotlib.pyplot as plt
import random
import matplotlib
import numpy as np
import pickle
from astropy.table import Table
import TheCannon_2
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

# read data
all_set = Table.read("/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/allStar-v304.fits")

# data


## 1. read the dat
## 2. add a mask (bitmask)
## 3. filter-mean inverse variance






"""
# choose some of them
choose = []
a=0
for i in range(0,100):
    a += np.random.randint(1,11)
    choose.append(a)

choose = np.array(choose)
all_set = all_set[choose]

"""
#Let's rock!

N = len(all_set)

choose = []
a=0
for i in range(0,N):
    a += 1
    choose.append(a)

choose = np.array(choose)
all_set = all_set[choose]




N = len(all_set)
fail = 0
star=0

# choose some flux and normalize them first, then filter the mean inverse variance.


# get the time of the program
start_time = time.time()


for i, row in enumerate(all_set):
    try:


        image_path = "/Users/caojunzhi/Desktop/Data/APOGEE_DR10_Apstar/apStar-s3-" + row["APOGEE_ID"] + ".fits"

        if not os.path.exists(image_path):
            print("{}/{} could not be found: {}".format(i + 1, N, image_path))
            fail += 1
            continue

        print("{}/{}: {}".format(i + 1, N, image_path))
        # Let's only use two extension of the data

        image = fits.open(image_path, ignore_missing_end=True)
        dat = Table.read(image_path)

        flux = image[1].data
        flux_err = image[2].data


    except IOError:
        print("opts. This one fail")
        fail+=1

    else:
        badpix = get_pixmask(flux, flux_err)
        ivar = 1.0 / flux_err ** 2
        error = flux_err
        # badpix is a array and the length is 8575
        flux = np.array(flux, dtype=np.float64)
        ivar = np.array(ivar, dtype=np.float64)

        flux[badpix] = 1.0
        ivar[badpix] = 0.0

        flux = np.array(flux)
        ivar = np.array(ivar)

        try:

            test_labels_all_i = [row["TEFF"], row["LOGG"], row["METALS"]]
            # print(test_labels_all_i)


        except ValueError:
            print("opts fail")
            fail+=1

        else:
            star += 1
            # Fiber number
            m = dat[0]["FIBER"]
            FiberID = m

            tr_ID = image_path

            ## This is the real successful ones.

            ##let's optimize them one by one

            # normalization

            # trick to deal with 1D array:
            # by jason

            try:
                star = len(flux[:, 0])
                one = np.ones(len(flux[0, :]))

            except IndexError:
                star = 1
                one = np.ones(len(flux))

            else:
                flux = np.vstack((one, flux))
                ivar = np.vstack((one, ivar))

                ####
                # value

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

                # get inf_labels
                inf_labels = model.fit(norm_tr_flux, norm_tr_ivar)
                v = model.vectorizer.get_label_vector(inf_labels)
                inf_flux = np.dot(v, model.theta.T)
                opt_flux, parameters = model.fitting_spectrum_parameters_single \
                    (norm_tr_flux, norm_tr_ivar, inf_flux)


                # by jason

                inf_flux = inf_flux[1:star, :]
                opt_flux = opt_flux[1:star, :]
                parameters = parameters[1:star, :]

                FiberID = np.array(FiberID)

                print(inf_flux.shape, opt_flux.shape, parameters.shape, FiberID.shape)

            ## save them

            # save inf_flux

                output = open("/Volumes/data-2T/Data/opt_apogee_dr10/inf_flux/" + row["APOGEE_ID"] + ".pkl", 'wb')
                pickle.dump(inf_flux, output)
                output.close()

                output = open("/Volumes/data-2T/Data/opt_apogee_dr10/opt_flux/" + row["APOGEE_ID"] + ".pkl", 'wb')
                pickle.dump(opt_flux, output)
                output.close()

                output = open("/Volumes/data-2T/Data/opt_apogee_dr10/parameters/" + row["APOGEE_ID"] + ".pkl", 'wb')
                pickle.dump(parameters, output)
                output.close()
















# count for time



print("--- %s seconds ---" % (time.time() - start_time))










