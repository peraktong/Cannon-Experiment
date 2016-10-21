import TheCannon_2
from TheCannon_2 import dataset,apogee
from astropy.table import Table
import numpy as np
import os
from astropy.io import fits
import pickle
import matplotlib.pyplot as plt
import matplotlib


# function
def get_pixmask(flux, err):
    bad_flux = ~np.isfinite(flux)
    bad_err = (~np.isfinite(err)) | (err <= 0)
    bad_pixels = bad_err | bad_flux
    return bad_pixels

# read data
all_set = Table.read("/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/allStar-v304.fits")
all_set_spectrum_dir = "/Volumes/data-2T/Data/APOGEE_DR10_ASPCAP/"

# data
all_set_flux = []
all_set_ivar = []
all_set_error = []
test_labels_all = []
fail = 0

# choose some of them
choose = []
a=0
for i in range(0,800):
    a += np.random.randint(1,10)
    choose.append(a)

choose = np.array(choose)
all_set = all_set[choose]


# save the chosen star
output = open('choose_500.pkl', 'wb')
pickle.dump(choose, output)
output.close()


N = len(all_set)
keep = np.ones(N, dtype=bool)
tr_ID = []
FiberID = []


for i, row in enumerate(all_set):
    try:
        image_path = "/Users/caojunzhi/Desktop/Data/APOGEE_DR10_ASPCAP/aspcapStar-" + row["APOGEE_ID"] + ".fits"
        # For apstar

        image_path_2 = "/Users/caojunzhi/Desktop/Data/APOGEE_DR10_Apstar/apStar-s3-" + row["APOGEE_ID"] + ".fits"

        if not os.path.exists(image_path):
            print("{}/{} could not be found: {}".format(i + 1, N, image_path))
            keep[i] = False
            continue

        print("{}/{}: {}".format(i + 1, N, image_path))
        # Let's only use two extension of the data

        image = fits.open(image_path, ignore_missing_end=True)
        dat = Table.read(image_path_2)

        flux = image[1].data
        flux_err = image[2].data

    except IOError:
        print("opts. This one fail")
        fail +=1
        keep[i] = False
    else:
        badpix = get_pixmask(flux, flux_err)
        ivar = 1.0 / flux_err ** 2
        error = flux_err
        # badpix is a array and the length is 8575
        flux[badpix] = 1.0
        ivar[badpix] = 0.0

        all_set_flux.append(flux)
        all_set_ivar.append(ivar)
        all_set_error.append(error)
        test_labels_all_i = [row["RV_TEFF"],row["RV_LOGG"],row["RV_FEH"]]
        test_labels_all.append(test_labels_all_i)


        # Fiber number
        m = dat[0]["FIBER"]
        m =np.array(m)

        print(m)
        print(type(m))

        try:
            FiberID.append(sum(m) / len(m))
        except TypeError:
            FiberID.append(m)



        tr_ID.append(image_path)


all_set_flux = np.array(all_set_flux)
all_set_ivar = np.array(all_set_ivar)
all_set_error = np.array(all_set_error)
test_labels_all =np.array(test_labels_all)
tr_ID = np.array(tr_ID)
FiberID = np.array(FiberID)

print(all_set_flux.shape,
      all_set_ivar.shape,all_set_error.shape,
      test_labels_all.shape,tr_ID.shape,FiberID.shape)


# now we have the data, let's normalize them

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()

# value
tr_flux = all_set_flux
tr_ivar = all_set_ivar
tr_label = test_labels_all

test_ID = tr_ID
test_flux =tr_flux
test_ivar =tr_ivar

ds = dataset.Dataset(wl, tr_ID, tr_flux, tr_ivar,
                     tr_label, test_ID, test_flux, test_ivar)

ds.ranges = [[371,3192], [3697,5997], [6461,8255]]

# set sudo-continuous spectrum
pseudo_tr_flux, pseudo_tr_ivar = ds.continuum_normalize_training_q\
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

#Let's plot them

p = np.random.randint(0,len(norm_tr_flux[:,0])-1)

font = {'weight': 'bold','size': 13}
matplotlib.rc('font', **font)
fig = plt.figure()


plt.plot(wl,tr_flux[p,:],"k",label="data")
plt.plot(wl,norm_tr_flux[p,:],"r",label = "nor")

axes = plt.gca()
#axes.set_xlim([15660,15780])
axes.set_xlim([16160,16280])
axes.set_ylim([0.8,1.21])
axes.set_yticks(np.arange(0.8,1.21,0.1))

plt.legend(loc="best")
plt.xlabel('reference labels', fontsize=18)
plt.ylabel('inferred labels', fontsize=18)

plt.show()

# save the testing set


output = open('testing_labels_500.pkl', 'wb')
pickle.dump(tr_label, output)
output.close()

output = open('tr_ID_500.pkl', 'wb')
pickle.dump(tr_ID, output)
output.close()


output = open('testing_flux_500.pkl', 'wb')
pickle.dump(tr_flux, output)
output.close()

output = open('testing_error_500.pkl', 'wb')
pickle.dump(all_set_error, output)
output.close()


output = open('nor_testing_flux_500.pkl', 'wb')
pickle.dump(norm_tr_flux, output)
output.close()


output = open('testing_ivar_500.pkl', 'wb')
pickle.dump(tr_ivar, output)
output.close()

output = open('nor_testing_ivar_500.pkl', 'wb')
pickle.dump(norm_tr_ivar, output)
output.close()

output = open('FiberID.pkl', 'wb')
pickle.dump(FiberID, output)
output.close()



