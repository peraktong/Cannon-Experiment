import TheCannon_2
from TheCannon_2 import dataset,apogee
from astropy.table import Table
import numpy as np
import os
from astropy.io import fits
import pickle
import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
from astropy.table import Table
from astropy.io import fits

import pickle
import matplotlib.pyplot as plt
import matplotlib
from astropy.io.fits import getdata

import AnniesLasso_2 as tc



# function
def get_pixmask(flux, err):
    bad_flux = ~np.isfinite(flux)
    bad_err = (~np.isfinite(err)) | (err <= 0)
    bad_pixels = bad_err | bad_flux
    return bad_pixels

# read data
all_set = Table.read("/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/allStar-v304.fits")




"""
# choose some of them
choose = []
a=0
for i in range(0,1000):
    a += np.random.randint(1,55)
    choose.append(a)

choose = np.array(choose)
all_set = all_set[choose]


# save the chosen star
output = open('choose_500.pkl', 'wb')
pickle.dump(choose, output)
output.close()

"""
# open the chosen star
pkl_file = open('choose_500.pkl', 'rb')
choose = pickle.load(pkl_file)
pkl_file.close()

all_set = all_set[choose]

# data
all_set_flux = []
all_set_ivar = []
all_set_error = []
test_labels_all = []
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
pickle.dump(norm_test_flux, output)
output.close()


output = open('testing_ivar_500.pkl', 'wb')
pickle.dump(tr_ivar, output)
output.close()

output = open('nor_testing_ivar_500.pkl', 'wb')
pickle.dump(norm_test_ivar, output)
output.close()

output = open('FiberID.pkl', 'wb')
pickle.dump(FiberID, output)
output.close()

#####################################
# Then Let's do it one by one

# training
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

########################
# optimize the spectrum

inf_label_500 = model.fit(norm_test_flux,norm_test_ivar)

v_500 = model.vectorizer.get_label_vector(inf_label_500)
# v_500 = model.vectorizer.get_label_vector(test_label)

inf_flux_500 = np.dot(v_500,model.theta.T)


opt_flux_500,parameters_500 = model.fitting_spectrum_parameters_single(norm_test_flux,norm_test_ivar,inf_flux_500)


#check
print("check")
print(opt_flux_500.shape,inf_flux_500.shape,FiberID.shape,parameters_500.shape)


###################
# Let's plot:

colors = ['b', 'g', 'r']
name = ["a","b","c"]

# Plot histogram
plt.hist(parameters_500,bins=20,stacked=True,color=colors,label=name)
plt.legend(prop={'size': 10})
plt.suptitle('Histogram of parameters a,b and c')
plt.xlabel('values of a, b and c', fontsize=20)
plt.ylabel('number of stars', fontsize=20)

plt.show()



# plot the parameters vs metadata


font = {'weight': 'bold','size': 30}
matplotlib.rc('font', **font)
fig = plt.figure()

plt.plot(FiberID,parameters_500[:,0], "bo", label="a", linewidth=3.0)
plt.plot(FiberID,parameters_500[:,1], "go", label="b", linewidth=3.0)
plt.plot(FiberID,parameters_500[:,2], "ro", label="c", linewidth=3.0)

plt.legend(loc="best")

fig.suptitle('a, b, c against mean fiber number', fontsize=40)
plt.xlabel('fiber number', fontsize=38)
plt.ylabel('parameters a, b, c', fontsize=36)
plt.show()







