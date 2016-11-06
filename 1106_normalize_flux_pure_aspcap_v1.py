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

# data

all_set_flux = np.zeros(8575)
all_set_ivar = np.zeros(8575)
all_set_error = np.zeros(8575)


# choose some of them
choose = []
a=0
for i in range(0,1600):
    a += np.random.randint(1,30)
    choose.append(a)

choose = np.array(choose)
all_set = all_set[choose]
# Now all set only have 1600 stars


# save the chosen star and the new all-set. Let's name it as testing set.
output = open('n_choose_800.pkl', 'wb')
pickle.dump(choose, output)
output.close()




N = len(all_set)
test_labels_all = []
tr_ID = []
FiberID = []
choose_800_final = []
fail = 0
star=0

# choose some flux and normalize them first, then filter the mean inverse variance.


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

        flux = image[1].data[0]
        flux_err = image[2].data[0]

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


        try:
            all_set_flux = np.vstack((all_set_flux, flux))
            all_set_ivar = np.vstack((all_set_ivar, ivar))
            all_set_error = np.vstack((all_set_error, error))
            test_labels_all_i = [row["TEFF"], row["LOGG"], row["METALS"]]
            test_labels_all.append(test_labels_all_i)
            choose_800_final.append(i)


        except ValueError:
            print("opts fail")
            fail+=1

        else:
            star += 1
            # Fiber number
            m = dat[0]["FIBER"]
            try:
                FiberID.append(sum(m) / len(m))
            except TypeError:
                FiberID.append(m)

            print(m)
            print(type(m))
            tr_ID.append(image_path)




all_set_flux = all_set_flux[1:star+1,:]
all_set_ivar = all_set_ivar[1:star+1,:]
all_set_error = all_set_error[1:star+1,:]
test_labels_all =np.array(test_labels_all)
tr_ID = np.array(tr_ID)
FiberID = np.array(FiberID)
choose_800_final = np.array(choose_800_final)


#print("check")
#print(all_set_flux.shape,
#      all_set_ivar.shape,all_set_error.shape,
#      test_labels_all.shape,tr_ID.shape,FiberID.shape,star)



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



### second filter for inverse mean square

keep = np.ones(star, dtype=bool)


for i in range(0,star):
    ivar = norm_tr_ivar[i,:]
    if np.mean(ivar) > 10000 and np.mean(ivar) < 40000:
        keep[i]=1

    else:
        keep[i] = 0
        print("mean ivar too big/small")
        fail += 1
        star -= 1

##re-calculate nor:
norm_flux = norm_tr_flux[keep]
norm_ivar =norm_test_ivar[keep]

# normalize the error
all_set_error = all_set_error[keep]
norm_error = all_set_error*all_set_flux[keep]/norm_flux

test_labels_all =test_labels_all[keep]
tr_ID = tr_ID[keep]
FiberID = FiberID[keep]
choose_800_final = choose_800_final[keep]

## check shape

print(norm_flux.shape,norm_ivar.shape,norm_error.shape,
      test_labels_all.shape,tr_ID.shape,
FiberID.shape,choose_800_final.shape,star,fail)






# save the testing set

print("number of fail",fail)

output = open('testing_set_900.pkl', 'wb')
pickle.dump(all_set[choose_800_final], output)
output.close()

print("There are %d star in the testing set"%len(choose_800_final))



output = open('n_testing_labels_900.pkl', 'wb')
pickle.dump(test_labels_all, output)
output.close()

output = open('n_tr_ID_900.pkl', 'wb')
pickle.dump(tr_ID, output)
output.close()

output = open('choose_900.pkl', 'wb')
pickle.dump(choose_800_final, output)
output.close()

output = open('n_FiberID_900.pkl', 'wb')
pickle.dump(FiberID, output)
output.close()


output = open('n_testing_flux_900.pkl', 'wb')
pickle.dump(norm_flux, output)
output.close()

output = open('n_testing_ivar_900.pkl', 'wb')
pickle.dump(norm_ivar, output)
output.close()

output = open('n_testing_error_900.pkl', 'wb')
pickle.dump(norm_error, output)
output.close()





##plot and check

#Let's plot them

p = np.random.randint(0,star)

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




