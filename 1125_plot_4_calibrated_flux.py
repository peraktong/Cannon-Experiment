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

import plot_class_4_star_v2_1107








# import data
pkl_file = open('n_delta_chi_900.pkl', 'rb')
delta_chi_500 = pickle.load(pkl_file)
pkl_file.close()



pkl_file = open('n_testing_labels_900.pkl', 'rb')
test_labels_500 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n_parameters_900.pkl', 'rb')
parameters_500 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n_FiberID_900.pkl', 'rb')
FiberID = pickle.load(pkl_file)
pkl_file.close()


pkl_file = open('n_tr_ID_900.pkl', 'rb')
tr_ID_500 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n_inf_label_900.pkl', 'rb')
inf_labels_900 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open("testing_set_900.pkl", 'rb')
testing_set = pickle.load(pkl_file)
pkl_file.close()


# flux,ivar and error

pkl_file = open("n_testing_flux_900.pkl", 'rb')
testing_flux = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open("n_testing_ivar_900.pkl", 'rb')
testing_ivar = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open("n_inf_flux_900.pkl", 'rb')
inf_flux_900 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open("n_testing_error_900.pkl", 'rb')
testing_error = pickle.load(pkl_file)
pkl_file.close()


## check error
"""
plt.plot(testing_error[1,:],"r")

axes = plt.gca()

axes.set_ylim([-0.01,0.01])


plt.show()

"""

## replace test_label with inf_label






#print(choose_800_final.shape)
#print(choose_800_final.shape,FiberID.shape)


# check shape
#print(parameters_500.shape,delta_chi_500.shape)

#### find the biggest four delta-chi and smallest four
N_star = len(delta_chi_500)


delta_chi_500 = np.array(delta_chi_500)

big_four = delta_chi_500.argsort()[N_star-4:N_star]
small_four = delta_chi_500.argsort()[0:4]

#print(big_four,small_four)
#print(delta_chi_500[big_four],delta_chi_500[small_four])

big_four_set = testing_set[big_four]
small_four_set = testing_set[small_four]

# print(delta_chi_500[big_four],delta_chi_500[small_four])

# find parameters closest to b =-1


a_500 = parameters_500[:,0]
a_500 = np.array(a_500)

b_500 = parameters_500[:,1]
b_500 = np.array(b_500)

c_500 = parameters_500[:,2]
c_500 = np.array(c_500)

plus = np.ones(N_star)
minus = np.ones(N_star)*(-1)

#print(plus,minus)


close_4_1 =((b_500-plus)**2).argsort()[0:4]
close_4_minus_1 =((b_500-minus)**2).argsort()[0:4]
biggest_a_c_b = (a_500+c_500-b_500).argsort()[N_star-4:N_star]

# find flux and save them
"""
print(close_4_1,close_4_minus_1)
print(b_500[close_4_1],b_500[close_4_minus_1])
print(tr_ID_500[close_4_1])
print(tr_ID_500[close_4_minus_1])
print(parameters_500[close_4_1],parameters_500[close_4_minus_1])

print(close_4_1,close_4_minus_1)
print(b_500[close_4_1],b_500[close_4_minus_1])
print(tr_ID_500[close_4_1])
print(tr_ID_500[close_4_minus_1])
print(parameters_500[close_4_1],parameters_500[close_4_minus_1])

"""




#print(testing_set[close_4_1]["APOGEE_ID"],testing_set[close_4_minus_1]["APOGEE_ID"])

close_1_set = testing_set[close_4_1]
close_minus_1_set = testing_set[close_4_minus_1]
biggest_a_c_b_set = testing_set[biggest_a_c_b]
## name

"""
print(big_four_set,small_four_set)

print(close_1_set,close_minus_1_set)

"""

# also choose four stars with biggest |a|-|c| and |c|-|a|

# a,c

"""
big_a_c = (np.absolute(parameters_500[:,0])-np.absolute(parameters_500[:,2])).argsort()[N_star-4:N_star]
big_a_c_set = testing_set[big_a_c]

big_c_a = (np.absolute(parameters_500[:,2])-np.absolute(parameters_500[:,0])).argsort()[N_star-4:N_star]
big_c_a_set = testing_set[big_c_a]

"""

"""
big_a_c = (parameters_500[:,0]-parameters_500[:,2]).argsort()[N_star-4:N_star]
big_a_c_set = testing_set[big_a_c]
"""

big_a_c = (parameters_500[:,2]-parameters_500[:,0]).argsort()[0:4]
big_a_c_set = testing_set[big_a_c]

big_c_a = (parameters_500[:,0]-parameters_500[:,2]).argsort()[0:4]
big_c_a_set = testing_set[big_c_a]









"""
print(big_a_c,big_c_a)
print(testing_set[big_a_c],testing_set[big_c_a])
print(big_a_c_set,big_c_a_set)

print(parameters_500[big_a_c],parameters_500[big_c_a])

print(testing_set[big_a_c]["APOGEE_ID"],testing_set[big_c_a]["APOGEE_ID"])


"""


"""
plt.plot(test_labels[:,1],"ro")
plt.show()

"""





#save

output = open('big_four_delta_set.pkl', 'wb')
pickle.dump(big_four_set, output)
output.close()

output = open('small_four_delta_set.pkl', 'wb')
pickle.dump(small_four_set, output)
output.close()

output = open('close_1_set.pkl', 'wb')
pickle.dump(close_1_set, output)
output.close()

output = open('close_minus_1_set.pkl', 'wb')
pickle.dump(close_minus_1_set, output)
output.close()

output = open('biggest_a_c_b_set.pkl', 'wb')
pickle.dump(biggest_a_c_b_set, output)
output.close()

output = open('big_a_c_set.pkl', 'wb')
pickle.dump(big_a_c_set, output)
output.close()

output = open('big_c_a_set.pkl', 'wb')
pickle.dump(big_c_a_set, output)
output.close()

print(tr_ID_500[big_a_c])

##################
# plot:
###################




"""
model = plot_class_4_star_v2_1107.plot_4()
model.train()

print("big4")

model.plot_v3(flux=testing_flux[big_four],
           ivar=testing_ivar[big_four],
           error=testing_error[big_four],
           star_name=testing_set[big_four]["APOGEE_ID"],
              inf_flux=inf_flux_900[big_four],inf_label = inf_labels_900[big_four])

# plot small delta-chi
print("small4")

model.plot_v3(flux=testing_flux[small_four],
           ivar=testing_ivar[small_four],
           error=testing_error[small_four],
           star_name=testing_set[small_four]["APOGEE_ID"],
              inf_flux=inf_flux_900[small_four],inf_label=inf_labels_900[small_four])

# plot b = 1

print("b=1")

model.plot_v3(flux=testing_flux[close_4_1],
           ivar=testing_ivar[close_4_1],
           error=testing_error[close_4_1],
           star_name=testing_set[close_4_1]["APOGEE_ID"],
              inf_flux=inf_flux_900[close_4_1],inf_label=inf_labels_900[close_4_1])


# plot biggest a+c-b

print("a+c>b")

model.plot_v3(flux=testing_flux[biggest_a_c_b],
           ivar=testing_ivar[biggest_a_c_b],
           error=testing_error[biggest_a_c_b],
           star_name=testing_set[biggest_a_c_b]["APOGEE_ID"],
              inf_flux=inf_flux_900[biggest_a_c_b],inf_label=inf_labels_900[biggest_a_c_b])



# plot a>c

print("a>c")

print(parameters_500[big_a_c])

model.plot_v3(flux=testing_flux[big_a_c],
           ivar=testing_ivar[big_a_c],
           error=testing_error[big_a_c],
           star_name=testing_set[big_a_c]["APOGEE_ID"],
              inf_flux= inf_flux_900[big_a_c],inf_label=inf_labels_900[big_a_c])


# plot c>a
print("c>a")
model.plot_v3(flux=testing_flux[big_c_a],
           ivar=testing_ivar[big_c_a],
           error=testing_error[big_c_a],
           star_name=testing_set[big_c_a]["APOGEE_ID"],
              inf_flux=inf_flux_900[big_c_a],inf_label=inf_labels_900[big_c_a])


model.plot_calibrated(flux=testing_flux[big_a_c],
           ivar=testing_ivar[big_a_c],
           error=testing_error[big_a_c],
           star_name=testing_set[big_a_c]["APOGEE_ID"],
              inf_flux= inf_flux_900[big_a_c],inf_label=inf_labels_900[big_a_c])


# Data/APOGEE_DR10_Apstar/apStar-s3-2M06323568+0408166.fits
"""

# plot the first star in big_a_c

image_path = tr_ID_500[big_a_c][0]


image = fits.open(image_path,ignore_missing_end=True)

dat = Table.read(image_path)

flux = image[1].data
print(flux.shape)

# normalize:


flux = image[1].data
flux_err = image[2].data


# read

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

tr_ID = image_path
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


flux = norm_tr_flux
#plot:
tran =0.5

font = {'weight': 'bold', 'size': 13}
matplotlib.rc('font', **font)
fig = plt.figure()

fig.suptitle('Comparison of the spectrum for each visit', fontsize=30, horizontalalignment="center")

#plt.plot(wl,flux[0,:],alpha=tran)
#plt.plot(wl,flux[1,:],alpha=tran)
plt.plot(wl,flux[2,:],alpha=tran)
plt.plot(wl,flux[3,:],alpha=tran)
plt.plot(wl,flux[4,:],alpha=tran)
plt.plot(wl,flux[5,:],alpha=tran)
plt.plot(wl,flux[6,:],alpha=tran)
plt.plot(wl,flux[7,:],alpha=tran)
plt.plot(wl,flux[8,:],alpha=tran)
plt.plot(wl,flux[9,:],alpha=tran)
plt.plot(wl,flux[10,:],alpha=tran)

plt.xlabel('wave length $\AA$', fontsize=30)
plt.ylabel('Spectrum', fontsize=30)

axes = plt.gca()
axes.set_xlim([15660,15780])

#axes.set_xlim([16160,16280])
axes.set_ylim([0.6,1.41])
axes.set_yticks(np.arange(0.6,1.41,0.1))

plt.show()

















































