
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
inf_labels = pickle.load(pkl_file)
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

pkl_file = open("n_testing_error_900.pkl", 'rb')
testing_error = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open("n_inf_flux_900.pkl", 'rb')
inf_flux = pickle.load(pkl_file)
pkl_file.close()



big_a_c = (parameters_500[:,2]-parameters_500[:,0]).argsort()[0:4]
big_a_c_set = testing_set[big_a_c]



##################
# plot:
###################
model = plot_class_4_star_v2_1107.plot_4()
model.train()



## make 4 fake stars with data spectrum equal to inf_flux moved 2 pixel left
inf_4 = inf_flux[big_a_c]

n_pixel = inf_4[0, :].size

one = np.ones(4)

# building matrix x,y,z

x_data = np.c_[inf_4,one]
x_data = x_data[:,1:n_pixel+1]


nor_fake_4 = x_data



# plot a>c

print("a>c")


model.plot_v3(flux=nor_fake_4,
           ivar=testing_ivar[big_a_c],
           error=testing_error[big_a_c],
           star_name=testing_set[big_a_c]["APOGEE_ID"],
              inf_flux=inf_4,inf_label=inf_labels[big_a_c])


























































