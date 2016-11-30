
import os
import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle
import plot_class_4_star_v2_1107

# load path

pkl_file = open('n_900_path_fits.pkl', 'rb')
path_fits = pickle.load(pkl_file)
pkl_file.close()


pkl_file = open('n_900_path_flux.pkl', 'rb')
path_flux = pickle.load(pkl_file)
pkl_file.close()


N = len(path_fits)
print(N)

parameters = []


for i in range(0,N):
    print("loading star %d"%(i+1))
    star_i = fits.open(path_fits[i])
    parameters.append(star_i[4].data[0])
    print(star_i[4].data[0].shape)

parameters = np.array(parameters)

parameters_500 = parameters

big_a_c = (parameters_500[:,2]-parameters_500[:,0]).argsort()[0:4]

big_c_a = (parameters_500[:,0]-parameters_500[:,2]).argsort()[0:4]


print(parameters_500[big_a_c])


# plot
model = plot_class_4_star_v2_1107.plot_4()
model.plot_name_4(path_fits[big_a_c])


