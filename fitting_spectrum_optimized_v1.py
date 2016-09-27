# import
from numpy import genfromtxt
import matplotlib.pyplot as plt
import random
import matplotlib
import numpy as np
import scipy.optimize as optimization
import csv
import pickle


# import the data

pkl_file = open('normalized_flux.pkl', 'rb')
norm_tr_flux = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('inferred_spectrum.pkl', 'rb')
inferred_flux = pickle.load(pkl_file)
pkl_file.close()


# fit a b c for a single pixel/wave length

# define a function for fitting
def func(m, a, b, c):
    (x,y,z) = m
    return a * x + b * y + c * z

# define the class for fitting spectrum by nearby spectrum
class fit_neighborhood:

    def __init__(self,n,i):
        self.n = np.array(n)
        self.i = np.array(i)

    def fitting_spectrum_parameters(self):
        nor = self.n
        inf = self.i
        n_pixel = nor[0, :].size
        n_star = inf[:, 0].size
        nor = norm_tr_flux
        inf = inferred_flux
        n_pixel = nor[0, :].size
        n_star = inf[:, 0].size

        x_data = inf
        y_data = inf
        z_data = inf

        one = np.ones(n_star)
        x_data = np.delete(x_data, (0), axis=1)
        x_data = np.c_[one, x_data]

        z_data = np.delete(z_data, (n_star - 1), axis=1)
        z_data = np.c_[z_data,one]

        x_data = x_data.ravel()
        y_data = y_data.ravel()
        z_data = z_data.ravel()
        nor = nor.ravel()

        # The trick is here. give x,y,z equal positions
        x0 = np.array([0, 0, 0])
            # fit
        popt, pcov= optimization.curve_fit(func, (x_data, y_data, z_data),nor, x0,method="lm")
        parameters = popt


        # infer flux from the parameters
        opt_flux = np.zeros([n_pixel])
        for star in range(0, n_star):
            star_f = []
            for i in range(0, n_pixel):
                if i > 0 and i< n_pixel - 1:
                    x_data = inf[star, i - 1]
                    y_data = inf[star, i]
                    z_data = inf[star, i+1]
                    opt_f = x_data * parameters[0] + y_data * parameters[1] + z_data * parameters[2]
                elif i == 0:
                    x_data = 1
                    y_data = inf[star, i]
                    z_data = inf[star, i + 1]
                    opt_f = parameters[0] + y_data * parameters[1] + z_data * \
                                                                               parameters[2]
                elif i == n_pixel-1:
                    x_data = inf[star, i-1]
                    y_data = inf[star, i]
                    z_data = 1
                    opt_f = x_data * parameters[0] + y_data * parameters[1] +parameters[2]
                star_f.append(opt_f)
            opt_flux = np.c_[opt_flux,star_f]
        opt_flux = opt_flux.T
        opt_flux = np.delete(opt_flux,0,axis=0)
        print(parameters[0],parameters[1],parameters[2])

        output = open('optimized_spectrum.pkl', 'wb')
        pickle.dump(opt_flux, output)
        output.close()
        return opt_flux



trial = fit_neighborhood(norm_tr_flux,inferred_flux)
print(trial.fitting_spectrum_parameters().shape)










