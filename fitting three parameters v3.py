# import
from numpy import genfromtxt
import matplotlib.pyplot as plt
import random
import matplotlib
import numpy as np
import scipy.optimize as optimization
import csv


# import the data
norm_tr_flux = genfromtxt('normalize flux.csv', delimiter=',')
inferred_flux = genfromtxt('inferred_spectrum.csv', delimiter=',')
wl = genfromtxt("wl")


# fit a b c for a single pixel/wave length

# define a function for fitting
def func((x, y, z), a, b, c):
    return a * x + b * y + c * z


class fit_neighborhood:

    def __init__(self,n,i):
        self.n = np.array(n)
        self.i = np.array(i)


    def fitting_spectrum_parameters(self):
        nor = self.n
        inf = self.i
        n_pixel = nor[0, :].size
        n_star = inf[:, 0].size

        parameters = np.array([0,1,0])
        for i in range(0,n_pixel-1):
            if i > 0 & i < n_pixel-1:
                x_data = inf[:, i - 1]
                y_data = inf[:, i]
                z_data = inf[:, i + 1]
            elif i == 0:
                x_data = np.ones(n_star)
                y_data = inf[:, i]
                z_data = inf[:, i + 1]
            else:
                x_data = inf[:, i - 1]
                y_data = inf[:, i]
                z_data = np.ones(n_star)
            x0 = np.array([0, 1, 0])
            # fit
            popt, pcov= optimization.curve_fit(func, (x_data, y_data, z_data),nor[:,i], x0)
            #least square fit
            #popt, pcov = optimization.leastsq(func, x0,args=((x_data, y_data, z_data),nor[:,i]))
            parameters = np.c_[parameters,popt]
        self.parameters = parameters
        with open("fitting parameters-3.csv", "wb") as f:
            writer = csv.writer(f)
            writer.writerows(parameters)

        # infer flux from the parameters
        opt_flux = np.zeros([n_pixel])
        for star in range(0, n_star):
            star_f = []
            for i in range(0, n_pixel):
                if i > 0 and i< n_pixel - 1:
                    x_data = inf[star, i - 1]
                    y_data = inf[star, i]
                    z_data = inf[star, i+1]
                    opt_f = x_data * parameters[0, i] + y_data * parameters[1, i] + z_data * parameters[2, i]
                elif i == 0:
                    x_data = 1
                    y_data = inf[star, i]
                    z_data = inf[star, i + 1]
                    opt_f = parameters[0, i] + y_data * parameters[1, i] + z_data * \
                                                                               parameters[2, i]
                elif i == n_pixel-1:
                    x_data = inf[star, i-1]
                    y_data = inf[star, i]
                    z_data = 1
                    opt_f = x_data * parameters[0, i] + y_data * parameters[1, i] +parameters[2,i]
                star_f.append(opt_f)
            opt_flux = np.c_[opt_flux,star_f]
        opt_flux = opt_flux.T
        opt_flux = np.delete(opt_flux,0,axis=0)



        with open("optimized_flux.csv", "wb") as f:
            writer = csv.writer(f)
            writer.writerows(opt_flux)
        return opt_flux







trial = fit_neighborhood(norm_tr_flux,inferred_flux)
print trial.fitting_spectrum_parameters().shape









