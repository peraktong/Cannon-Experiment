# In this version, we choose another objective function and improve the fitting by using a better optimization method.



# import
import matplotlib
import scipy.optimize as optimization
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import pickle
import AnniesLasso as tc
from scipy.stats import chisquare
from numpy.linalg import inv

#import the data

training_set = Table.read("reference_labels.csv")

pkl_file = open('training_set_flux.pkl', 'rb')
training_set_flux = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('training_set_ivar.pkl', 'rb')
training_set_ivar = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('training_set_error.pkl', 'rb')
error = pickle.load(pkl_file)
pkl_file.close()





# train the model
model = tc.L1RegularizedCannonModel(
    training_set, training_set_flux, training_set_ivar,threads=8)

model.s2 = 0
model.regularization = 0
model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(training_set,
    tc.vectorizer.polynomial.terminator(("Teff_{corr}", "logg_{corr}", "[M/H]_{corr}"), 2))

model.train()

# obtain the spectrum data
norm_tr_flux = training_set_flux

inferred_labels = model.fit_labelled_set()
v = model.vectorizer.get_label_vector(inferred_labels)
inferred_flux = np.dot(v,model.theta.T)
print(v.shape)



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


    def fitting_spectrum_parameters(self,ivar):
        nor = self.n
        inf = self.i
        n_pixel = nor[0, :].size
        n_star = inf[:, 0].size
        nor = norm_tr_flux
        inf = inferred_flux
        n_pixel = nor[0, :].size
        n_star = inf[:, 0].size
        one = np.ones(n_star)

        for pixel in range(0, n_pixel):
            if pixel ==0 :
                x_data = one
                y_data = inf[:,pixel]
                z_data = inf[:,pixel+1]

            elif pixel>0 and pixel < n_pixel-1:
                x_data = np.c_[x_data,inf[:,pixel-1]]
                y_data = np.c_[y_data,inf[:,pixel]]
                z_data = np.c_[z_data,inf[:,pixel+1]]

            elif pixel == n_pixel-1:
                x_data = np.c_[x_data, inf[:, pixel - 1]]
                y_data = np.c_[y_data, inf[:, pixel]]
                z_data = np.c_[z_data, one]




        """
                output = open('x_data.pkl', 'wb')
        pickle.dump(x_data, output)
        output.close()

        output = open('y_data.pkl', 'wb')
        pickle.dump(y_data, output)
        output.close()

        output = open('z_data.pkl', 'wb')
        pickle.dump(z_data, output)
        output.close()

        """


        # fit
        x_data_r = x_data.ravel()
        y_data_r = y_data.ravel()
        z_data_r = z_data.ravel()
        nor_r = nor.ravel()
        error_r = error.ravel()

        # It's not good. let's do it one star each time.

        n_star = len(nor[:, 0])

        para = [0, 0, 0]

        for p in range(0, n_star):

            x_data_p = x_data[p, :]
            y_data_p = y_data[p, :]
            z_data_p = z_data[p, :]
            nor_p = nor[p, :]
            ivar_p = ivar[p, :]

            # construct
            ivar_r = ivar_p.ravel()
            ni = len(ivar_r)
            print(p)
            c = np.zeros((ni, ni))

            for i in range(0, ni):
                c[i, i] = ivar_r[i]

            y = nor_p.ravel()
            a = np.c_[np.c_[x_data_p.ravel(), y_data_p.ravel()], z_data_p.ravel()]

            para_p = np.dot(inv(np.dot(np.dot(a.T, c), a)), np.dot(np.dot(a.T, c), y))
            para += para_p

        para = para / n_star

        parameters = para


        # infer flux from the parameters
        # add scale factor

        opt_flux = (parameters[0]*x_data+parameters[1]*y_data+parameters[2]*z_data)

        self.optimized_spectrum = opt_flux

        print(parameters[0],parameters[1],parameters[2],parameters[0]+parameters[1]+parameters[2])
        #save the optimized spectrum

        output = open('optimized_spectrum_v5.pkl', 'wb')
        pickle.dump(opt_flux, output)
        output.close()



        return opt_flux

    def plot_spectrum(self):
        nor = self.n
        inf_old = self.i
        inf_opt = self.optimized_spectrum

        p = np.random.randint(0, 547)
        # plot
        font = {'weight': 'bold', 'size': 30}
        matplotlib.rc('font', **font)
        fig = plt.figure()

        plt.plot(wl, nor[p, :], "b",label="Observed", linewidth=3.0)
        plt.plot(wl, inf_old[p, :], "g",label="Synthesized", linewidth=3.0)
        plt.plot(wl, inf_opt[p, :], "r",label="Optimized", linewidth=3.0)
        plt.legend()
        fig.suptitle('Comparison of the spectrum', fontsize=40)
        plt.xlabel('wave length', fontsize=38)
        plt.ylabel('Spectrum', fontsize=36)

        axes = plt.gca()
        axes.set_xlim([15660, 15780])
        axes.set_ylim([0.8, 1.1])

        plt.show()

    def plot_chi_square(self,ivar):

        nor = self.n
        inf_old = self.i
        inf_opt = self.optimized_spectrum

        chi_old = np.array([])
        p_old = np.array([])

        chi_opt = np.array([])
        p_opt = np.array([])

        for p in range(0, nor[:, 0].size):
            ivar_r = ivar[p, :]
            ni = len(ivar_r)

            c = np.zeros((ni, ni))

            for i in range(0, ni):
                c[i, i] = ivar_r[i]

            # correct chi-squared
            a_old = np.dot(np.dot(nor[p, :] - inf_old[p, :], c), (nor[p, :] - inf_old[p, :]).T)
            a_opt = np.dot(np.dot(nor[p, :] - inf_opt[p, :], c), (nor[p, :] - inf_opt[p, :]).T)
            chi_old = np.append(chi_old, a_old)
            # p_old.append(b_old)
            chi_opt = np.append(chi_opt, a_opt)
            # p_old.append(b_opt)

        # plot


        font = {'weight': 'bold', 'size': 30}
        matplotlib.rc('font', **font)
        fig = plt.figure()

        plt.plot(chi_old, "g",label="Synthesized", linewidth=3.0)
        plt.plot(chi_opt, "r", label="Optimized", linewidth=3.0)
        plt.legend()
        fig.suptitle('compare chi-square', fontsize=40)
        plt.xlabel('stellar', fontsize=38)
        plt.ylabel('chi-square', fontsize=36)


        plt.show()

    def plot_chi_square_single(self,ivar):

        nor = self.n
        inf_old = self.i
        inf_opt = self.optimized_spectrum

        p =np.random.randint(0,547)

        ivar_r = ivar[p, :]
        ni = len(ivar_r)

        c = np.zeros((ni, ni))

        for i in range(0, ni):
            c[i, i] = ivar_r[i]

        # correct chi-squared
        a_old = np.dot(np.dot(nor[p, :] - inf_old[p, :], c), (nor[p, :] - inf_old[p, :]).T)
        a_opt = np.dot(np.dot(nor[p, :] - inf_opt[p, :], c), (nor[p, :] - inf_opt[p, :]).T)

        print(a_old,a_opt)




trial = fit_neighborhood(norm_tr_flux,inferred_flux)
print(trial.fitting_spectrum_parameters(training_set_ivar).shape)

trial.plot_spectrum()
#trial.plot_chi_square(training_set_ivar)
trial.plot_chi_square_single(training_set_ivar)
print(trial.optimized_spectrum.shape)