# import
import matplotlib
import scipy.optimize as optimization
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import pickle
import AnniesLasso as tc
from scipy.stats import chisquare

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

        self.optimized_spectrum = opt_flux

        print(parameters[0],parameters[1],parameters[2])
        #save the optimized spectrum
        output = open('optimized_spectrum_v2.pkl', 'wb')
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

        plt.plot(wl, nor[p, :], "b", wl, inf_old[p, :], "g", wl, inf_opt[p, :], "r", linewidth=3.0)
        fig.suptitle('Comparison of the spectrum', fontsize=40)
        plt.xlabel('wave length', fontsize=38)
        plt.ylabel('Spectrum', fontsize=36)

        axes = plt.gca()
        axes.set_xlim([15660, 15780])
        axes.set_ylim([0.8, 1.1])

        plt.show()

    def plot_chi_square(self):

        nor = self.n
        inf_old = self.i
        inf_opt = self.optimized_spectrum

        chi_old = np.array([])
        p_old = np.array([])

        chi_opt = np.array([])
        p_opt = np.array([])

        for p in range(0, nor[:, 0].size):
            nor_i = nor[p, :]
            inf_old_i = inf_old[p, :]
            inf_opt_i = inf_opt[p, :]
            a_old, b_old = chisquare(nor_i, inf_old_i)
            a_opt, b_opt = chisquare(nor_i, inf_opt_i)
            chi_old = np.append(chi_old, a_old)
            # p_old.append(b_old)
            chi_opt = np.append(chi_opt, a_opt)
            # p_old.append(b_opt)

        # plot


        font = {'weight': 'bold', 'size': 30}
        matplotlib.rc('font', **font)
        fig = plt.figure()

        plt.plot(chi_old, "r", chi_opt, "g", linewidth=5.0)

        fig.suptitle('compare chi-square', fontsize=40)
        plt.xlabel('stellar', fontsize=38)
        plt.ylabel('chi-square', fontsize=36)

        # axes = plt.gca()
        # axes.set_xlim([15660,15780])
        # axes.set_xlim([16160,16280])
        # axes.set_ylim([-200,200])

        plt.show()



trial = fit_neighborhood(norm_tr_flux,inferred_flux)
print(trial.fitting_spectrum_parameters().shape)

trial.plot_spectrum()
trial.plot_chi_square()
print(trial.optimized_spectrum.shape)