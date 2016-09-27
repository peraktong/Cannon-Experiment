from numpy import genfromtxt
import matplotlib.pyplot as plt
import random
import matplotlib
import numpy as np
import pickle


# read data
pkl_file = open('normalized_flux.pkl', 'rb')
nor = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('inferred_spectrum.pkl', 'rb')
inf_old = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('optimized_spectrum.pkl', 'rb')
inf_opt = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()


p = np.random.randint(0,547)

# plot
font = {'weight': 'bold','size': 30}
matplotlib.rc('font', **font)
fig = plt.figure()

plt.plot(wl,nor[p,:],"b",wl,inf_old[p,:],"g",wl,inf_opt[p,:],"r",linewidth =3.0)
fig.suptitle('Comparison of the spectrum', fontsize=40)
plt.xlabel('wave length', fontsize=38)
plt.ylabel('Spectrum', fontsize=36)

axes = plt.gca()
axes.set_xlim([15660,15780])
axes.set_ylim([0.8,1.1])

plt.show()
