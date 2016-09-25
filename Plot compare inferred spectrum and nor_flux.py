from numpy import genfromtxt
import matplotlib.pyplot as plt
import random
import matplotlib
import numpy as np

norm_tr_flux = genfromtxt('normalize flux.csv', delimiter=',')
inferred_flux = genfromtxt('inferred_spectrum.csv', delimiter=',')
inferred_flux_opt = genfromtxt('optimized_flux.csv', delimiter=',')
wl = genfromtxt("wl")

pixel =np.arange(8575)
# choose one
p = np.random.randint(0,547)

# plot
font = {'weight': 'bold',
            'size': 30}

matplotlib.rc('font', **font)
fig = plt.figure()

plt.plot(wl,norm_tr_flux[p,:],wl,inferred_flux[p,:],wl,inferred_flux_opt[p,:],linewidth =5.0)
fig.suptitle('Comparison of the spectrum', fontsize=40)
plt.xlabel('wave length', fontsize=38)
plt.ylabel('Spectrum', fontsize=36)

axes = plt.gca()
axes.set_xlim([15660,15780])
axes.set_ylim([0.8,1.1])

plt.show()
