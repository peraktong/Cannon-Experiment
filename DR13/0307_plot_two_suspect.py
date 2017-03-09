
import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
import pickle

from TheCannon_2 import dataset,apogee
from TheCannon_2 import model


pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()



name = ["/Users/caojunzhi/Desktop/Data/suspect_2_flux/2M08160368+3055091.fits","/Users/caojunzhi/Desktop/Data/suspect_2_flux/2M08165743+3256475.fits"]



"""


star = fits.open(name[0])

print(star[3].data)
print(star[6].data)
print(star[8].data[:,3:6])



"""





def plot_visit_sim_old_customized(name):
    # set font size

    font = {'family': 'normal',
            'weight': 'bold',
            'size': 8}

    matplotlib.rc('font', **font)

    # only choose individual visits:
    star = fits.open(name)

    flux = star[0].data
    ivar = star[1].data
    inf_flux = star[2].data
    parameters =star[3].data
    inf_labels = star[4].data
    inf_flux_sim = star[7].data
    inf_labels_sim = star[8].data[:,0:3]
    parameters_sim = star[8].data[:,3:6]

    parameters_new = star[6].data

    N = len(flux[:,0])
    for i in range(0, N):

        # left
        plt.subplot(N, 2, 2 * i + 1)

        name = name.replace("/Users/caojunzhi/Desktop/Data/suspect_2_flux/", "")
        name = name.replace(".fits", "")


        ## Add inf labels

        plt.step(wl, flux[i,:], "k", label="Data flux", linewidth=0.7, alpha=0.8)
        #plt.errorbar(wl,ivar[i,:]**(-0.5),alpha = 0.1)

        plt.plot(wl, inf_flux[i,:], "b", label="Inferred flux from fitting separately Teff=%.2f K logg=%.2f Fe/H =%.2f" % (
        inf_labels[i,0], inf_labels[i,1], inf_labels[i,2]), linewidth=0.7, alpha=0.4)

        plt.plot(wl, inf_flux_sim[i,:], "g", label="Inferred flux from fitting simultaneously Teff=%.2f K logg=%.2f Fe/H =%.2f" % (
        inf_labels_sim[i,0], inf_labels_sim[i,1], inf_labels_sim[i,2]), linewidth=0.7, alpha=0.4)

        plt.xlabel("Wave length $\AA$")
        plt.ylabel("Flux", fontsize=20)
        axes = plt.gca()
        axes.set_xlim([15660,15780])
        # share x axis

        if i == N - 1:
            non = 1
        else:
            axes.set_xticklabels([])

        # axes.set_xlim([16160,16280])
        axes.set_ylim([0.5, 2])
        axes.set_yticks(np.arange(0.5, 1.99, 0.5))
        plt.legend()

        # inf

        plt.subplot(N, 2, 2 * i + 2)

        ## Add inf labels

        plt.step(wl, flux[i,:], "k", linewidth=0.7, alpha=0.8)
        #plt.errorbar(wl,ivar[i,:]**(-0.5),alpha = 0.1)

        ai = parameters[i,0]
        bi = parameters[i,1]
        ci = parameters[i,2]

        plt.plot(wl, inf_flux[i,:], "b", label="Velocity from fitting separately a=%.2f b=%.2f c =%.2f $RV_{shift}$=%.2f $m/s$" % (
        ai,bi,ci,(ci-ai)/(ai+bi+ci)*4144.68), linewidth=0.7, alpha=0.4)

        ai = parameters_sim[i,0]
        bi = parameters_sim[i,1]
        ci = parameters_sim[i,2]

        plt.plot(wl, inf_flux_sim[i,:], "g", label="Velocity from fitting simultaneously a=%.2f b=%.2f c =%.2f $RV_{shift}$=%.2f $m/s$" % (
        ai, bi, ci, (ci - ai) / (ai + bi + ci) * 4144.68), linewidth=0.7, alpha=0.4)

        # add absorption lines:

        ai = parameters_new[i,0]
        bi = parameters_new[i,1]
        ci = parameters_new[i,2]

        plt.plot(0,0, "ko",markersize=0.5,
                 label="Velocity from fitting absorption lines a=%.2f b=%.2f c =%.2f $RV_{shift}$=%.2f $m/s$" % (
                     ai, bi, ci, (ci - ai) / (ai + bi + ci) * 4144.68), linewidth=0.7, alpha=0.4)

        plt.xlabel("Wave length $\AA$")
        plt.ylabel("Flux", fontsize=20)
        axes = plt.gca()
        axes.set_xlim([16160,16280])
        # share x axis

        if i == N - 1:
            non = 1
        else:
            axes.set_xticklabels([])

        # axes.set_xlim([16160,16280])
        axes.set_ylim([0.5, 2])
        axes.set_yticks(np.arange(0.5, 1.99, 0.5))
        plt.legend()


    plt.suptitle("The fitting result of star %s"%name, fontsize=20)

    # share x
    plt.subplots_adjust(hspace=.0)

    plt.show()

plot_visit_sim_old_customized(name[1])

