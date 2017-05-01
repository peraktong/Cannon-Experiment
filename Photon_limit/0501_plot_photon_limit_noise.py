import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
import pickle
from os.path import isfile, join
from os import listdir

import astropy.io.fits as ts
from TheCannon import dataset,apogee
from astropy.table import Table
import numpy as np
import os
from astropy.io import fits
import pickle
import AnniesLasso_2 as tc
import time
from os import listdir
from os.path import isfile, join

from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time


from astropy.time import Time


# do not need to train the Cannon first:
training_set_path = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/reference_labels.csv"
training_set_spectrum_dir = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/Data/"

training_set = Table.read("reference_labels.csv")


def get_pixmask(flux, err):
    bad_flux = ~np.isfinite(flux)
    bad_err = (~np.isfinite(err)) | (err <= 0)
    bad_pixels = bad_err | bad_flux
    return bad_pixels


N = len(training_set)
keep = np.ones(N, dtype=bool)

training_set_flux = []
training_set_ivar = []
training_set_error = []

for i, row in enumerate(training_set):
    image_path = os.path.join(training_set_spectrum_dir, row["ID"])
    if not os.path.exists(image_path):
        print("{}/{} could not be found: {}".format(i + 1, N, image_path))
        keep[i] = False
        continue
    print("{}/{}: {}".format(i + 1, N, image_path))
    image = fits.open(image_path)
    flux = image[1].data
    flux_err = image[2].data
    badpix = get_pixmask(flux, flux_err)
    ivar = 1.0 / flux_err ** 2
    error = flux_err
    # badpix is a array and the length is 8575
    flux[badpix] = np.median(flux)
    ivar[badpix] = 0.0
    training_set_flux.append(flux)
    training_set_ivar.append(ivar)
    training_set_error.append(error)

training_set_flux = np.array(training_set_flux)
training_set_ivar = np.array(training_set_ivar)
training_set_error = np.array(training_set_error)

assert all(keep)

# Construct model.
model = tc.L1RegularizedCannonModel(
    training_set, training_set_flux, training_set_ivar, threads=8)
model.s2 = 0
model.regularization = 0
model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(training_set,
                                                                tc.vectorizer.polynomial.terminator(
                                                                    ("Teff_{corr}", "logg_{corr}", "[M/H]_{corr}"),
                                                                    2))

# model.train()



# read data


pkl_file = open('/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/dr13_fits.pkl',
                'rb')
path_fits = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/dr13_ori.pkl',
                'rb')
path_ori = pickle.load(pkl_file)
pkl_file.close()


pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()



# choose a number



def plot_star(path,ori):

    # read it
    star = fits.open(path)

    star_ori = fits.open(ori)

    #SNR

    SNR = star_ori[0].header["SNR"]

    nor_flux = star[4].data
    nor_ivar = star[5].data
    inf_flux = star[6].data

    VHELIO = star[11].data
    shift = star[1].data

    model.fitting_spectrum_parameters_single_add_photon_noise(normalized_flux=nor_flux, normalized_ivar=nor_ivar,
                                                              inf_flux=inf_flux, SNR=SNR)

    # attention! The result is in m/s
    un_SNR = model.uncertainty_SNR[2:]


    #Let's plot it:



    font = {'family': 'normal',
            'weight': 'bold',
            'size': 8}

    matplotlib.rc('font', **font)

    # only choose individual visits:

    lens = len(star[0].data[:,0])

    flux = star[4].data[2:lens,:]
    ivar = star[5].data[2:lens,:]
    inf_flux = star[6].data[2:lens,:]
    mix_flux = star[7].data[2:lens,:]

    parameters =star[0].data[2:lens,:]
    inf_labels = star[8].data[2:lens,:]

    #print("check shape")
    #print(flux.shape,ivar.shape,parameters.shape,parameters_sim.shape)

    # select visits with 2b<a+c

    a = parameters[:,0]
    b = parameters[:,1]
    c = parameters[:,2]

    mask = []

    for l in range(0,len(a)):
        if 2*b[l]>a[l]+c[l]:
            mask.append(1)

        else:
            mask.append(0)



    N = len(flux[:,0])
    name=path

    name = name.replace("/Volumes/Data_2TB/Data/dr13_all/", "")
    name = name.replace(".fits", "")

    # If you add combined spectra, change 0 to 2

    for i in range(0, N):

        # left
        plt.subplot(N, 2, 2 * i + 1)


        ## Add notifications for 2b<a+c
        if mask[i]==1:
            data_label = "Data flux"
            color = None
        else:
            #data_label = colored("Data flux of the visit with 2b < a+c","red")
            data_label = "Data flux"
            plt.plot(0,0,"ro",label="Data flux of the visit with $2b<a+c$")
            color = "red"


        plt.step(wl, flux[i,:], "k", label="Data flux", linewidth=0.7, alpha=0.8)
        plt.errorbar(wl,flux[i,:],yerr=ivar[i,:]**(-0.5),alpha = 0.2,ecolor='k')

        #plt.plot(wl, inf_flux[i,:], "g", label="Inferred flux from fitting separately", linewidth=0.7, alpha=0.4)


        # inferred flux
        plt.plot(wl, inf_flux[i,:], "b", label="Inferred flux  Teff=%.2f K logg=%.2f Fe/H =%.2f " % (
        inf_labels[i,0], inf_labels[i,1], inf_labels[i,2]), linewidth=0.7, alpha=0.4)


        # opt flux

        plt.plot(wl, mix_flux[i, :], "r", label="Mixed flux $RV_{shift}$=%.2f $Km/s$ Photon Limit = %.3f $Km/s$" % (
        (c[i]-a[i])/(a[i]+b[i]+c[i])*4.14468,un_SNR[i]/1000), linewidth=0.7, alpha=0.4)




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


        ## Add notifications for 2b<a+c
        if mask[i]==1:
            data_label = "Data flux"
            color = None
        else:
            #data_label = colored("Data flux of the visit with 2b < a+c","red")
            data_label = "Data flux"
            plt.plot(0,0,"ro",label="Data flux of the visit with $2b<a+c$")
            color = "red"




        plt.step(wl, flux[i,:], "k", label=data_label, linewidth=0.7, alpha=0.8)
        plt.errorbar(wl,flux[i,:],yerr=ivar[i,:]**(-0.5),alpha = 0.2,ecolor='k')

        #plt.plot(wl, inf_flux[i,:], "g", label="Inferred flux from fitting separately", linewidth=0.7, alpha=0.4)


        # inferred flux
        plt.plot(wl, inf_flux[i,:], "b", label="Inferred flux  Teff=%.2f K logg=%.2f Fe/H =%.2f" % (
        inf_labels[i,0], inf_labels[i,1], inf_labels[i,2]), linewidth=0.7, alpha=0.4)

        # opt flux

        plt.plot(wl, mix_flux[i, :], "r", label="Mixed flux $RV_{shift}$=%.2f $Km/s$ Photon Limit = %.3f $Km/s$" % (
        (c[i]-a[i])/(a[i]+b[i]+c[i])*4.14468,un_SNR[i]/1000), linewidth=0.7, alpha=0.4)




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

    #plt.show()

    fig = matplotlib.pyplot.gcf()

    # adjust the size based on the number of visit
    visit = len(flux[:,0])-2

    height = 2+2*visit

    fig.set_size_inches(18.5, height)

    save_path = "/Users/caojunzhi/Downloads/upload_20170501_David"+"/"+str(name)+".png"
    fig.savefig(save_path, dpi=500)

    plt.close()

#######################################
    # plot RV
    print("plot RV")

    path_fit = path



    font = {'family': 'normal',

            'size': 24}

    matplotlib.rc('font', **font)


    dat = Table.read(ori)
    image = fits.open(ori)

    VHELIO = dat["VHELIO"]


    VHELIO = VHELIO.ravel()



    # read HJDs

    # Read BJD for individual visits
    HJD = []
    VHELIO_un = []


    for k in range(0, len(VHELIO)):
        name_k = "HJD" + str(k + 1)
        name_k2 = "VERR"+ str(k + 1)
        #print(name_k)

        HJD_i = image[0].header[name_k]
        HJD.append(HJD_i)
        VHELIO_un.append(image[0].header[name_k2])

    HJD = np.array(HJD)
    VHELIO_un = np.array(VHELIO_un)



    # Your data from new
    image = fits.open(path_fit)

    """

    a =image[6].data[:,3]
    b = image[6].data[:,4]
    c = image[6].data[:,5]

    ve_sim = (c-a)/(a+b+c)*4.14468

    ve_sim_un = image[7].data


    """



    ve = image[1].data/1000

    # Only show visits:

    le = len(ve)

    ve = ve[2:le]

    a =image[0].data[2:le,0]
    b = image[0].data[2:le,1]
    c = image[0].data[2:le,2]

    mask = []
    mask_data = []
    for m in range(0,len(a)):
        if (2*b[m] > a[m]+c[m]):
            mask.append(0)
            mask_data.append(1)

        else:
            mask.append(1)
            mask_data.append(0)

    mask = np.array(mask,dtype=bool)
    mask_data = np.array(mask_data,dtype=bool)


    plt.subplot(1,1,1)

    std = np.std(VHELIO)
    std_sim = np.std(ve+VHELIO)

    plt.plot(HJD[mask_data],VHELIO[mask_data],"ko",markersize=5,label="DR13 RVs scatter = %.3f $km/s$ Photon limit = %.3f $Km/s$"%(std,np.nanmean(un_SNR)/1000))
    plt.plot(HJD[mask], VHELIO[mask], "ro", markersize=5, label="Visits with $2b<a+c$")




    plt.plot(HJD[mask_data], (VHELIO+ve)[mask_data], "bo", markersize=5, label="After correction scatter = %.3f $km/s$ Photon limit = %.3f $Km/s$"%(std_sim,np.nanmean(un_SNR)/1000))

    plt.legend(loc='upper center')

    plt.plot(HJD[mask], (VHELIO + ve)[mask], "ro", markersize=5, label="Visits with $2b<a+c$")



    plt.errorbar(HJD[mask_data], VHELIO[mask_data], yerr=VHELIO_un[mask_data], fmt='ko',alpha = 0.5)
    plt.errorbar(HJD[mask_data], (VHELIO + ve)[mask_data], yerr=VHELIO_un[mask_data], fmt='bo',alpha = 0.5)

    plt.errorbar(HJD[mask], VHELIO[mask], yerr=VHELIO_un[mask], fmt='ro',alpha = 0.5)
    plt.errorbar(HJD[mask], (VHELIO + ve)[mask], yerr=VHELIO_un[mask], fmt='ro',alpha = 0.5)



    name= ori

    name = name.replace("/Volumes/Data_2TB/Data/APOGEE_DR13_Apstar/","")
    name = name.replace(".fits","")

    plt.xlabel("BJD")
    plt.ylabel("Radial velocity $Km/s$")
    plt.suptitle("%s"%name)



    # plt.show()
    # save them:

    fig = matplotlib.pyplot.gcf()

    # adjust the size based on the number of visit

    fig.set_size_inches(14.5, 8.5)

    save_path = "/Users/caojunzhi/Downloads/upload_20170501_David/rv/"+"velocity_"+name+".png"
    fig.savefig(save_path, dpi=500)

    plt.close()





for i in range(10250,10260):
    plot_star(path=path_fits[i], ori=path_ori[i])

