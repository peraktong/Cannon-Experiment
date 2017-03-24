import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
import pickle
from os.path import isfile, join
from os import listdir


from astropy.time import Time





pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()



#Now you only need to plot old method

def plot_visit_spectra(name):
    # set font size

    font = {'family': 'normal',
            'weight': 'bold',
            'size': 8}

    matplotlib.rc('font', **font)

    # only choose individual visits:
    star = fits.open(name)

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

    # If you add combined spectra, change 0 to 2

    for i in range(0, N):

        # left
        plt.subplot(N, 2, 2 * i + 1)

        name = name.replace("/Users/caojunzhi/Desktop/Data/d_15/", "")
        name = name.replace(".fits", "")

        chi_old = (flux[i,:]-inf_flux[i,:])**2*ivar[i,:]
        chi_old = np.sum(chi_old)

        chi_new = (flux[i, :] - mix_flux[i, :]) ** 2 * ivar[i, :]
        chi_new = np.sum(chi_new)




        ## Add notifications for 2b<a+c
        if mask[i]==1:
            data_label = "Data flux"
            color = None
        else:
            #data_label = colored("Data flux of the visit with 2b < a+c","red")
            data_label = "Data flux"
            plt.plot(0,0,"ro",label="Data flux of the visit with 2b < a+c")
            color = "red"


        plt.step(wl, flux[i,:], "k", label="Data flux", linewidth=0.7, alpha=0.8)
        plt.errorbar(wl,flux[i,:],yerr=ivar[i,:]**(-0.5),alpha = 0.2,ecolor='k')

        #plt.plot(wl, inf_flux[i,:], "g", label="Inferred flux from fitting separately", linewidth=0.7, alpha=0.4)


        # inferred flux
        plt.plot(wl, inf_flux[i,:], "b", label="Inferred flux  Teff=%.2f K logg=%.2f Fe/H =%.2f chi-squared = %.2f" % (
        inf_labels[i,0], inf_labels[i,1], inf_labels[i,2],chi_old), linewidth=0.7, alpha=0.4)

        # opt flux

        plt.plot(wl, mix_flux[i, :], "r", label="Mixed flux  Velocity shift a=%.2f b=%.2f c =%.2f $RV_{shift}$=%.2f $m/s$ chi-squared = %.2f" % (
        a[i],b[i],c[i],(c[i]-a[i])/(a[i]+b[i]+c[i])*4144.68,chi_new), linewidth=0.7, alpha=0.4)










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
            plt.plot(0,0,"ro",label="Data flux of the visit with 2b < a+c")
            color = "red"




        plt.step(wl, flux[i,:], "k", label=data_label, linewidth=0.7, alpha=0.8)
        plt.errorbar(wl,flux[i,:],yerr=ivar[i,:]**(-0.5),alpha = 0.2,ecolor='k')

        #plt.plot(wl, inf_flux[i,:], "g", label="Inferred flux from fitting separately", linewidth=0.7, alpha=0.4)


        # inferred flux
        plt.plot(wl, inf_flux[i,:], "b", label="Inferred flux  Teff=%.2f K logg=%.2f Fe/H =%.2f chi-squared = %.2f" % (
        inf_labels[i,0], inf_labels[i,1], inf_labels[i,2],chi_old), linewidth=0.7, alpha=0.4)

        # opt flux

        plt.plot(wl, mix_flux[i, :], "r", label="Mixed flux  Velocity shift a=%.2f b=%.2f c =%.2f $RV_{shift}$=%.2f $m/s$ chi-squared = %.2f" % (
        a[i],b[i],c[i],(c[i]-a[i])/(a[i]+b[i]+c[i])*4144.68,chi_new), linewidth=0.7, alpha=0.4)




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

    name = name.replace("/Users/caojunzhi/Desktop/Data/dr13_red_clump/","")

    save_path = "/Users/caojunzhi/Downloads/upload_20170322/flux"+"/"+str(name)+".png"
    fig.savefig(save_path, dpi=500)

    plt.close()


# input apstar
def Plot_RV_before_after(path_fit,ori):


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

    plt.plot(HJD[mask_data],VHELIO[mask_data],"ko",markersize=5,label="DR13 RVs and the standard deviation = %.3f $km/s$"%std)
    plt.plot(HJD[mask], VHELIO[mask], "ro", markersize=5, label="Visits with 2b<a+c")




    plt.plot(HJD[mask_data], (VHELIO+ve)[mask_data], "bo", markersize=5, label="After shift RVs and the standard deviation = %.3f $km/s$"%std_sim)

    plt.legend(loc='upper center')

    plt.plot(HJD[mask], (VHELIO + ve)[mask], "ro", markersize=5, label="Visits with 2b<a+c")



    plt.errorbar(HJD[mask_data], VHELIO[mask_data], yerr=VHELIO_un[mask_data], fmt='ko',alpha = 0.5)
    plt.errorbar(HJD[mask_data], (VHELIO + ve)[mask_data], yerr=VHELIO_un[mask_data], fmt='bo',alpha = 0.5)

    plt.errorbar(HJD[mask], VHELIO[mask], yerr=VHELIO_un[mask], fmt='ro',alpha = 0.5)
    plt.errorbar(HJD[mask], (VHELIO + ve)[mask], yerr=VHELIO_un[mask], fmt='ro',alpha = 0.5)



    name= ori

    name = name.split("-")[2]
    name = name.replace(".fits","")

    plt.xlabel("HJD")
    plt.ylabel("Radial velocity $Km/s$")
    plt.suptitle("%s"%name)



    # plt.show()
    # save them:

    fig = matplotlib.pyplot.gcf()

    # adjust the size based on the number of visit

    fig.set_size_inches(14.5, 8.5)

    save_path = "/Users/caojunzhi/Downloads/upload_20170322/rv/"+"velocity_"+name+".png"
    fig.savefig(save_path, dpi=500)

    plt.close()

















    #print(HJD)




## Plot fluxes




pkl_file = open('dr13_rc_fits_atleast_4_epoch.pkl', 'rb')
my_path = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('dr13_rc_ori_atleast_4_epoch.pkl', 'rb')
path_ori = pickle.load(pkl_file)
pkl_file.close()


# randomly choose 20 stars

index = []
a = 0
number = 20

for r in range(0,number):
    a = a+np.random.randint(1,len(my_path)//number-1)
    index.append(a)
index = np.array(index)

my_path = my_path[index]
path_ori = path_ori[index]



# plot flux

for i in range(0,len(my_path)):

    try:
        print("Doing star %d" % (i + 1))
        name_i = my_path[i]

        plot_visit_spectra(name_i)


    except IndexError:
        print("This one fails")



    except OSError:
        print("This one fails")




## Plot RVs



for i in range(0,len(path_ori)):


    try:

        print("doing star %d"%i)

        name_i = my_path[i]
        ori_i = path_ori[i]

        Plot_RV_before_after(name_i, ori_i)




    except IndexError:
        print("This one fails")



    except OSError:
        print("This one fails")
        print("os")
        print(name_i, ori_i)