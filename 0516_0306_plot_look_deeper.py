from scipy import stats
import numpy as np

import matplotlib.pyplot as plt
import matplotlib
import pickle
import scipy

from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time

# read data:



pkl_file = open(
    '/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/dr13_nick_fits.pkl',
    'rb')
nick_path = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(
    '/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/dr13_mine_fits.pkl',
    'rb')
jason_path = pickle.load(pkl_file)
pkl_file.close()


def MJD2BJD(mjd, target_coord, site_location=(-105.820417, 32.780361)):
    """do the conversion
mjd -- input mjd, scale is utc
target_coord -- the coordinate of the target, in astropy.coord format,
to caculate the light travel time
site_location -- location of the telescope, to make accurate calcualtion of
light travel time and tdb conversion. Not very important here. The default value
is for Greenwich Observatory.
"""
    t = Time(mjd, format='mjd', scale='utc', location=site_location)
    # calculate light travel time
    ltt = t.light_travel_time(target_coord)
    # print(t, ltt)
    # convert t to tdb, and add light travel time
    t_out = (t.tdb + ltt).jd

    return t_out


"""
c = SkyCoord(RA,DEC, frame='icrs', unit='deg')

BJD = MJD2BJD(MJD,c)

"""


# plot:

def plot_jason_aspcap(jason):
    #######################################
    # plot RV
    print("plot RV")

    star_jason = fits.open(jason)


    font = {'family': 'normal',

            'size': 24}

    matplotlib.rc('font', **font)

    VHELIO = star_jason[11].data
    # read BJDs

    BJD = star_jason[13].data

    ve = star_jason[1].data / 1000

    # Only show visits:

    ve = ve[2:]

    a = star_jason[0].data[2:, 0]
    b = star_jason[0].data[2:, 1]
    c = star_jason[0].data[2:, 2]

    mask = []
    mask_data = []
    for m in range(0, len(a)):
        if (2 * b[m] > a[m] + c[m]):
            mask.append(0)
            mask_data.append(1)

        else:
            mask.append(1)
            mask_data.append(0)

    mask = np.array(mask, dtype=bool)
    mask_data = np.array(mask_data, dtype=bool)

    plt.subplot(1, 1, 1)

    std = np.std(VHELIO)
    std_sim = np.std(ve + VHELIO)

    plt.plot(BJD, VHELIO, "ko", markersize=5,
             label="DR13 RVs scatter = %.3f $km/s$" % (std))

    plt.plot(BJD[mask_data], (VHELIO + ve)[mask_data], "bo", markersize=5,
             label="After correction scatter = %.3f $km/s$" % (std_sim))

    plt.plot(BJD[mask], (VHELIO + ve)[mask], "ro", markersize=5, label="Visits with $2b<a+c$")


    plt.legend(loc='upper center')


    name = jason

    name = name.replace("/Volumes/Data_2TB/Data/dr13_all/", "")
    name = name.replace(" .fits", "")

    plt.xlabel("BJD")
    plt.ylabel("Radial velocity $Km/s$")
    plt.suptitle("%s" % name)

    # plt.show()
    # save them:

    fig = matplotlib.pyplot.gcf()

    # adjust the size based on the number of visit

    fig.set_size_inches(16.5, 8.5)

    # select some plots:

    if (std<0.3)&(std_sim>0.6):
        save_path = "/Users/caojunzhi/Downloads/upload_0516_David/" + "velocity_" + name + ".png"
        fig.savefig(save_path, dpi=500)
        print("yes")
        plt.close()
        return 1
    else:
        plt.close()
        print("Nope")
        return 0

    plt.close()



def plot_jason_nick(jason,nick):

    #######################################
    # plot RV
    print("plot RV")

    star_jason = fits.open(jason)


    font = {'family': 'normal',

            'size': 24}

    matplotlib.rc('font', **font)

    VHELIO = star_jason[11].data
    # read BJDs

    BJD = star_jason[13].data

    ve = star_jason[1].data / 1000

    # Only show visits:

    ve = ve[2:]

    a = star_jason[0].data[2:, 0]
    b = star_jason[0].data[2:, 1]
    c = star_jason[0].data[2:, 2]

    mask = []
    mask_data = []
    for m in range(0, len(a)):
        if (2 * b[m] > a[m] + c[m]):
            mask.append(0)
            mask_data.append(1)

        else:
            mask.append(1)
            mask_data.append(0)

    mask = np.array(mask, dtype=bool)
    mask_data = np.array(mask_data, dtype=bool)

    # read nick's data

    table = Table.read(nick)

    VHELIO_nick = table["VHELIO"]

    MJD = table["MJD"]
    RA = table["RA"][0]
    DEC = table["DEC"][0]
    c = SkyCoord(RA, DEC, frame='icrs', unit='deg')

    ni = len(MJD)

    BJD_nick = []
    for j in range(0, ni):
        BJD_nick.append(MJD2BJD(MJD[j], c))



    plt.subplot(1, 1, 1)

    std = np.std(VHELIO)
    std_sim = np.std(ve + VHELIO)
    std_nick = np.std(VHELIO_nick)

    # plot nick's


    plt.plot(BJD_nick,VHELIO_nick, "ko", markersize=5,
             label="Nick's RVs scatter = %.3f $km/s$" % (std_nick))



    plt.plot(BJD[mask_data], (VHELIO + ve)[mask_data], "bo", markersize=5,
             label="After correction scatter = %.3f $km/s$" % (std_sim))

    plt.plot(BJD[mask], (VHELIO + ve)[mask], "ro", markersize=5, label="Visits with $2b<a+c$")


    plt.legend(loc='upper center')


    name = jason

    name = name.replace("/Volumes/Data_2TB/Data/dr13_all/", "")
    name = name.replace(" .fits", "")

    plt.xlabel("BJD")
    plt.ylabel("Radial velocity $Km/s$")
    plt.suptitle("%s" % name)

    # plt.show()
    # save them:

    fig = matplotlib.pyplot.gcf()

    # adjust the size based on the number of visit

    fig.set_size_inches(16.5, 8.5)

    # select some plots:

    if (std_nick<0.3)&(std_sim>0.6):
        save_path = "/Users/caojunzhi/Downloads/upload_0516_David/" + "velocity_jason_nick_" + name + ".png"
        fig.savefig(save_path, dpi=500)
        print("yes")
        plt.close()
        return 1
    else:
        plt.close()
        print("Nope")
        return 0


    plt.close()


pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()


def plot_flux(jason):
    # set font size

    font = {'family': 'normal',
            'weight': 'bold',
            'size': 8}

    matplotlib.rc('font', **font)

    # only choose individual visits:
    star = fits.open(jason)

    flux = star[4].data
    ivar = star[5].data
    inf_flux = star[6].data
    parameters = star[0].data
    inf_labels = star[8].data
    chi_inf = star[9].data

    max_chi = chi_inf.argsort()[-1]

    name = jason

    N = len(flux[:, 0])
    for i in range(0, N):

        # left
        plt.subplot(N, 2, 2 * i + 1)

        name = name.replace("/Volumes/Data_2TB/Data/dr13_all/", "")
        name = name.replace(".fits", "")

        ## Add inf labels

        chi_i = chi_inf[i]

        if i==max_chi:
            plt.plot([],[],"r",label="Epoch with the biggest chi-squared")

        plt.step(wl, flux[i, :], "k", label="Data flux, chi-squared = %.2f"%(chi_i), linewidth=0.7, alpha=0.8)
        plt.errorbar(wl,ivar[i,:]**(-0.5),alpha = 0.1)

        plt.plot(wl, inf_flux[i, :], "b",
                 label="Inferred flux Teff=%.2f K logg=%.2f Fe/H =%.2f" % (
                     inf_labels[i, 0], inf_labels[i, 1], inf_labels[i, 2]), linewidth=0.7, alpha=0.4)

        plt.xlabel("Wave length $\AA$")
        plt.ylabel("Flux", fontsize=20)
        axes = plt.gca()
        axes.set_xlim([15660, 15780])

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


        if i==max_chi:
            plt.plot([],[],"r",label="Epoch with the biggest chi-squared")


        plt.step(wl, flux[i, :], "k", linewidth=0.7, alpha=0.8)
        # plt.errorbar(wl,ivar[i,:]**(-0.5),alpha = 0.1)

        ai = parameters[i, 0]
        bi = parameters[i, 1]
        ci = parameters[i, 2]

        plt.plot(wl, inf_flux[i, :], "b",
                 label="Velocity shift a=%.2f b=%.2f c =%.2f $RV_{shift}$=%.2f $m/s$" % (
                     ai, bi, ci, (ci - ai) / (ai + bi + ci) * 4144.68), linewidth=0.7, alpha=0.4)

        plt.xlabel("Wave length $\AA$")
        plt.ylabel("Flux", fontsize=20)
        axes = plt.gca()
        axes.set_xlim([16160, 16280])
        # share x axis

        if i == N - 1:
            non = 1
        else:
            axes.set_xticklabels([])

        # axes.set_xlim([16160,16280])
        axes.set_ylim([0.5, 2])
        axes.set_yticks(np.arange(0.5, 1.99, 0.5))
        plt.legend()

    plt.suptitle("The fitting result of star %s" % name, fontsize=20)

    # share x
    plt.subplots_adjust(hspace=.0)


    # plt.show()
    # save them:

    fig = matplotlib.pyplot.gcf()

    # adjust the size based on the number of visit

    height = 4+2*N

    fig.set_size_inches(16.5, height)

    # select some plots:


    save_path = "/Users/caojunzhi/Downloads/upload_0516_David/flux/" + "flux_" + name + ".png"
    fig.savefig(save_path, dpi=500)

    plt.close()





for i in range(0,len(nick_path)):
    print("doing star %d"%i)
    m1 = plot_jason_aspcap(jason=jason_path[i])
    m2 = plot_jason_nick(jason=jason_path[i],nick=nick_path[i])

    if (m1==1) | (m2==1):
        plot_flux(jason_path[i])
