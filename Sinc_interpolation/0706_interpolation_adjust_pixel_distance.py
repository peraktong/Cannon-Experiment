import matplotlib.pyplot as plt
import matplotlib
import pickle
import math
import numpy as np
import os
from astropy.io import fits
import pickle
from astropy.table import Table
import AnniesLasso_2 as tc
import sincinterpol


def log(x):
    if x > 0:
        return math.log10(x)
    else:
        return -np.inf


log = np.vectorize(log)

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()

wl_log = log(wl)

# value of a
a_value = "a=2"

# Both a=2 and 3 are good.

a = 2.5

a_value = "a=2".replace("2", str(a))


def sinc_interp(x, s, u, a):
    # Your x is the raw flux
    # Your s is the log wave length of the raw flux. s should have equal distance

    # Your u is the log wl, which has equal distance between the neighborhoods.
    # The length of u is 8575 and you can use log wl

    if len(x) != len(s):
        print("len(x) should be equal to len(s")

    # I don't think these two methods have a big different.

    # Can we use this?
    # parameter a:

    # T is very small

    N = len(s)

    # T is a matrix with the same shape with sincM

    #s2 = np.append(s,s[-1])
    #s1 = np.append(0,s)


    #delta_s = (s2-s1)[0:-1]

    T = (s[N - 1] - s[0]) / N

    # T = np.tile(delta_s[:, np.newaxis], (1, len(u)))

    sincM = np.tile(u, (len(s), 1)) - np.tile(s[:, np.newaxis], (1, len(u)))

    # a=2 or 3

    mask = abs(sincM / T) > a

    windows = np.sinc(sincM / T) * np.sinc(sincM / T / a)
    windows[mask] = 0

    y = np.dot(x, windows)

    return y


log = np.vectorize(log)

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()

data_path = "/Users/caojunzhi/Desktop/Data_example/"

# s = wl

apstar = fits.open(data_path + "apStar-r6-2M00005143+5615568.fits")

apstar_table = Table.read(data_path + "apStar-r6-2M00005143+5615568.fits")

jd_array = np.array(apstar_table[0]["JD"])
print(jd_array)

## Let's save everything into a fits file.


# import data
image_path = np.array(
    ["apVisit-r6-5094-55874-088.fits", "apVisit-r6-5094-56643-088.fits", "apVisit-r6-5094-56651-082.fits"])

# index!! For which visit

for index in range(0, 3):


    print("doing a=%.2f index=%d" % (a, index))

    image = fits.open(data_path + image_path[index], ignore_missing_end=True)

    dat = Table.read(data_path + image_path[index])

    print(dat[0]["JD"])

    # flux at 1 and wl at 4
    flux_raw = image[1].data.ravel()[::-1]

    # Three chips

    XSHIFT = image[0].header["XSHIFT"]

    # Dither
    # red green blue 5, 4.25 and 3.5

    """

    beta_red = (XSHIFT+5)*4.144/(3*10**5)
    beta_green = (XSHIFT+4.25)*4.144/(3*10**5)
    beta_blue = (XSHIFT+3.5)*4.144/(3*10**5)



    """

    beta_red = (XSHIFT + 3.5) * 4.144 / (3 * 10 ** 5)
    beta_green = (XSHIFT + 3.5) * 4.144 / (3 * 10 ** 5)
    beta_blue = (XSHIFT + 3.5) * 4.144 / (3 * 10 ** 5)

    red = image[4].data[0, :]
    green = image[4].data[1, :]
    blue = image[4].data[2, :]

    # Dither:

    red = (-beta_red + 1) * red

    green = (-beta_green + 1) * green

    blue = (-beta_blue + 1) * blue

    wl_raw = np.append(red, [green, blue])
    wl_raw = wl_raw[::-1]

    # Let's only focus on 15660 to 15780

    flux = apstar[1].data[2 + index, :]
    # Choose part:

    
    
    mask = (wl_raw > 15660) & (wl_raw < 15780)
    mask_wl = (wl > 15660) & (wl < 15780)

    



    wl_raw_part = wl_raw[mask]
    flux_raw_part = flux_raw[mask]
    wl_log = log(wl)
    wl_raw_log = log(wl_raw)

    print("doing interpolation")
    y_inter_p1 = sinc_interp(x=flux_raw, s=wl_raw_log, u=wl_log, a=a)

    # normalize
    y_inter_p1 = y_inter_p1 * np.nanmean(flux[mask_wl]) / np.nanmean(y_inter_p1[mask_wl])


    
    mask = (wl_raw > 16160) & (wl_raw < 16280)
    mask_wl = (wl > 16160) & (wl < 16280)



    wl_raw_part = wl_raw[mask]
    flux_raw_part = flux_raw[mask]
    wl_log = log(wl)
    wl_raw_log = log(wl_raw)

    print("doing interpolation")
    y_inter_p2 = sinc_interp(x=flux_raw, s=wl_raw_log, u=wl_log, a=a)

    # normalize
    y_inter_p2 = y_inter_p2 * np.nanmean(flux[mask_wl]) / np.nanmean(y_inter_p2[mask_wl])

    ## Do interpolation:


    # Add velocity

    # labels
    array = np.array([[4.80458566e+03, 2.57179854e+00, -2.01372206e-01],
                      [4.79263171e+03, 2.64864274e+00, -1.82310447e-01],
                      [4.78393682e+03, 2.65073084e+00, -1.81901661e-01],
                      [4.80550440e+03, 2.57941485e+00, -2.04721940e-01],
                      [4.79156434e+03, 2.63379621e+00, -1.83733041e-01],
                      [4.77938538e+03, 2.65496682e+00, -1.83967319e-01]])

    ### Normalize y spectra:

    # y_inter = y_inter*(np.nansum(apstar[1].data[2+index,:]))/np.nansum(y_inter)


    ############ plots

    font = {'family': 'normal',
            'weight': 'bold',
            'size': 15}

    matplotlib.rc('font', **font)

    plt.subplot(2, 1, 1)

    plt.step(wl, apstar[1].data[2 + index, :], "k", label="$APOGEE\quad team\quad spectra$", linewidth=0.7, alpha=1)
    plt.plot(wl, 20 * apstar[2].data[2 + index, :], "r", label="$20\quad times\quad spectra\quad error$",
             linewidth=0.7,
             alpha=0.5)
    plt.plot(wl, -20 * apstar[2].data[2 + index, :], "r", linewidth=0.7,
             alpha=0.5)

    plt.plot(wl, 20 * (apstar[1].data[2 + index, :] - y_inter_p1), "g", label="$20\quad times\quad Residual\quad$",
             linewidth=0.7, alpha=0.5)

    plt.xlabel("$Wave\quad length\quad \AA$", fontsize=20)
    plt.ylabel("$Flux$", fontsize=20)
    plt.suptitle(
        "$Comparison\quad of\quad the\quad spectra\quad for\quad one\quad epoch\quad %s$" % (str(index + 1)),
        fontsize=30)
    plt.legend()

    axes = plt.gca()
    axes.set_xlim([15660, 15780])
    axes.set_ylim([-500, 1000])

    plt.subplot(2, 1, 2)

    plt.step(wl, apstar[1].data[2 + index, :], "k", label="$APOGEE\quad team\quad spectra$", linewidth=0.7, alpha=1)
    plt.plot(wl, 20 * apstar[2].data[2 + index, :], "r", label="$20\quad times\quad spectra\quad error$",
             linewidth=0.7,
             alpha=0.5)
    plt.plot(wl, -20 * apstar[2].data[2 + index, :], "r", linewidth=0.7,
             alpha=0.5)

    plt.plot(wl, 20 * (apstar[1].data[2 + index, :] - y_inter_p2), "g", label="$20\quad times\quad Residual\quad$",
             linewidth=0.7, alpha=0.5)

    plt.ylabel("$Flux$", fontsize=20)
    plt.suptitle(
        "$Comparison\quad of\quad the\quad spectra\quad for\quad one\quad epoch\quad %s$" % (str(index + 1)),
        fontsize=30)
    plt.legend()

    axes = plt.gca()
    axes.set_xlim([16160, 16280])
    axes.set_ylim([-500, 1000])

    # save it:

    fig = matplotlib.pyplot.gcf()

    fig.set_size_inches(14.5, 11.5)

    save_path = "/Users/caojunzhi/Downloads/upload_201707_interpolation_David/" + "New_flux_residual_" + a_value + "_index=" + str(
        index) + ".png"
    fig.savefig(save_path, dpi=500)

    plt.close()

    #### Compare spectra:



    font = {'family': 'normal',
            'weight': 'bold',
            'size': 15}

    matplotlib.rc('font', **font)

    plt.subplot(2, 1, 1)

    plt.step(wl, apstar[1].data[2 + index, :], "k",
             label="$APOGEE\quad team\quad Teff=%.2fK \quad Logg=%.2f dex \quad FeH=%.2f dex$" % (
                 array[index, 0], array[index, 1], array[index, 2]), linewidth=0.7, alpha=1)
    plt.plot(wl, y_inter_p1, "g", label="$My\quad code\quad Teff=%.2f K \quad Logg=%.2f dex\quad FeH=%.2f dex$" % (
        array[index + 3, 0], array[index + 3, 1], array[index + 3, 2]), linewidth=0.7, alpha=0.5)

    plt.xlabel("$Wave\quad length\quad \AA$", fontsize=20)
    plt.ylabel("$Flux$", fontsize=20)
    plt.suptitle(
        "$Comparison\quad of\quad the\quad spectra\quad for\quad one\quad epoch\quad %s$" % (str(index + 1)),
        fontsize=30)
    plt.legend()

    axes = plt.gca()
    axes.set_xlim([15660, 15780])

    plt.subplot(2, 1, 2)

    plt.step(wl, apstar[1].data[2 + index, :], "k",
             label="$APOGEE\quad team\quad Teff=%.2fK \quad Logg=%.2f dex \quad FeH=%.2f dex$" % (
                 array[index, 0], array[index, 1], array[index, 2]), linewidth=0.7, alpha=1)
    plt.plot(wl, y_inter_p2, "g", label="$My\quad code\quad Teff=%.2f K \quad Logg=%.2f dex\quad FeH=%.2f dex$" % (
        array[index + 3, 0], array[index + 3, 1], array[index + 3, 2]), linewidth=0.7, alpha=0.5)

    plt.ylabel("$Flux$", fontsize=20)
    plt.suptitle(
        "$Comparison\quad of\quad the\quad spectra\quad for\quad one\quad epoch\quad %s$" % (str(index + 1)),
        fontsize=30)
    plt.legend()

    axes = plt.gca()
    axes.set_xlim([16160, 16280])

    # save it:

    fig = matplotlib.pyplot.gcf()

    fig.set_size_inches(14.5, 11.5)

    save_path = "/Users/caojunzhi/Downloads/upload_201707_interpolation_David/" + "New_flux_" + a_value + "_index=" + str(
        index) + ".png"
    fig.savefig(save_path, dpi=500)

    plt.close()














