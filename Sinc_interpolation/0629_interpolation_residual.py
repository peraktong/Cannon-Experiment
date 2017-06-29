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
    if x>0:
        return math.log10(x)
    else:
        return -np.inf





pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()



data_path = "/Users/caojunzhi/Desktop/Data_example/"

# s = wl

apstar = fits.open(data_path+"apStar-r6-2M00005143+5615568.fits")

apstar_table = Table.read(data_path+"apStar-r6-2M00005143+5615568.fits")

jd_array = np.array(apstar_table[0]["JD"])
print(jd_array)

## Let's save everything into a fits file.


# import data
image_path = np.array(["apVisit-r6-5094-55874-088.fits","apVisit-r6-5094-56643-088.fits","apVisit-r6-5094-56651-082.fits"])



# index!!
index=0

image = fits.open(data_path+image_path[index],ignore_missing_end=True)

dat = Table.read(data_path+image_path[index])


print(dat[0]["JD"])


# flux at 1 and wl at 4
flux_raw =image[1].data.ravel()[::-1]

# Three chips

XSHIFT = image[0].header["XSHIFT"]

# Dither
# red green blue 5, 4.25 and 3.5
beta_red = (XSHIFT+5)*4.144/(3*10**5)
beta_green = (XSHIFT+4.25)*4.144/(3*10**5)
beta_blue = (XSHIFT+3.5)*4.144/(3*10**5)


red = image[4].data[0,:]
green = image[4].data[1,:]
blue = image[4].data[2,:]

# Dither:

red = (-beta_red+1)*red

green = (-beta_green+1)*green

blue = (-beta_blue+1)*blue

wl_raw =np.append(red,[green,blue])
wl_raw = wl_raw[::-1]


## Do interpolation:

# Add velocity

# 3.5 dither for red!!!


y_inter = np.interp(x=wl,xp=wl_raw,fp=flux_raw)

print(y_inter.shape)


# labels
array = np.array([[  4.80458566e+03 ,  2.57179854e+00 , -2.01372206e-01],
 [  4.79263171e+03  , 2.64864274e+00,  -1.82310447e-01],
 [  4.78393682e+03 ,  2.65073084e+00 , -1.81901661e-01],
 [  4.80550440e+03  , 2.57941485e+00 , -2.04721940e-01],
 [  4.79156434e+03 ,  2.63379621e+00 , -1.83733041e-01],
 [  4.77938538e+03 ,  2.65496682e+00 , -1.83967319e-01]])




font = {'family': 'normal',
        'weight': 'bold',
        'size': 15}


matplotlib.rc('font', **font)


plt.subplot(2,1,1)

plt.step(wl,apstar[1].data[2+index,:],"k", label = "$APOGEE\quad team\quad spectra$",linewidth=0.7, alpha=1)
plt.plot(wl,20*apstar[2].data[2+index,:],"r", label = "$20\quad times\quad spectra\quad error$",linewidth=0.7, alpha=0.5)
plt.plot(wl,20*(apstar[1].data[2+index,:]-y_inter),"g",label = "$20\quad times\quad Residual\quad$", linewidth=0.7, alpha=0.5)

plt.xlabel("$Wave\quad length\quad \AA$", fontsize=20)
plt.ylabel("$Flux$", fontsize=20)
plt.suptitle("$Comparison\quad of\quad the\quad spectra\quad for\quad one\quad epoch\quad %s$"%(str(index+1)), fontsize=30)
plt.legend()

axes = plt.gca()
axes.set_xlim([15660, 15780])
axes.set_ylim([-500,1000])

plt.subplot(2,1,2)

plt.step(wl,apstar[1].data[2+index,:],"k", label = "$APOGEE\quad team\quad spectra$",linewidth=0.7, alpha=1)
plt.plot(wl,20*apstar[2].data[2+index,:],"r", label = "$20\quad times\quad spectra\quad error$",linewidth=0.7, alpha=0.5)
plt.plot(wl,20*(apstar[1].data[2+index,:]-y_inter),"g",label = "$20\quad times\quad Residual\quad$", linewidth=0.7, alpha=0.5)

plt.ylabel("$Flux$", fontsize=20)
plt.suptitle("$Comparison\quad of\quad the\quad spectra\quad for\quad one\quad epoch\quad %s$"%(str(index+1)), fontsize=30)
plt.legend()

axes = plt.gca()
axes.set_xlim([16160, 16280])
axes.set_ylim([-500,1000])


# save it:

fig = matplotlib.pyplot.gcf()

fig.set_size_inches(14.5, 11.5)

save_path = "/Users/caojunzhi/Downloads/upload_20170621_David/" + "New_flux_residual_"+ str(index)+ ".png"
fig.savefig(save_path, dpi=500)

plt.close()

#### Compare spectra:



font = {'family': 'normal',
        'weight': 'bold',
        'size': 15}


matplotlib.rc('font', **font)


plt.subplot(2,1,1)

plt.step(wl,apstar[1].data[2+index,:],"k", label = "$APOGEE\quad team\quad Teff=%.2fK \quad Logg=%.2f dex \quad FeH=%.2f dex$"%(array[index,0],array[index,1],array[index,2]),linewidth=0.7, alpha=1)
plt.plot(wl,y_inter,"g",label = "$My\quad code\quad Teff=%.2f K \quad Logg=%.2f dex\quad FeH=%.2f dex$"%(array[index+3,0],array[index+3,1],array[index+3,2]), linewidth=0.7, alpha=0.5)

plt.xlabel("$Wave\quad length\quad \AA$", fontsize=20)
plt.ylabel("$Flux$", fontsize=20)
plt.suptitle("$Comparison\quad of\quad the\quad spectra\quad for\quad one\quad epoch\quad %s$"%(str(index+1)), fontsize=30)
plt.legend()

axes = plt.gca()
axes.set_xlim([15660, 15780])

plt.subplot(2,1,2)


plt.step(wl,apstar[1].data[2+index,:],"k", label = "$APOGEE\quad team\quad Teff=%.2fK \quad Logg=%.2f dex \quad FeH=%.2f dex$"%(array[index,0],array[index,1],array[index,2]),linewidth=0.7, alpha=1)
plt.plot(wl,y_inter,"g",label = "$My\quad code\quad Teff=%.2f K \quad Logg=%.2f dex\quad FeH=%.2f dex$"%(array[index+3,0],array[index+3,1],array[index+3,2]), linewidth=0.7, alpha=0.5)

plt.ylabel("$Flux$", fontsize=20)
plt.suptitle("$Comparison\quad of\quad the\quad spectra\quad for\quad one\quad epoch\quad %s$"%(str(index+1)), fontsize=30)
plt.legend()

axes = plt.gca()
axes.set_xlim([16160, 16280])


# save it:

fig = matplotlib.pyplot.gcf()

fig.set_size_inches(14.5, 11.5)

save_path = "/Users/caojunzhi/Downloads/upload_20170621_David/" + "New_flux_"+ str(index)+ ".png"
fig.savefig(save_path, dpi=500)

plt.close()













