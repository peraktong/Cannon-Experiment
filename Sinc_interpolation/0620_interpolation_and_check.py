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

## Let's save everything into a fits file.


# import data
image_path = np.array(["apVisit-r6-5094-55874-088.fits","apVisit-r6-5094-56643-088.fits","apVisit-r6-5094-56651-082.fits"])
image = fits.open(data_path+image_path[0],ignore_missing_end=True)

dat = Table.read(data_path+image_path[0])


i=0
print(dat[0]["JD"])

# flux at 1 and wl at 4
flux_raw =image[1].data.ravel()[::-1]
wl_raw = image[4].data.ravel()[::-1]

# Let's interpolate:
print("doing interpolation")
# https://gist.github.com/peraktong/ce566e47b069a05bfe8b85d2a97da71f

raw = np.arange(len(wl_raw))
target = np.arange(len(wl))

# x: raw data xp:target xt:raw

new_flux = sincinterpol.interp3 (x=flux_raw,xt=raw,xp=target,workers=8)

plt.step(wl,apstar[1].data[2+0,:],"k", linewidth=0.7, alpha=0.7)
plt.plot(wl,new_flux,"g", linewidth=0.7, alpha=0.2)
plt.xlabel("Wave length $\AA$")
plt.ylabel("Flux", fontsize=20)

axes = plt.gca()
axes.set_xlim([15660, 15780])

fig = matplotlib.pyplot.gcf()

# save it:

fig.set_size_inches(14.5, 11.5)

save_path = "/Users/caojunzhi/Downloads/upload_20170621_David/" + "New_flux" + ".png"
fig.savefig(save_path, dpi=500)

plt.close()


