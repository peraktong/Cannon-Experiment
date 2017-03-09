import numpy as np


import astropy.io.fits as ts
from TheCannon_2 import dataset,apogee
from astropy.table import Table
import numpy as np
import os
from astropy.io import fits
import pickle
import AnniesLasso_2 as tc
import time




# read data
all_set = Table.read("allStar-l30e.2.fits",hdu=1)
#image=fits.open("allStar-v603.fits")

print(all_set[2])


x= []
N = len(all_set)
for i,row in enumerate(all_set):
    b = str(row["FILE"])
    #path = "curl -o " + b + " http://data.sdss3.org/sas/dr13/apogee/spectro/redux/r6/stars/apo25m/"+str(row["LOCATION_ID"])+"/" + b
    path = "http://data.sdss3.org/sas/dr13/apogee/spectro/redux/r6/stars/apo25m/" + str(row["LOCATION_ID"]) + "/" + b
    x.append(path)
    print("doing star %i"%i)


print(type(x))


thefile = open("APOGEE_dr13_fits_web.txt","w")
for item in x:
  thefile.write("%s\n" % item)





