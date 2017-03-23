from astropy.io import fits
from astropy.table import Table
import numpy as np
import pickle
import os

import astropy.io.fits as ts
from TheCannon_2 import dataset,apogee
from astropy.table import Table
import numpy as np
import os
from astropy.io import fits
import pickle
import AnniesLasso_2 as tc
import time
from os import listdir
from os.path import isfile, join


mypath = "/Volumes/Data_2TB/Data/DR13_rc/"

big_path = "/Volumes/Data_2TB/Data/APOGEE_DR13_Apstar/"

# read data

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

N = len(onlyfiles)
onlyfiles = onlyfiles[1:N]
print(N)


pkl_file = open('Red_clump_dr13_name.pkl', 'rb')
rc_name = pickle.load(pkl_file)
pkl_file.close()

red_clump_ori = []

# What???

for i in range(0,50):
    print("doing star %d"%i)


    image_path = mypath+rc_name[i]
    print(image_path)

    if not os.path.exists(image_path):
        print("{}/{} could not be found: {}".format(i + 1, N, image_path))


    big_i = big_path+rc_name[i]
    red_clump_ori.append(big_i)

    print(big_i)

    if not os.path.exists(big_i):
        print("there is no original one %d"%i)

big_i = np.array(big_i)


output = open('Red_clump_dr13.pkl', 'wb')
pickle.dump(big_i, output)
output.close()




