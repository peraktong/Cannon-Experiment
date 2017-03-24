import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
import pickle
from os.path import isfile, join
from os import listdir
import os


from astropy.time import Time


mypath = "/Users/caojunzhi/Desktop/Data/dr13_red_clump/"
mypath_ori = "/Volumes/Data_2TB/Data/DR13_rc/apStar-r6-"

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]


rc_path = []
rc_ori = []

for i in range(1,len(onlyfiles)):

    print("doing star %d"%(i+1))

    path_i = mypath+onlyfiles[i]
    ori_i = mypath_ori+onlyfiles[i]

    em=1

    if not os.path.exists(path_i):
        print("{}/{} could not be found: {}".format(i + 1, len(onlyfiles), path_i))
        em=0

    """

    if not os.path.exists(ori_i):
        print("{}/{} could not be found: {}".format(i + 1, len(onlyfiles), ori_i))
        em=0


    """

    try:
        star = fits.open(path_i)
        visit = len(star[1].data)
    except IndexError:
        print("This one failed")
        em=0

    except OSError:
        print("This one failed")
        em=0

    if visit>3:
        nothing=1
    else:
        em = 0


    if em==1:
        rc_path.append(mypath + onlyfiles[i])
        rc_ori.append(mypath_ori + onlyfiles[i])
    else:
        print("This one fail")

rc_path = np.array(rc_path)
rc_ori = np.array(rc_ori)


print(len(rc_path))
print(len(rc_ori))


output = open('dr13_rc_fits_atleast_4_epoch.pkl', 'wb')
pickle.dump(rc_path, output)
output.close()

output = open('dr13_rc_ori_atleast_4_epoch.pkl', 'wb')
pickle.dump(rc_ori, output)
output.close()





