from astropy.io import fits
from astropy.table import Table
import numpy as np
import pickle
import os


# read data
red_clump = Table.read("/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/apogee-rc-DR13.fits")


N = len(red_clump)
print(N)

red_clump = red_clump

rc_path = []

for i,row in enumerate(red_clump):
    b = str(row["APOGEE_ID"])

    b = b+".fits"
    b = "/Volumes/Data_2TB/Data/DR13_rc/apStar-r6-"+b
    print(b)

    path = b

    """

    if not os.path.exists(path):
        print("{}/{} could not be found: {}".format(i + 1, N, path))

    """


    rc_path.append(path)

rc_path = np.array(rc_path)

print(len(rc_path))


output = open('Red_clump_dr13_ori.pkl', 'wb')
pickle.dump(rc_path, output)
output.close()




