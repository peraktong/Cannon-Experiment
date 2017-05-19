import numpy as np
from astropy.table import Table
from astropy.io import fits
import pickle


path = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/allStar-l30e.2.fits"

table = Table.read(path,hdu=1)


pkl_file = open('/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/dr13_fits.pkl',
                'rb')
path_fits = pickle.load(pkl_file)
pkl_file.close()

# print(path_fits)


# print(table[0])

print(table["VSINI"])

"""

print(table[100]["VSINI"])



for i in range(0,len(path_fits)):

    path_fits_i = path_fits[i].replace("/Volumes/Data_2TB/Data/dr13_all/","")
    path_fits_i = path_fits_i.replace(".fits ","")

    print(path_fits_i)

    index = np.where(table["APOGEE_ID"] == path_fits_i)

    print(table[index]["VSINI"])
    print(table[index]["VMICRO"])
    print(table[index]["VMACRO"])



"""




