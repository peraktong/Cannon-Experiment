from astropy.table import Table
from astropy.io.fits import getdata
from astropy.io import fits

image_path = "/Users/caojunzhi/Desktop/Data/APOGEE_DR10_Apstar/apStar-s3-2M00000068+5710233.fits"
image = fits.open(image_path,ignore_missing_end=True)

dat = Table.read(image_path)

#print(dat[0]["FIBER"],sum(dat[0]["FIBER"])/len(dat[0]["FIBER"]))
print(dat[0]["FIBER"])
