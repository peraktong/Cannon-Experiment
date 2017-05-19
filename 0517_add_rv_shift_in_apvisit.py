import numpy as np
from astropy.io import fits
from astropy.table import Table

alls = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/allStar-l30e.2.fits"

path = "/Users/caojunzhi/Downloads/apVisit-r6-5094-55874-088.fits"

ap_star = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/apStar-r5-2M00005143+5615568.fits"

star_visit = fits.open(path)

star_ap =  fits.open(ap_star)
table = Table.read(ap_star)

VHELIO = star_visit[0].header["VHELIO"]

print(star_visit[0].header["VRAD"])

alls = Table.read(alls)


print(alls[22]["ALL_VISITS"])

