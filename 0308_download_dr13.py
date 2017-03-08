
from astropy.table import Table

import urllib.request



# read data
all_set = Table.read("allStar-l30e.2.fits",hdu=1)
#image=fits.open("allStar-v603.fits")

print(all_set[2])



N = len(all_set)
all_set = all_set[0:100]

for i,row in enumerate(all_set):
    b = str(row["FILE"])
    #path = "curl -o " + b + " http://data.sdss3.org/sas/dr13/apogee/spectro/redux/r6/stars/apo25m/"+str(row["LOCATION_ID"])+"/" + b
    path = "http://data.sdss3.org/sas/dr13/apogee/spectro/redux/r6/stars/apo25m/" + str(row["LOCATION_ID"]) + "/" + b
    save = "/Volumes/Data_2TB/Data/APOGEE_DR13_Apstar/"
    loc = save+b


    try:
        urllib.request.urlretrieve(path, loc)
        print("doing star %i" % i)

    except urllib.error.HTTPError:
        print("%s not found"%b)








