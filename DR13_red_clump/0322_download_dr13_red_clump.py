
from astropy.table import Table

import urllib.request



# read data
red_clump = Table.read("/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/apogee-rc-DR13.fits")
#image=fits.open("allStar-v603.fits")





N = len(red_clump)
print(N)

red_clump = red_clump[0:100]


for i,row in enumerate(red_clump):
    b = str(row["APOGEE_ID"])
    b = b+".fits"
    b = "apStar-r6-" +b
    print(b)

    #path = "curl -o " + b + " http://data.sdss3.org/sas/dr13/apogee/spectro/redux/r6/stars/apo25m/"+str(row["LOCATION_ID"])+"/" + b
    path = "http://data.sdss3.org/sas/dr13/apogee/spectro/redux/r6/stars/apo25m/" + str(row["LOCATION_ID"]) + "/" + b
    save = "/Volumes/Data_2TB/Data/DR13_rc/"
    loc = save+b


    try:
        urllib.request.urlretrieve(path, loc)
        print("doing star %i" % i)

    except urllib.error.HTTPError:
        print("%s not found"%b)








