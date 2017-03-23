
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time




## attention!!
# The observer locations are different for different telescopes

"""
APOGEE-Mexicon

    lon = -105.820417
    lat = 32.780361

('-0.0077d', '51.4826')

"""


def MJD2BJD(mjd, target_coord, site_location=(-105.820417,32.780361)):
    """do the conversion
mjd -- input mjd, scale is utc
target_coord -- the coordinate of the target, in astropy.coord format,
to caculate the light travel time
site_location -- location of the telescope, to make accurate calcualtion of
light travel time and tdb conversion. Not very important here. The default value
is for Greenwich Observatory.
"""
    t = Time(mjd, format='mjd', scale='utc', location=site_location)
    # calculate light travel time
    ltt = t.light_travel_time(target_coord)
    # print(t, ltt)
    # convert t to tdb, and add light travel time
    t_out = (t.tdb + ltt).jd

    return t_out


star = fits.open("apStar-r5-2M00005143+5615568.fits")
dat = Table.read("apStar-r5-2M00005143+5615568.fits")

RA = star[0].header["RA"]
DEC = star[0].header["DEC"]
MJD = dat[0]["MJD"]
JD = dat[0]["JD"]

c = SkyCoord(RA,DEC, frame='icrs', unit='deg')

BJD = MJD2BJD(MJD,c)


print(BJD)
print(JD)

