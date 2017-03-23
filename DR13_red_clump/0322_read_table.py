import numpy as np
from astropy.table import Table
from astropy.io import fits

# table path

path = "/Users/caojunzhi/Downloads/upload_20170322/table.fits"


star = fits.open(path)
table = Table.read(path)

"""
There are 13 columns in the table:

1. 'APOGEEID' -- The name of the star
2. 'VISIT' -- The name of the visit file
3. BJD -- Barycentric JD

Inferred labels are from the Cannon. The spectra we use are from the first combined spectra
(There are two combined spectra for each star, which are obtained by two different methods)

: (1) global weighting, where each visit spectrum is weighted by its (S/N)2, and
(2) pixel-by-pixel weighting, where each pixel is weighted by its (S/N)2.

4. TEFF
5. LOGG
6. FEH

The abc parameters for each visit:

7. A -- parameter a
8. B -- parameter b
9. C -- parameter c

10. CHIINF -- chi-squared for the inferred flux from the cannon (a=0,b=1,c=0)
11. CHIMIX -- chi-squared for the mixed flux from the abc fit.

12. VBARY -- The barycentric Velocity(km/s) from the APOGEE team.
13. VSHIFT -- The velocity shift from the abc fit(km/s)

####
The covariance matrix of the abc fit is in HDU0 data, which is
a 3*3*N 3-d matrix. N is the number of visits.
###
"""

# read covariance matrix from the abc fit:

un_cov = star[0].data[:,:,0]

print(un_cov)


# read the velocity shift from the abc fit
v_shift = table["VSHIFT"]
print(v_shift.shape)

###