import os
import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
import pickle

import astropy.io.fits as ts
import AnniesLasso_2 as tc


import math

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from math import sin, asin, cos, acos, pi

# calculate altaz from radectime


def radec_to_azalt(time, ra,dec):

    lon = -105.820417
    lat = 32.780361

    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')

    # in degree
    ra = c.ra.hour/24*360

    ha = time.sidereal_time("mean").hour/24*360 - ra;
    print(ra,time.sidereal_time("mean").hour/24*360,ha)

    if (ha < 0):
        ha = ha + 360


        #    print "ha:", ha
        # convert degrees to radians
    ha = ha * pi / 180
    dec = dec * pi / 180
    lat = lat * pi / 180

    # compute altitude in radians
    sin_alt = sin(dec) * sin(lat) + cos(dec) * cos(lat) * cos(ha);
    alt = asin(sin_alt)

    # compute azimuth in radians
    # divide by zero error at poles or if alt = 90 deg
    cos_az = (sin(dec) - sin(alt) * sin(lat)) / (cos(alt) * cos(lat))
    az = acos(cos_az)

    # convert radians to degrees
    hrz_altitude = alt * 180 / pi;
    hrz_azimuth = az * 180 / pi;

    # choose hemisphere
    if (sin(ha) > 0):
        hrz_azimuth = 360 - hrz_azimuth

    hrz_azimuth = hrz_azimuth * 24.0 / 360.0

    # return alt and az
    return hrz_altitude ,hrz_azimuth


def alt2airmass(alt):

    zenithAngle = 90.0 - alt
    za_radians = zenithAngle/180.0*math.pi
    airmass = 1.0/math.cos(za_radians)

    return airmass



class airmass():

    def cal_airmass(self,ra,dec,time):

        alt,az = radec_to_azalt(time,ra,dec)

        airmass = alt2airmass(alt)
        self.airmass = airmass
        print("alt = %.2f;az=%.2f;airmass=%.2f"%(alt,az,airmass))
        return airmass


pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()

"""

# import data
pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()



# load path

pkl_file = open('n_600_path_fits.pkl', 'rb')
path_fits = pickle.load(pkl_file)
pkl_file.close()


pkl_file = open('n_600_path_flux.pkl', 'rb')
path_flux = pickle.load(pkl_file)
pkl_file.close()


# mean_ivar

pkl_file = open('n_600_mean_ivar.pkl', 'rb')
mean_ivar = pickle.load(pkl_file)
pkl_file.close()




print(len(path_fits))


## read fluxes

nor_flux = np.ones((1,8575))
nor_ivar = np.ones((1,8575))
inf_flux = np.ones((1,8575))
parameters = np.array([0,1,0])


N = len(path_fits)
for i in range(0,N):
    print("loading star %d" % (i + 1))

    star_i = fits.open(path_fits[i])
    table_i = Table.read(path_fits[i])

    print(star_i[0].data.shape,star_i[1].data.shape,star_i[2].data.shape,star_i[3].data.shape)


    nor_flux = np.vstack((nor_flux,star_i[0].data))
    nor_ivar = np.vstack((nor_ivar, star_i[1].data))
    inf_flux = np.vstack((inf_flux, star_i[2].data))



    # load table part:
    parameters_i = np.array(list(table_i["parameters"]))

    parameters = np.vstack((parameters,parameters_i))
    print(parameters_i.shape)




print(nor_flux.shape,nor_ivar.shape,inf_flux.shape,parameters.shape)




"""

#########################



# load path

pkl_file = open('n_900_r_path_fits.pkl', 'rb')
path_fits = pickle.load(pkl_file)
pkl_file.close()


pkl_file = open('n_900_r_path_flux.pkl', 'rb')
path_flux = pickle.load(pkl_file)
pkl_file.close()


# mean_ivar

pkl_file = open('n_900_r_mean_ivar.pkl', 'rb')
mi = pickle.load(pkl_file)
pkl_file.close()

# choose some of them
path_fits = path_fits

N = len(path_fits)

# define a function which can derive the radial
# Generate a mask matrix, which has the very same shape with flux and ivar.

# set flux = 1 and ivar = 0 for masked pixels



def get_pixmask(flux, err):
    bad_flux = ~np.isfinite(flux)
    bad_err = (~np.isfinite(err)) | (err <= 0)
    bad_pixels = bad_err | bad_flux
    return bad_pixels



class fit_absorption_line():

    # Oh, mask is important
    # set pixel to be 10 and flux<0.95

    def find_mask(self,flux):

        # use 10 pixel:
        N_pixel = len(flux[0,:])
        N_star = len(flux[:,0])

        mask = []


        for i in range(0,N_star):

            print("Dealing with star %d"%(i+1))

            # one pixel is approximately 0.216\AA

            mask_i = np.zeros(N_pixel)

            width = 5
            limit = 0.95

            for j in range(15,N_pixel-15):
                flux_j = flux[i,j]
                min_j = np.min(flux[i,j-width:j+width])

                if flux_j <= min_j and flux_j < 1:
                    mask_i[j-width:j+width] = 1
                else:
                    nothing=1

            for j in range(0,N_pixel):
                flux_j = flux[i, j]
                if flux_j > limit:
                    mask_i[j] = 0
                else:
                    nothing=1

            mask.append(mask_i)
        mask = np.array(mask)
        print("mask shape")
        self.mask = mask

        print(mask.shape)
        return mask



    def train_cannon(self):

        ## training the cannon

        pkl_file = open('wl.pkl', 'rb')
        wl = pickle.load(pkl_file)
        pkl_file.close()

        self.wl = wl

        training_set_path = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/reference_labels.csv"
        training_set_spectrum_dir = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/Data/"

        training_set = Table.read("reference_labels.csv")

        N = len(training_set)
        keep = np.ones(N, dtype=bool)

        training_set_flux = []
        training_set_ivar = []
        training_set_error = []

        for i, row in enumerate(training_set):
            image_path = os.path.join(training_set_spectrum_dir, row["ID"])
            if not os.path.exists(image_path):
                print("{}/{} could not be found: {}".format(i + 1, N, image_path))
                keep[i] = False
                continue
            print("{}/{}: {}".format(i + 1, N, image_path))
            image = fits.open(image_path)
            flux = image[1].data
            flux_err = image[2].data
            badpix = get_pixmask(flux, flux_err)
            ivar = 1.0 / flux_err ** 2
            error = flux_err
            # badpix is a array and the length is 8575
            flux[badpix] = 1.0
            ivar[badpix] = 0.0
            training_set_flux.append(flux)
            training_set_ivar.append(ivar)
            training_set_error.append(error)

        training_set_flux = np.array(training_set_flux)
        training_set_ivar = np.array(training_set_ivar)
        training_set_error = np.array(training_set_error)

        assert all(keep)

        # Construct model.
        model = tc.L1RegularizedCannonModel(
            training_set, training_set_flux, training_set_ivar, threads=8)
        model.s2 = 0
        model.regularization = 0
        model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(training_set,
                                                                        tc.vectorizer.polynomial.terminator((
                                                                                                            "Teff_{corr}",
                                                                                                            "logg_{corr}",
                                                                                                            "[M/H]_{corr}"),
                                                                                                            2))

        model.train()

        self.model = model

    def model_cannon(self):


        ## training the cannon

        pkl_file = open('wl.pkl', 'rb')
        wl = pickle.load(pkl_file)
        pkl_file.close()

        self.wl = wl

        training_set_path = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/reference_labels.csv"
        training_set_spectrum_dir = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/Data/"

        training_set = Table.read("reference_labels.csv")

        N = len(training_set)
        keep = np.ones(N, dtype=bool)

        training_set_flux = []
        training_set_ivar = []
        training_set_error = []

        for i, row in enumerate(training_set):
            image_path = os.path.join(training_set_spectrum_dir, row["ID"])
            if not os.path.exists(image_path):
                print("{}/{} could not be found: {}".format(i + 1, N, image_path))
                keep[i] = False
                continue
            print("{}/{}: {}".format(i + 1, N, image_path))
            image = fits.open(image_path)
            flux = image[1].data
            flux_err = image[2].data
            badpix = get_pixmask(flux, flux_err)
            ivar = 1.0 / flux_err ** 2
            error = flux_err
            # badpix is a array and the length is 8575
            flux[badpix] = 1.0
            ivar[badpix] = 0.0
            training_set_flux.append(flux)
            training_set_ivar.append(ivar)
            training_set_error.append(error)

        training_set_flux = np.array(training_set_flux)
        training_set_ivar = np.array(training_set_ivar)
        training_set_error = np.array(training_set_error)

        assert all(keep)

        # Construct model.
        model = tc.L1RegularizedCannonModel(
            training_set, training_set_flux, training_set_ivar, threads=8)
        model.s2 = 0
        model.regularization = 0
        model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(training_set,
                                                                        tc.vectorizer.polynomial.terminator((
                                                                                                            "Teff_{corr}",
                                                                                                            "logg_{corr}",
                                                                                                            "[M/H]_{corr}"),
                                                                                                            2))

        # Do not need to train the model
        #model.train()



        self.model = model


    # opt is a all-in-one program
    # return flux_masked(call it flux_m),ivar_m, inf_m, parameters_m

    def opt(self,flux,ivar,inf_flux):

        self.flux = flux
        self.ivar = ivar
        self.inf_flux = inf_flux

        mask = self.find_mask(flux=flux)

        self.mask = mask

        N_star = len(mask[:,0])

        # simple code!


        one = flux/flux


        # flux m maintain the shape of flux. However, some elements are 0.

        flux_m = flux[mask]

        ivar_m = ivar[mask]

        inf_flux_m = inf_flux[mask]

        self.flux_m = flux_m
        self.ivar_m = ivar_m
        self.inf_flux_m = inf_flux_m

        ### Let's optimize:


        # get inf_labels
        # Shall we use inf_flux from flux_m or just old inf_flux_m??

        # When you optimize, only use part of the spectrum!!


        #inf_labels = self.model.fit(flux_m, ivar_m)
        #v = self.model.vectorizer.get_label_vector(inf_labels)
        #inf_flux = np.dot(v, self.model.theta.T)
        # attention! only change ivar. Do no touch inf_flux/flux!!
        opt_flux_m, parameters_m = self.model.fitting_spectrum_parameters_single \
            (flux, ivar_m, inf_flux)

        self.parameters_m = parameters_m
        self.opt_flux_m = opt_flux_m

        return parameters_m



    def test(self):

        flux = self.flux
        flux_m = self.flux_m

        N = len(flux[:,0])

        wl = self.wl

        for i in range(0,N):

            plt.plot(wl,flux[i,:],"k")
            plt.plot(wl,flux_m[i,:],"r",alpha=0.5)
            plt.xlim([15660,15780])
            plt.show()

    def test_para(self,parameters):
        parameters_m =self.parameters_m
        plt.plot(np.ravel(parameters_m),np.ravel(parameters),"ro")
        plt.show()

    def save_single_file(self,path,path_origin):
        N = len(path)
        for i in range(0,N):

            """

            star_i = fits.open(path[i])
            data_i = fits.open(path_origin[i])

            flux = star_i[0].data
            ivar = star_i[1].data
            inf_flux = star_i[2].data

            mask = self.find_mask(flux=flux)

            # simple code!

            one = flux / flux

            flux_m = (flux - one) * mask + one

            ivar_m = ivar * mask

            inf_flux_m = (inf_flux - one) * mask + one

            # Don't use the inf_flux_m...

            # optimize

            opt_flux_m, parameters_m = self.model.fitting_spectrum_parameters_single \
                (flux_m, ivar_m, inf_flux_m)


            """

            star_i = fits.open(path[i])
            data_i = fits.open(path_origin[i])

            flux = star_i[0].data
            ivar = star_i[1].data
            inf_flux = star_i[2].data
            ra = star_i[0].header["RA"]
            dec = star_i[0].header["DEC"]

            parameters_m = self.opt(flux=flux,ivar=ivar,inf_flux=inf_flux)




            # calculate velocity_m


            velocity_m = []
            velocity_uncertainty = self.model.uncertainty
            n_star = len(parameters_m[:, 0])

            # for i in range(0, n_star):
            #    velocity_correction.append((parameters[i, 2] - parameters[i, 0]) * 4144.68)
            # velocity_correction = np.array(velocity_correction)
            for j in range(0, n_star):
                a = parameters_m[j, 0]
                b = parameters_m[j, 1]
                c = parameters_m[j, 2]

                # jason
                v1 = 4144.68 * (c - a) / (a + b + c)
                # david
                v2 = velocity_uncertainty[j]
                # uncertainty is v2

                velocity_m.append([v1, v2])


                # velocity_correction.append((parameters[i, 2] - parameters[i, 0]) * 4144.68/(4*(parameters[i,0]+parameters[i,2]-2*parameters[i,1])))
            velocity_m = np.array(velocity_m)

            # Read BJD for individual visits
            BJD = []
            AirMass = []
            ai = airmass()
            for k in range(0, len(data_i[1].data[:, 0]) - 2):
                name = "HJD" + str(k + 1)

                BJD_i = data_i[0].header[name]
                BJD.append(BJD_i)

                # meanwhile calculate the airmass
                t_jd = Time(BJD_i, format='jd', scale='utc', location=(-105.820417, 32.780361))
                airmass_i =ai.cal_airmass(ra,dec,t_jd)
                AirMass.append(airmass_i)


            BJD = np.array(BJD)
            AirMass = np.array(AirMass)

            # calculate airmass:
            # waiting

            # check shape

            mask = self.mask

            print(flux.shape,mask.shape,parameters_m.shape,velocity_m.shape,BJD.shape)


            # save
            print("saving path %d"%(i+1))
            print(path[i])


            ts.append(path[i], mask)
            ts.append(path[i], parameters_m)
            ts.append(path[i], velocity_m)
            ts.append(path[i], BJD)
            ts.append(path[i], AirMass)




my_model = fit_absorption_line()


# opt
my_model.train_cannon()

my_model.save_single_file(path=path_fits,path_origin=path_flux)

# my_model.test_para(parameters=parameters)



#my_model.train_cannon()
#my_model.save_single_file(path=path_fits,path_origin=path_flux)



"""

# train
my_model.train_cannon()

# opt
my_model.opt(flux=nor_flux,ivar=nor_ivar,inf_flux=inf_flux)

# test
# my_model.test()

# my_model.test_para(parameters=parameters)

# save parameters and compare


parameters_m = my_model.parameters_m

mask = my_model.mask

output = open('parameters_m_10_pixels.pkl', 'wb')
pickle.dump(parameters_m, output)
output.close()


output = open('mask_10_pixels.pkl', 'wb')
pickle.dump(mask, output)
output.close()


"""