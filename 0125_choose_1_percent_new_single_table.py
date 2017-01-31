
import astropy.io.fits as ts
from TheCannon_2 import dataset,apogee
from astropy.table import Table
import numpy as np
import os
from astropy.io import fits
import pickle
import AnniesLasso_2 as tc
import time

## training the cannon

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()



training_set_path = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/reference_labels.csv"
training_set_spectrum_dir = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/Data/"


training_set = Table.read("reference_labels.csv")



def get_pixmask(flux, err):
    bad_flux = ~np.isfinite(flux)
    bad_err = (~np.isfinite(err)) | (err <= 0)
    bad_pixels = bad_err | bad_flux
    return bad_pixels


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
    ivar = 1.0/flux_err**2
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
    training_set, training_set_flux, training_set_ivar,threads=8)
model.s2 = 0
model.regularization = 0
model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(training_set,
    tc.vectorizer.polynomial.terminator(("Teff_{corr}", "logg_{corr}", "[M/H]_{corr}"), 2))

model.train()



# function
def get_pixmask(flux, err):
    bad_flux = ~np.isfinite(flux)
    bad_err = (~np.isfinite(err)) | (err <= 0)
    bad_pixels = bad_err | bad_flux
    return bad_pixels

def get_error_mask(err):
    bad_err = (~np.isfinite(err))
    bad_pixels = bad_err
    return bad_pixels

# read data
all_set = Table.read("/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/allStar-v304.fits")

# data


## 1. read the dat
## 2. add a mask (bitmask)
## 3. filter-mean inverse variance







# choose some of them
choose = []
a=0
for i in range(0,4000):
    a += np.random.randint(1,14)
    choose.append(a)

choose = np.array(choose)
all_set = all_set[choose]


#Let's rock!

"""

N = len(all_set)

choose = []
a=0
for i in range(a,N):
    a += 1
    choose.append(a)

choose = np.array(choose)
all_set = all_set[choose]
"""




N = len(all_set)
fail = 0
star=0
success = 0
# choose some flux and normalize them first, then filter the mean inverse variance.


# get the time of the program
start_time = time.time()

# save the path of the flux fits and the path of the opt parameters fits

path_flux = np.array(["start"])
path_fits = np.array(["start"])


# mean ivar
mi = []


for i, row in enumerate(all_set):
    try:


        image_path = "/Volumes/Data_2TB/Data/APOGEE_DR10_Apstar/apStar-s3-" + row["APOGEE_ID"] + ".fits"


        RA = row["RA"]
        DEC = row["DEC"]
        SNR = row["SNR"]

        if not os.path.exists(image_path):
            print("{}/{} could not be found: {}".format(i + 1, N, image_path))
            fail += 1
            continue

        print("{}/{}: {}".format(i + 1, N, image_path))
        # Let's only use two extension of the data

        image = fits.open(image_path, ignore_missing_end=True)
        dat = Table.read(image_path)

        flux = image[1].data
        flux_err = image[2].data


    except IOError:
        print("opts. This one fail")
        fail+=1

    else:
        badpix = get_pixmask(flux, flux_err)
        ivar = 1.0 / flux_err ** 2
        error = flux_err
        # badpix is a array and the length is 8575
        flux = np.array(flux, dtype=np.float64)
        ivar = np.array(ivar, dtype=np.float64)

        flux[badpix] = 1.0
        ivar[badpix] = 0.0

        flux = np.array(flux)
        ivar = np.array(ivar)

        try:

            test_labels_all_i = [row["TEFF"], row["LOGG"], row["METALS"]]
            # print(test_labels_all_i)


        except ValueError:
            print("opts fail")
            fail+=1

        else:
            star += 1
            # Fiber number
            m = dat[0]["FIBER"]
            FiberID = m

            tr_ID = image_path
            # data
            MJD = dat[0]["MJD"]

            # APOGEE velocity

            VHELIO = dat[0]["VHELIO"]

            ## This is the real successful ones.

            ##let's optimize them one by one

            # normalization

            # trick to deal with 1D array:
            # by jason

            try:
                ni = len(flux[:, 0])
                one = np.ones(len(flux[0, :]))

            except IndexError:
                ni = 1
                one = np.ones(len(flux))

            else:
                flux = np.vstack((one, flux))
                ivar = np.vstack((one, ivar))

                ####
                # value

                ds = dataset.Dataset(wl, tr_ID, flux, ivar,
                                     test_labels_all_i, tr_ID, flux, ivar)

                ds.ranges = [[371, 3192], [3697, 5997], [6461, 8255]]

                # set sudo-continuous spectrum
                pseudo_tr_flux, pseudo_tr_ivar = ds.continuum_normalize_training_q \
                    (q=0.90, delta_lambda=50)

                # set mask
                contmask = ds.make_contmask(pseudo_tr_flux, pseudo_tr_ivar, frac=0.07)

                # get continuous mask

                ds.set_continuum(contmask)

                # fit the normalized-spectrum in the continuous region

                cont = ds.fit_continuum(3, "sinusoid")

                # Obtain the normalized flux
                norm_tr_flux, norm_tr_ivar, norm_test_flux, norm_test_ivar = \
                    ds.continuum_normalize(cont)

                # get inf_labels
                inf_labels = model.fit(norm_tr_flux, norm_tr_ivar)
                v = model.vectorizer.get_label_vector(inf_labels)
                inf_flux = np.dot(v, model.theta.T)
                opt_flux, parameters = model.fitting_spectrum_parameters_single \
                    (norm_tr_flux, norm_tr_ivar, inf_flux)

                opt_flux_5 , parameters_5 = model.fitting_spectrum_parameters_single_5(norm_tr_flux, norm_tr_ivar, inf_flux)

                parameters = np.c_[parameters,parameters_5]


                delta_chi_squared = model.delta_chi_squared(norm_tr_flux,norm_tr_ivar,inf_flux)

                chi_squared = model.chi_squared


                # by jason

                norm_tr_flux = norm_tr_flux[1:ni+1, :]
                norm_tr_ivar = norm_tr_ivar[1:ni+1, :]
                inf_labels = inf_labels [1:ni+1, :]
                delta_chi_squared = delta_chi_squared[1:ni+1]
                chi_squared = chi_squared[1:ni+1]

                inf_flux = inf_flux[1:ni+1, :]
                opt_flux = opt_flux[1:ni+1, :]
                parameters = parameters[1:ni+1, :]

                # others

                FiberID = np.array(FiberID)
                labels = np.array(test_labels_all_i)

                # calculate velocity


                velocity_correction = []
                n_star = len(parameters[:, 0])

                #for i in range(0, n_star):
                #    velocity_correction.append((parameters[i, 2] - parameters[i, 0]) * 4144.68)
                #velocity_correction = np.array(velocity_correction)
                for i in range(0, n_star):

                    # jason
                    v1 = 4144.68*(parameters[i,2]-parameters[i,0])/(parameters[i,0]+parameters[i,1]+parameters[i,2])
                    # david
                    v2 = 4144.68*(parameters[i,2]-parameters[i,0])/(-parameters[i,0]+2*parameters[i,1]-parameters[i,2])/2
                    # Jason_5_point
                    a = parameters[i,0]
                    b = parameters[i,1]
                    c = parameters[i,2]
                    a1 = parameters[i,3]
                    a2 = parameters[i, 4]
                    a3 = parameters[i, 5]
                    a4 = parameters[i, 6]
                    a5 = parameters[i, 7]
                    v3 = 4144.68*(2*a5+1*a4-1*a2-2*a1)/(a1+a2+a3+a4+a5)

                    if 2*b<a+c:
                        ve = v3
                    else:
                        ve = v1


                    #velocity_i = np.array([v1, v2, v3])
                    velocity_correction.append([v1,v3,ve])


                    #velocity_correction.append((parameters[i, 2] - parameters[i, 0]) * 4144.68/(4*(parameters[i,0]+parameters[i,2]-2*parameters[i,1])))
                velocity_correction = np.array(velocity_correction)

                # check shape
                print(inf_labels.shape,norm_tr_flux.shape,norm_tr_ivar.shape,
                      inf_flux.shape, opt_flux.shape, parameters.shape,
                      FiberID.shape,delta_chi_squared.shape,chi_squared.shape
                      ,velocity_correction.shape,labels.shape)

                ## save them

                ### remember to save them in single files.


                # nor flux 1
                a1 = norm_tr_flux

                # nor_ivar 2
                a2 = norm_tr_ivar

                #inf_flux 3

                a3 = inf_flux

                #opt_flux 4

                a4 = opt_flux

                #parameters 5

                a5 = parameters

                #chi-squared -6

                a6 = chi_squared

                # delta-chi-squared - 7

                a7 = delta_chi_squared

                # fiber_id - 8

                a8 = FiberID

                # labels -9

                a9 = labels

                #inf_labels -10

                a10 = inf_labels

                # velocity_shift -11

                a11 = velocity_correction



                # check mean ivar and labels

                ivar = norm_tr_ivar[0,:]

                mean_ivar = np.mean(ivar)
                teff = labels[0]
                logg = labels[1]
                metal = labels[2]

                if mean_ivar > 10000 and mean_ivar < 40000 and teff>0 and logg> -15 and metal > -15:
                #if mean_ivar > 10000 and mean_ivar < 40000:
                    ## save the fits files and the path files.

                    success += 1
                    mi.append(mean_ivar)

                    path_fits_i = "/Users/caojunzhi/Desktop/Data/n_600/" + row["APOGEE_ID"] + ".fits"

                    path_fits = np.append(path_fits,path_fits_i)

                    path_flux = np.append(path_flux,image_path)


                    # save fits

                    print("saving files"+path_fits_i)

                    prihdr = fits.Header()
                    prihdr['COMMENT'] = "The first four HDU lists are nor_flux, nor_ivar, inf_flux and opt_flux" \
                                    "Other information is in header"

                    prihdu = fits.PrimaryHDU(data=a1, header=prihdr)

                    # add header info
                    prihdu.header['SNR'] = SNR
                    prihdu.header['RA'] = RA
                    prihdu.header['DEC'] = DEC



                    # Let's write some columns

                    #MJD

                    c1 = fits.Column(name='MJD', format='PJ()',
                                 array=np.array([MJD],
                                                dtype=np.object))

                    #VHELIO
                    c2 = fits.Column(name='VHELIO', format="PJ()", array=np.array([VHELIO], dtype=np.object))


                    # a5-parameters

                    c3 = fits.Column(name='parameters', format='PJ()',
                                 array=np.array(parameters,
                                                dtype=np.object))

                    # a6-chi-squared

                    c4 = fits.Column(name='chi_squared', format='PJ()',
                                 array=np.array([chi_squared],
                                                dtype=np.object))

                    # a7-delta-chi_squared

                    c5 = fits.Column(name='d_chi_squared', format='PJ()',
                                 array=np.array([delta_chi_squared],
                                                dtype=np.object))

                    # a8-fiberID

                    c6 = fits.Column(name='FiberID', format='PJ()',
                                 array=np.array([FiberID],
                                                dtype=np.object))

                    # a9-labels

                    c7 = fits.Column(name='labels', format='PJ()',
                                 array=np.array([labels],
                                                dtype=np.object))

                    # a10 inf_labels

                    c8 = fits.Column(name='inf_labels', format='PJ()',
                                 array=np.array(inf_labels,
                                                dtype=np.object))

                    # a11 velocity_correction

                    c9 = fits.Column(name='vc', format='PJ()',
                                 array=np.array(velocity_correction,
                                                dtype=np.object))

                    tbhdu = fits.BinTableHDU.from_columns([c1, c2,c3,c4,c5,c6,c7,c8,c9])

                    thdulist = fits.HDUList([prihdu, tbhdu])


                    thdulist.writeto(path_fits_i, clobber=True)

                    ts.append(path_fits_i, a2)
                    ts.append(path_fits_i, a3)
                    ts.append(path_fits_i, a4)






                else:
                    print("mean ivar too big or too small")













N = len(path_fits)
path_fits = path_fits[1:N]
path_flux = path_flux[1:N]

path_fits = np.array(path_fits)

path_flux = np.array(path_flux)

mi = np.array(mi)


print(success)


# save them

output = open('n_600_path_fits.pkl', 'wb')
pickle.dump(path_fits, output)
output.close()

output = open('n_600_path_flux.pkl', 'wb')
pickle.dump(path_flux, output)
output.close()

output = open('n_600_choose.pkl', 'wb')
pickle.dump(choose, output)
output.close()

output = open('n_600_mean_ivar.pkl', 'wb')
pickle.dump(mi, output)
output.close()







# count for time



print("--- %s seconds ---" % (time.time() - start_time))










