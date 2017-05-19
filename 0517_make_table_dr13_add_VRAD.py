import numpy as np
import pickle
import time

from astropy.io import fits
from astropy.table import Table

import math

from astropy import units as u
from astropy.coordinates import SkyCoord
from math import sin, asin, cos, acos, pi
from astropy.time import Time


def radec_to_azalt(time, ra, dec):
    lon = -105.820417
    lat = 32.780361

    c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')

    # in degree
    ra = c.ra.hour / 24 * 360

    ha = time.sidereal_time("mean").hour / 24 * 360 - ra;
    print(ra, time.sidereal_time("mean").hour / 24 * 360, ha)

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
    return hrz_altitude, hrz_azimuth


def alt2airmass(alt):
    zenithAngle = 90.0 - alt
    za_radians = zenithAngle / 180.0 * math.pi
    airmass = 1.0 / math.cos(za_radians)

    return airmass


# t_jd = Time(2455197.5, format='jd', scale='utc', location=(-105.820417,32.780361))


class airmass():
    def cal_airmass(self, ra, dec, time):
        alt, az = radec_to_azalt(time, ra, dec)

        airmass = alt2airmass(alt)
        self.airmass = airmass
        print("alt = %.2f;az=%.2f;airmass=%.2f" % (alt, az, airmass))
        return airmass


fail = 0


class make_table():
    def read_data(self):

        pkl_file = open(
            '/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/dr13_fits.pkl', 'rb')
        path_fits = pickle.load(pkl_file)
        pkl_file.close()

        pkl_file = open(
            '/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/dr13_ori.pkl', 'rb')
        path_ori = pickle.load(pkl_file)
        pkl_file.close()

        # Diagnostic


        
        """
        
        NN = 1000
        path_fits = path_fits[1:NN]
        path_ori = path_ori[1:NN]

        
        
        """



        







        path_all = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/allStar-l30e.2.fits"

        table_all = Table.read(path_all, hdu=1)

        N = len(path_fits)
        print(N)

        star_name = []
        star_visit = []

        vsini = []
        vmicro = []
        vmacro = []

        BJD = []
        inf_labels = np.array([0, 0, 0])
        chi_inf = []
        chi_mix = []

        parameters = np.array([0, 1, 0])

        un_cov = np.zeros((3, 3))

        fiber = []
        VBARY = []
        SHIFT = []

        SNR = []
        RA = []
        DEC = []
        jd = []

        ## I don't think you have every visit, thus you can set VRAD = -9999 if the file doesn't exist.
        vrad = []
        xshift = []

        for i in range(0, N):

            try:

                print("reading star %d  %.2f %%" % (i + 1, (i + 1) / N * 100))

                star = fits.open(path_fits[i])

                dat = Table.read(path_ori[i])

                ni = len(star[1].data)

                star_name_i = path_fits[i].replace(".fits ", "")
                star_name_i = star_name_i.replace("/Volumes/Data_2TB/Data/dr13_all/", "")

                index_name = np.where(table_all["APOGEE_ID"] == star_name_i)

                apvisit = str(table_all[index_name]["ALL_VISITS"]).split(",")
                #print(apvisit)

                vsini_i = table_all[index_name]["VSINI"]
                vmicro_i = table_all[index_name]["VMICRO"]
                vmacro_i = table_all[index_name]["VMACRO"]

                ap_path = "/Volumes/Data_2TB/Data/APOGEE_DR13_Apvisit/apVisit-"



                SNR_i = []
                RA_i = []
                DEC_i = []
                star_name_ii = []
                vsini_ii = []
                vmicro_ii = []
                vmacro_ii = []
                vrad_ii = []
                xshift_ii = []

                for j in range(0, ni - 2):
                    star_name_ii = np.append(star_name_ii, star_name_i)
                    vsini_ii = np.append(vsini_ii, vsini_i)
                    vmicro_ii = np.append(vmicro_ii, vmicro_i)
                    vmacro_ii = np.append(vmacro_ii, vmacro_i)

                    SNR_i = np.append(SNR_i, star[0].header["SNR"])
                    # add RA DEC
                    RA_i = np.append(RA_i, star[0].header["RA"])
                    DEC_i = np.append(DEC_i, star[0].header["DEC"])

                    # Try apvisit:

                    #print("check path")
                    #print(ap_path + str(apvisit[j]) + ".fits")

                    try:
                        visit_i = fits.open(ap_path+str(apvisit[j])+".fits")

                        vrad_ii.append(visit_i[0].header["VRAD"])
                        xshift_ii.append(visit_i[0].header["XSHIFT"])
                    except:
                        vrad_ii.append(-9999)
                        xshift_ii.append(-9999)
                        #print("fail")




                star_visit_i = np.array(dat[0]["FILE"])
                jd_i = np.array(dat[0]["JD"])

                chi_inf_i = star[9].data[2:ni]

                chi_mix_i = star[10].data[2:ni]

                BJD_i = star[13].data

                fiber_i = star[12].data

                VBARY_i = star[11].data

                SHIFT_i = star[1].data[2:ni]

                inf_labels_i = star[8].data[2:ni, :]

                parameters_i = star[0].data[2:ni, :]

                un_cov_i = star[3].data[:, :, 2:ni]



            except IndexError:
                print("This one fails")



            except OSError:
                print("This one fails")




            else:
                SNR = np.append(SNR, SNR_i)
                RA = np.append(RA, RA_i)
                DEC = np.append(DEC, DEC_i)

                star_visit = np.append(star_visit, star_visit_i)
                jd = np.append(jd, jd_i)

                star_name = np.append(star_name, star_name_ii)
                vsini = np.append(vsini, vsini_ii)
                vmicro = np.append(vmicro, vmicro_ii)
                vmacro = np.append(vmacro, vmacro_ii)
                vrad = np.append(vrad,vrad_ii)
                xshift = np.append(xshift,xshift_ii)


                chi_inf = np.append(chi_inf, chi_inf_i)

                chi_mix = np.append(chi_mix, chi_mix_i)

                BJD = np.append(BJD, BJD_i)

                fiber = np.append(fiber, fiber_i)

                VBARY = np.append(VBARY, VBARY_i)

                SHIFT = np.append(SHIFT, SHIFT_i)

                inf_labels = np.vstack((inf_labels, inf_labels_i))

                parameters = np.vstack((parameters, parameters_i))

                un_cov = np.dstack((un_cov, un_cov_i))

        n_visit = len(un_cov[0, 0, :])

        star_name = np.array(star_name)
        self.star_name = star_name

        vsini = np.array(vsini)
        self.vsini = vsini

        vmicro = np.array(vmicro)
        self.vmicro = vmicro

        vmacro = np.array(vmacro)
        self.vmacro = vmacro

        vrad = np.array(vrad)
        self.vrad = vrad

        xshift = np.array(xshift)
        self.xshift = xshift

        SNR = np.array(SNR)
        self.SNR = SNR

        star_visit = np.array(star_visit)
        star_visit = star_visit.ravel()
        self.star_visit = star_visit

        chi_inf = np.array(chi_inf)
        chi_inf = chi_inf.ravel()
        self.chi_inf = chi_inf

        chi_mix = np.array(chi_mix)
        chi_mix = chi_mix.ravel()
        self.chi_mix = chi_mix

        BJD = np.array(BJD).ravel()
        self.BJD = BJD

        fiber = np.array(fiber).ravel()
        self.fiber = fiber

        VBARY = np.array(VBARY).ravel()
        self.VBARY = VBARY

        SHIFT = np.array(SHIFT).ravel()
        self.SHIFT = SHIFT

        inf_labels = inf_labels[1:n_visit, :]
        self.inf_labels = inf_labels

        parameters = parameters[1:n_visit, :]
        self.parameters = parameters

        un_cov = un_cov[:, :, 1:n_visit]
        self.un_cov = un_cov

        # Add RA DEC JD and Airmass

        self.RA = np.array(RA)
        self.DEC = np.array(DEC)
        self.jd = np.array(jd)

        md = airmass()

        am = []

        n_visit = n_visit - 1

        for vi in range(0, n_visit):
            t_jd = Time(jd[vi], format='jd', scale='utc', location=(-105.820417, 32.780361))
            am_i = md.cal_airmass(RA[vi], DEC[vi], t_jd)
            am.append(am_i)
            print("calculating airmass %.2f percent" % (vi / n_visit * 100))

        self.am = np.array(am)

        print("check shape")

        print(n_visit)

        print(star_name.shape, SNR.shape, star_visit.shape, chi_inf.shape, chi_mix.shape)
        print(BJD.shape, fiber.shape, inf_labels.shape, parameters.shape, un_cov.shape, VBARY.shape, SHIFT.shape)

        print(np.array(RA).shape, np.array(DEC).shape, np.array(jd).shape, np.array(am).shape)

        print(vsini.shape, vmicro.shape, vmacro.shape,vrad.shape,xshift.shape)

    # Input the path of the table:

    def write_table(self, path):

        # save them in the header




        self.table_path = path

        prihdr = fits.Header()
        prihdr['COMMENT'] = "id MJD HJD VELIO VELIOUN Shift ShiftUN Shiftabs ShiftabsUN"

        prihdu = fits.PrimaryHDU(data=self.un_cov, header=prihdr)

        # Table list


        col1 = fits.Column(name='APOGEEID', format="25A", array=self.star_name)
        col2 = fits.Column(name='VISIT', format="25A", array=self.star_visit)
        col3 = fits.Column(name='BJD', format="E", array=self.BJD)

        col4 = fits.Column(name='TEFF', format='E', array=self.inf_labels[:, 0])
        col5 = fits.Column(name='LOGG', format="E", array=self.inf_labels[:, 1])
        col6 = fits.Column(name='FEH', format='E', array=self.inf_labels[:, 2])

        col7 = fits.Column(name='A', format="E", array=self.parameters[:, 0])
        col8 = fits.Column(name='B', format="E", array=self.parameters[:, 1])
        col9 = fits.Column(name='C', format="E", array=self.parameters[:, 2])

        col10 = fits.Column(name='CHIINF', format="E", array=self.chi_inf)
        col11 = fits.Column(name='CHIMIX', format="E", array=self.chi_mix)

        col12 = fits.Column(name='VBARY', format="E", array=self.VBARY)
        col13 = fits.Column(name='VSHIFT', format="E", array=self.SHIFT / 1000)

        col14 = fits.Column(name='FIBER', format="E", array=self.fiber)
        col15 = fits.Column(name='SNR', format="E", array=self.SNR)

        # add RA, DEC and AirMass

        col16 = fits.Column(name='RA', format="E", array=self.RA)

        col17 = fits.Column(name='DEC', format="E", array=self.DEC)
        col18 = fits.Column(name='AirMass', format="E", array=self.am)

        # Add vsini, vmicro and vmacro


        col19 = fits.Column(name='VSINI', format="E", array=self.vsini)

        col20 = fits.Column(name='VMICRO', format="E", array=self.vmicro)
        col21 = fits.Column(name='VMACRO', format="E", array=self.vmacro)
        col22 = fits.Column(name='VRAD', format="E", array=self.vrad)
        col23 = fits.Column(name='XSHIFT', format="E", array=self.xshift)

        cols = fits.ColDefs(
            [col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16,
             col17, col18, col19, col20, col21,col22,col23])

        tbhdu = fits.BinTableHDU.from_columns(cols)

        thdulist = fits.HDUList([prihdu, tbhdu])

        print("saving table")

        thdulist.writeto(path, clobber=True)
        thdulist.close()

    def check(self, path):

        star = fits.open(path)
        table = Table.read(path)
        self.star = star
        self.table = table

        print(table[0])
        print(star[0].data.shape)


# construct table:

start_time = time.time()

model = make_table()
model.read_data()

path = "/Users/caojunzhi/Downloads/upload_0516_David/dr13.fits"
model.write_table(path=path)
# model.check(path)

print("The number of fails %d" % fail)

stop_time = time.time()

print("The time we use is %.2f" % (stop_time - start_time))

