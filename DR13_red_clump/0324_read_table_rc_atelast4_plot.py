import numpy as np
from astropy.table import Table
from astropy.io import fits

import matplotlib.pyplot as plt
import matplotlib
import pickle
from matplotlib import cm
from numpy.random import randn

# table path

path = "/Users/caojunzhi/Downloads/upload_20170322/red_clump_dr13.fits"


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

14. FIBER -- Fiber ID
15. SNR  -- SNR of the visit
####
The covariance matrix of the abc fit is in HDU0 data, which is
a 3*3*N 3-d matrix. N is the number of visits.
###
"""

# read covariance matrix from the abc fit:

un_cov = star[0].data[:,:,0]

#print(un_cov)


# read the velocity shift from the abc fit
v_shift = table["VSHIFT"]
#print(v_shift.shape)

########################
#Read table and plot to check.

class plot():

    def read_table(self):
        path = "/Users/caojunzhi/Downloads/upload_20170322/red_clump_atleast_4.fits"

        star = fits.open(path)
        table = Table.read(path)

        # read it:

        name = table["APOGEEID"]
        self.name = name

        un_cov = star[0].data
        self.un_cov = un_cov

        a = table["A"]
        b = table["B"]
        c = table["C"]
        self.a = a
        self.b = b
        self.c = c

        mask = 2*b>a+c
        self.mask = mask

        VBARY = table["VBARY"]
        self.VBARY =VBARY

        SHIFT = table["VSHIFT"]
        self.shift = SHIFT

        teff = table["TEFF"]
        self.teff = teff

        logg = table["LOGG"]
        self.logg = logg


        feh = table["FEH"]
        self.feh = feh

        self.chi_inf = table["CHIINF"]
        self.chi_mix = table["CHIMIX"]

        self.BJD = table["BJD"]

        self.fiber = table["FIBER"]

        self.SNR =table["SNR"]


    def plot_std_before_after(self):
        # From the average (c+a)/(a+b+c)
        # Do put a mask here
        mask = self.mask

        # add points with the same fiberid together
        name = self.name[mask]
        target = list(set(name))

        VBARY = self.VBARY[mask]
        shift =self.shift[mask]

        SNR = self.SNR[mask]

        fusion_new = []

        # name+std_old and std_new + SNR

        for i in range(0,len(target)):

            print("Doing %.2f %%"%(i/len(target)*100))

            index = np.where(name == target[i])
            index = np.array(index)
            index = index.ravel()

            std_old_i = np.std(VBARY[index])

            std_new_i = np.std(VBARY[index]+shift[index])

            SNR_i = np.nanmedian(SNR[index])

            fusion_new.append([target[i],std_old_i,std_new_i,SNR_i])

        fusion_new = np.array(fusion_new)



        # portion+fiber+rv

        name = fusion_new[:, 0]
        std_old = fusion_new[:,1]
        std_new = fusion_new[:,2]

        self.std_old = std_old
        self.std_new = std_new

        SNR = fusion_new[:,3]

        print("check shape")
        print(name.shape,std_old.shape,std_new.shape,SNR.shape)


        # let's plot

        font = {'family': 'normal',
                'weight': 'bold',
                'size': 14}

        matplotlib.rc('font', **font)

        plt.subplot(1,1,1)

        plt.plot(std_old,"ko",label="Before the correction",markersize=3,alpha=0.3)

        plt.xlabel("Visit")
        plt.ylabel("Standard deviation $km/s$")

        plt.plot(std_new, "rx", label="After the correction", markersize=3,alpha=0.3)

        plt.legend()


        # ax1

        """

        ax1.scatter(std_old, marker='x', c=SNR,
                    vmin=np.min(SNR), vmax=np.max(SNR), alpha=alpha, cmap=cm.coolwarm)



        ax1.scatter(std_new, marker='o', c=SNR,
                    vmin=np.min(SNR), vmax=np.max(SNR), alpha=alpha, cmap=cm.coolwarm)


        """

        plt.suptitle("Standard deviations of RV before and after the correction", fontsize=30)

        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(14.5, 8.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170322/" + "std_before_after_rc" + ".png"
        fig.savefig(save_path, dpi=500)

        plt.close()



    def diagnostic_plot(self):


        VBARY = self.std_old
        V_c = self.std_new

        plt.plot(VBARY,"ro")
        plt.plot(V_c,"ro")
        plt.show()

    def hist_rv_std(self):

        VBARY = self.std_old
        V_c = self.std_new

        VBARY = np.array(VBARY,dtype=float).ravel()

        V_c = np.array(V_c, dtype=float).ravel()


        font = {'weight': 'bold', 'size': 15}
        matplotlib.rc('font', **font)


        f, (ax1,ax2) = \
            plt.subplots(1,2)

        colors = ["cyan","r"]
        name = ["RV std before correction","RV std after correction"]

        rms_old = (np.sum(VBARY*VBARY)/len(VBARY))**0.5
        rms_new = (np.sum(V_c*V_c)/len(V_c)) ** 0.5

        ax1.hist(VBARY, bins=60,range=[-0.2,3], color=colors[0], label="%s RMS = %.2f $km/s$"%(name[0],rms_old))

        ax2.hist(V_c,bins=60,range=[-0.2,3], color=colors[1], label="%s RMS = %.2f $km/s$" % (name[1], rms_new))


        #ax1.set_title('Histogram of Radial velocity shifts', fontsize=30)
        ax1.set_xlabel('values of radial velocity std $ Km/s$', fontsize=15)
        ax1.set_ylabel('Number', fontsize=15)

        ax2.set_xlabel('values of radial velocity std $ Km/s$', fontsize=15)
        ax2.set_ylabel('Number', fontsize=15)

        ax1.legend(prop={'size': 15})
        ax2.legend(prop={'size': 15})

        # add vertical grey line
        # ax1.plot((wl[index], wl[index]), (0.5, 1 + 0.5 * N), 'k-', linewidth=1.5)

        # add vertical grey line
        # ax1.plot((wl[index], wl[index]), (0.5, 1 + 0.5 * N), 'k-', linewidth=1.5)
        f.suptitle("Histogram of RV std before and after the correction",fontsize=25)
        f.legends
        #f.suptitle("Histogram of RV shifts, a, b and c by using the absorption lines")


        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(18.5, 8.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170322/" + "histogram_rv_std_before_after" + ".png"
        fig.savefig(save_path, dpi=500)

        plt.close()












model = plot()
model.read_table()
model.plot_std_before_after()
model.hist_rv_std()
#model.diagnostic_plot()

