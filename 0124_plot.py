import os
import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
import pickle

from TheCannon_2 import dataset,apogee
from TheCannon_2 import model


pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()


# load path

pkl_file = open('n_900_path_fits.pkl', 'rb')
path_fits = pickle.load(pkl_file)
pkl_file.close()


pkl_file = open('n_900_path_flux.pkl', 'rb')
path_flux = pickle.load(pkl_file)
pkl_file.close()


# mean_ivar

pkl_file = open('n_900_mean_ivar.pkl', 'rb')
mi = pickle.load(pkl_file)
pkl_file.close()



N = len(path_fits)
print(N)

velocity = []
fiber_id = []
mean_ivar = []
parameters = []
a = []
b = []
c= []

nor = []
ivar=[]
inf =[]
inf_label= []
dchi = []

MJD = []
meanivar = []

RA = []
DEC = []
SNR = []


for i in range(0,N):
    print("loading star %d"%(i+1))

    star_i = fits.open(path_fits[i])


    ni = len(star_i[4].data[:, 0])

    #mean ivar
    one = np.ones(ni-2)
    meanivar = np.append(meanivar,one*mi[i])

    dchi = np.append(dchi,star_i[6].data[2:ni])
    # SNR RA DEC
    SNR = np.append(SNR,(star_i[0].header["SNR"]*one))
    RA = np.append(RA, (star_i[0].header["RA"] * one))
    DEC = np.append(DEC, (star_i[0].header["DEC"] * one))



    velocity = np.append(velocity,star_i[10].data[2:ni,0])
    fiber_id = np.append(fiber_id,star_i[7].data)
    mean_ivar.append(np.mean(star_i[1].data[0]))
    parameters = np.append(parameters,star_i[4].data[2:ni,2])

    a = np.append(a, star_i[4].data[2:ni, 0])
    b = np.append(b, star_i[4].data[2:ni, 1])
    c = np.append(c, star_i[4].data[2:ni, 2])


    nor.append(star_i[0].data[0])
    ivar.append(star_i[1].data[0])
    inf.append(star_i[2].data[0])
    inf_label.append(star_i[9].data[0])

    MJD = np.append(MJD,star_i[11].data)

    print(star_i[4].data[:,0].shape,star_i[0].data.shape,star_i[11].data.shape,star_i[12].data.shape)
    print(star_i[11].data)

velocity =np.array(velocity)
fiber_id = np.array(fiber_id)
parameters = np.array(parameters)
nor_flux=np.array(nor)
ivar_flux = np.array(ivar)
inf_flux = np.array(inf)
inf_label = np.array(inf_label)
dchi = np.array(dchi)


MJD = np.array(MJD)
SNR = np.array(SNR)
RA = np.array(RA)
DEC = np.array(DEC)

meanivar = np.array(meanivar)

# check shape

print(MJD.shape,SNR.shape,RA.shape,DEC.shape)

# plot velocity vs MJD

class plot():


    # for ve


    def ve_subplot_4(self,velocity,MJD,Fiber,RA,DEC,mean_ivar):




        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)

        fig = plt.figure()


        f, ((ax1, ax2), (ax3, ax4)) = \
            plt.subplots(2, 2)

        #ax1

        ax1.scatter(MJD, velocity, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax1.set_title('RV shifts vs MJD', fontsize=24,y=0.90)

        ax1.set_xlabel('MJD', fontsize=30)
        ax1.set_ylabel('RV shifts $m/s$', fontsize=30)


        ax1.set_ylim([-6000,8000])
        ax1.set_yticks(np.arange(-6000,8001,3500))

        #ax1.set_position([0,0.6,0.4,0.4])


        #ax2

        ax2.scatter(RA, velocity, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax2.set_title('RV shifts vs RA', fontsize=24,y=0.90)

        ax2.set_xlabel('RA', fontsize=30)
        ax2.set_ylabel('RV shifts $m/s$', fontsize=30)

        #ax2.set_position([0.5, 0.6, 0.4, 0.4])

        ax2.set_ylim([-6000, 8000])
        ax2.set_yticks(np.arange(-6000, 8001, 3500))



        #ax3

        ax3.scatter(DEC, velocity, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax3.set_title('RV shifts vs DEC', fontsize=24,y=0.90)

        ax3.set_xlabel('DEC', fontsize=30)
        ax3.set_ylabel('RV shifts $m/s$', fontsize=30)

        ax3.set_ylim([-6000, 8000])
        ax3.set_yticks(np.arange(-6000, 8001, 3500))

        #ax3.set_position([0, 0, 0.4, 0.4])



        #ax4

        ax4.scatter(Fiber, velocity, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax4.set_title('RV shifts vs FiberID', fontsize=24,y=0.90)

        ax4.set_xlabel('FIberID', fontsize=30)
        ax4.set_ylabel('RV shifts $m/s$', fontsize=30)

        ax4.set_ylim([-6000, 8000])
        ax4.set_yticks(np.arange(-6000, 8001, 3500))

        #ax4.set_position([0.5,0, 0.4, 0.4])



        f.subplots_adjust(right=0.8)

        pl = ax1.scatter(MJD, velocity, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)


        cb.set_label("Mean inverse variance", fontsize=30)

        plt.show()





    def MJD_ve(self,MJD,velocity):
        print(MJD.shape, velocity.shape)


        plt.plot(MJD, velocity, "ro")
        plt.title("Radial velocity shift vs MJD",fontsize =20)
        plt.xlabel("MJD",fontsize =20)
        plt.ylabel("Radial velocity shift $m/s$",fontsize =20)
        plt.show()

    def MJD_ve_colorbar(self,MJD,velocity,mean_ivar):

        print(MJD.shape,velocity.shape,mean_ivar.shape)
        # Plot
        font = {'weight': 'bold', 'size': 25}
        matplotlib.rc('font', **font)
        fig = plt.figure()

        ax1 = fig.add_axes([0.10, 0.05, 0.85, 0.85])
        pl = ax1.scatter(MJD, velocity, marker='x', c=mean_ivar,
                         vmin=10000, vmax=40000, alpha=1)

        plt.xlabel('MJD', fontsize=30)
        plt.ylabel('RV shifts $m/s$', fontsize=30)
        fig.suptitle('RV shifts vs MJD', fontsize=24)

        #axes = plt.gca()
        # axes.set_xlim([-4,4])
        #axes.set_ylim([0, 50000])

        cb = plt.colorbar(pl, ax=ax1, orientation='horizontal')
        cb.set_label("Mean inverse variance", fontsize=30)

        plt.show()


    def Fiber_ve_colorbar(self,Fiber,velocity,mean_ivar):

        print(Fiber.shape,velocity.shape,mean_ivar.shape)
        # Plot
        font = {'weight': 'bold', 'size': 25}
        matplotlib.rc('font', **font)
        fig = plt.figure()

        ax1 = fig.add_axes([0.10, 0.05, 0.85, 0.85])
        pl = ax1.scatter(Fiber, velocity, marker='x', c=mean_ivar,
                         vmin=10000, vmax=40000, alpha=1)

        plt.xlabel('FiberID', fontsize=30)
        plt.ylabel('RV shifts $m/s$', fontsize=30)
        fig.suptitle('RV shifts vs FiberID', fontsize=24)

        #axes = plt.gca()
        # axes.set_xlim([-4,4])
        #axes.set_ylim([0, 50000])

        cb = plt.colorbar(pl, ax=ax1, orientation='horizontal')
        cb.set_label("Mean inverse variance", fontsize=30)

        plt.show()

    def RA_ve_colorbar(self,RA,velocity,mean_ivar):

        print(RA.shape,velocity.shape,mean_ivar.shape)
        # Plot
        font = {'weight': 'bold', 'size': 25}
        matplotlib.rc('font', **font)
        fig = plt.figure()

        ax1 = fig.add_axes([0.10, 0.05, 0.85, 0.85])
        pl = ax1.scatter(RA, velocity, marker='x', c=mean_ivar,
                         vmin=10000, vmax=40000, alpha=1)

        plt.xlabel('RA', fontsize=30)
        plt.ylabel('RV shifts $m/s$', fontsize=30)
        fig.suptitle('RV shifts vs RA', fontsize=24)

        #axes = plt.gca()
        # axes.set_xlim([-4,4])
        #axes.set_ylim([0, 50000])

        cb = plt.colorbar(pl, ax=ax1, orientation='horizontal')
        cb.set_label("Mean inverse variance", fontsize=30)

        plt.show()

    def DEC_ve_colorbar(self,DEC,velocity,mean_ivar):

        print(DEC.shape,velocity.shape,mean_ivar.shape)
        # Plot
        font = {'weight': 'bold', 'size': 25}
        matplotlib.rc('font', **font)
        fig = plt.figure()

        ax1 = fig.add_axes([0.10, 0.05, 0.85, 0.85])
        pl = ax1.scatter(DEC, velocity, marker='x', c=mean_ivar,
                         vmin=10000, vmax=40000, alpha=1)

        plt.xlabel('DEC', fontsize=30)
        plt.ylabel('RV shifts $m/s$', fontsize=30)
        fig.suptitle('RV shifts vs DEC', fontsize=24)

        #axes = plt.gca()
        # axes.set_xlim([-4,4])
        #axes.set_ylim([0, 50000])

        cb = plt.colorbar(pl, ax=ax1, orientation='horizontal')
        cb.set_label("Mean inverse variance", fontsize=30)

        plt.show()

    def Fiber_ve(self,fiber,velocity):
        print(fiber.shape, velocity.shape)

        plt.plot(MJD, velocity, "bo")

        plt.title("Radial velocity shift vs FiberID",fontsize =20)
        plt.xlabel("FiberID",fontsize =20)
        plt.ylabel("Radial velocity shift $m/s$",fontsize =20)

        plt.show()



    ###############
    # for parameter a


    def a_subplot_4(self,a,MJD,Fiber,RA,DEC,mean_ivar):




        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)

        fig = plt.figure()


        f, ((ax1, ax2), (ax3, ax4)) = \
            plt.subplots(2, 2)

        #ax1

        ax1.scatter(MJD, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax1.set_title('a vs MJD', fontsize=24,y=0.90)

        ax1.set_xlabel('MJD', fontsize=30)
        ax1.set_ylabel('a', fontsize=30)



        #ax1.set_position([0,0.6,0.4,0.4])


        #ax2

        ax2.scatter(RA, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax2.set_title('a vs RA', fontsize=24,y=0.90)

        ax2.set_xlabel('RA', fontsize=30)
        ax2.set_ylabel('a', fontsize=30)

        #ax2.set_position([0.5, 0.6, 0.4, 0.4])





        #ax3

        ax3.scatter(DEC, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax3.set_title('a vs DEC', fontsize=24,y=0.90)

        ax3.set_xlabel('DEC', fontsize=30)
        ax3.set_ylabel('a', fontsize=30)



        #ax3.set_position([0, 0, 0.4, 0.4])



        #ax4

        ax4.scatter(Fiber, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax4.set_title('a vs FiberID', fontsize=24,y=0.90)

        ax4.set_xlabel('FIberID', fontsize=30)
        ax4.set_ylabel('a', fontsize=30)



        #ax4.set_position([0.5,0, 0.4, 0.4])



        f.subplots_adjust(right=0.8)

        pl = ax1.scatter(MJD, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)


        cb.set_label("Mean inverse variance", fontsize=30)

        plt.show()




    ###############
    # for parameter b


    def b_subplot_4(self,b,MJD,Fiber,RA,DEC,mean_ivar):




        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)

        f, ((ax1, ax2), (ax3, ax4)) = \
            plt.subplots(2, 2)

        #ax1

        ax1.scatter(MJD, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax1.set_title('b vs MJD', fontsize=24,y=0.90)

        ax1.set_xlabel('MJD', fontsize=30)
        ax1.set_ylabel('b', fontsize=30)

        ax1.set_ylim([-3,4])
        ax1.set_yticks(np.arange(-3,4.1,1))




        #ax1.set_position([0,0.6,0.4,0.4])


        #ax2

        ax2.scatter(RA, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax2.set_title('b vs RA', fontsize=24,y=0.90)

        ax2.set_xlabel('RA', fontsize=30)
        ax2.set_ylabel('b', fontsize=30)

        ax2.set_ylim([-3,4])
        ax2.set_yticks(np.arange(-3,4.1,1))

        #ax2.set_position([0.5, 0.6, 0.4, 0.4])






        #ax3

        ax3.scatter(DEC, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax3.set_title('b vs DEC', fontsize=24,y=0.90)

        ax3.set_xlabel('DEC', fontsize=30)
        ax3.set_ylabel('b', fontsize=30)

        ax3.set_ylim([-3,4])
        ax3.set_yticks(np.arange(-3,4.1,1))



        #ax3.set_position([0, 0, 0.4, 0.4])



        #ax4

        ax4.scatter(Fiber, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax4.set_title('b vs FiberID', fontsize=24,y=0.90)

        ax4.set_xlabel('FIberID', fontsize=30)
        ax4.set_ylabel('b', fontsize=30)

        ax4.set_ylim([-3,4])
        ax4.set_yticks(np.arange(-3,4.1,1))




        #ax4.set_position([0.5,0, 0.4, 0.4])



        f.subplots_adjust(right=0.8)

        pl = ax1.scatter(MJD, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)


        cb.set_label("Mean inverse variance", fontsize=30)

        plt.show()




    ###############
    # for parameter c


    def c_subplot_4(self,c,MJD,Fiber,RA,DEC,mean_ivar):




        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)

        fig = plt.figure()


        f, ((ax1, ax2), (ax3, ax4)) = \
            plt.subplots(2, 2)

        #ax1

        ax1.scatter(MJD, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax1.set_title('c vs MJD', fontsize=24,y=0.90)

        ax1.set_xlabel('MJD', fontsize=30)
        ax1.set_ylabel('c', fontsize=30)




        #ax1.set_position([0,0.6,0.4,0.4])


        #ax2

        ax2.scatter(RA, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax2.set_title('c vs RA', fontsize=24,y=0.90)

        ax2.set_xlabel('RA', fontsize=30)
        ax2.set_ylabel('c', fontsize=30)

        #ax2.set_position([0.5, 0.6, 0.4, 0.4])





        #ax3

        ax3.scatter(DEC, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax3.set_title('c vs DEC', fontsize=24,y=0.90)

        ax3.set_xlabel('DEC', fontsize=30)
        ax3.set_ylabel('c', fontsize=30)



        #ax3.set_position([0, 0, 0.4, 0.4])



        #ax4

        ax4.scatter(Fiber, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax4.set_title('c vs FiberID', fontsize=24,y=0.90)

        ax4.set_xlabel('FIberID', fontsize=30)
        ax4.set_ylabel('c', fontsize=30)


        #ax4.set_position([0.5,0, 0.4, 0.4])



        f.subplots_adjust(right=0.8)

        pl = ax1.scatter(MJD, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)


        cb.set_label("Mean inverse variance", fontsize=30)

        plt.show()


    ###############
    # for delta-chi-squared


    def dchi_subplot_4(self,dchi,MJD,Fiber,RA,DEC,mean_ivar):




        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)

        fig = plt.figure()


        f, ((ax1, ax2), (ax3, ax4)) = \
            plt.subplots(2, 2)

        #ax1

        ax1.scatter(MJD, dchi, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax1.set_title('delta-chi-squared vs MJD', fontsize=24,y=0.90)

        ax1.set_xlabel('MJD', fontsize=30)
        ax1.set_ylabel('delta-chi-squared', fontsize=30)




        #ax1.set_position([0,0.6,0.4,0.4])


        #ax2

        ax2.scatter(RA, dchi, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax2.set_title('delta-chi-squared vs RA', fontsize=24,y=0.90)

        ax2.set_xlabel('RA', fontsize=30)
        ax2.set_ylabel('delta-chi-squared', fontsize=30)

        #ax2.set_position([0.5, 0.6, 0.4, 0.4])





        #ax3

        ax3.scatter(DEC, dchi, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax3.set_title('delta-chi-squared vs DEC', fontsize=24,y=0.90)

        ax3.set_xlabel('DEC', fontsize=30)
        ax3.set_ylabel('delta-chi-squared', fontsize=30)



        #ax3.set_position([0, 0, 0.4, 0.4])



        #ax4

        ax4.scatter(Fiber, dchi, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        ax4.set_title('delta-chi-squared vs FiberID', fontsize=24,y=0.90)

        ax4.set_xlabel('FIberID', fontsize=30)
        ax4.set_ylabel('delta-chi-squared', fontsize=30)


        #ax4.set_position([0.5,0, 0.4, 0.4])



        f.subplots_adjust(right=0.8)

        pl = ax1.scatter(MJD, dchi, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=1)

        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)


        cb.set_label("Mean inverse variance", fontsize=30)

        plt.show()



class plot_histogram():


    def histogram_abc_visits(self,a,b,c):

        parameters_500 = np.c_[a,b]
        parameters_500 = np.c_[parameters_500,c]

        colors = ['b', 'g', 'r']
        name = ["a", "b", "c"]

        # Plot histogram
        fig = plt.figure()
        plt.hist(parameters_500, bins=15, color=colors, label=name)
        plt.legend(prop={'size': 30})
        fig.suptitle('Histogram of parameters a,b and c', fontsize=30)
        plt.xlabel('values of a, b and c', fontsize=30)
        plt.ylabel('Number of stars', fontsize=30)
        plt.show()

        print("Plot histogram abc-visits")



    def histogram_delta_chi_visits(self,dchi):

        colors = ['b']
        name = ["Delta_chi_squared"]

        # Plot histogram
        fig = plt.figure()
        plt.hist(dchi, bins=15, color=colors, label=name)
        plt.legend(prop={'size': 30})
        fig.suptitle('Histogram of Delta-chi-squared', fontsize=30)
        plt.xlabel('values of Delta-chi-squared', fontsize=30)
        plt.ylabel('Number of stars', fontsize=30)
        plt.show()

        print("Plot histogram delta-chi")

    def histogram_ve_visits(self,ve):

        colors = ['b']
        name = ["RV shifts"]

        # Plot histogram
        fig = plt.figure()
        plt.hist(ve, bins=15, color=colors, label=name)
        plt.legend(prop={'size': 30})
        fig.suptitle('Histogram of RV shifts', fontsize=30)
        plt.xlabel('values of RV shifts $m/s$', fontsize=30)
        plt.ylabel('Number of stars', fontsize=30)
        plt.show()

        print("Plot histogram RVs")




class plot_continuum_pixel():

    def plot_single_star(self,flux, ivar):


        # obtain contmask



        tr_ID = "biggest_c_a"

        test_labels_all_i = ["Teff", "Logg", "Fe/H"]

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

        ## diagnostic
        # contmask is bool
        contmask = np.array(contmask)


        N = len(flux[:, 0])

        name = ["combined spectrum", "combined spectrum"]
        for j in range(0, N - 2):
            name.append("individual visit")

        plt.figure()

        for i in range(N):
            plt.plot(wl, flux[i] + 0.5 * i, label=name[i])
            plt.plot(wl, (flux[i] + 0.5 * i) * contmask, "ko", label=name[i])

            # plt.errorbar(wl,flux[i] + 0.3 * i, ecolor='k', alpha=0.02, capthick=0.2, yerr=ivar[i]**(-0.5))

        axes = plt.gca()
        axes.set_xlim([15660, 15780])


        axes.set_xlim([16160,16280])
        axes.set_ylim([0.5, 1 + 0.5 * N])
        # axes.set_yticks(np.arange(0.8,1.21,0.1))

        plt.xlabel("Wave length $\AA$", fontsize=20)
        plt.ylabel("Flux", fontsize=20)
        plt.title("The fluxes of individual visits for one star", fontsize=20)

        plt.show()


    def find_peak(self,flux):

        # return the index of peak location

        N = len(flux)

        peak = []
        for i in range(50,N-50):

            m = flux[i]
            l1 = flux[i-1]
            l2 = flux[i-2]
            l3 = flux[i-10]


            r1 = flux[i+1]
            r2 = flux[i + 2]
            r3 = flux[i + 10]


            if m<min(l1,l2,l3) and m<min(r1,r2,r3):
                peak.append(i)
            elif m>max(l1,l2,l3) and m>max(r1,r2,r3):
                peak.append(i)

        peak = np.array(peak)
        return peak

    def plot_peak_single_star(self,flux, inf_flux):



        N = len(flux[:, 0])

        name = ["combined spectrum", "combined spectrum"]
        for j in range(0, N - 2):
            name.append("individual visit")

        plt.figure(1)

        plt.subplot(2,1,1)


        for i in range(N):
            plt.plot(wl, flux[i] + 0.5 * i, label=name[i])


            # plt.errorbar(wl,flux[i] + 0.3 * i, ecolor='k', alpha=0.02, capthick=0.2, yerr=ivar[i]**(-0.5))




        ### Add vertical lines

        peak_i = self.find_peak(flux[0,:])

        n_p = len(peak_i)

        for i in range(0,n_p):
            index = peak_i[i]

            plt.plot((wl[index],wl[index]), (0.5,1 + 0.5 * N), 'k-')

        axes = plt.gca()
        #axes.set_xlim([15660, 15780])

        axes.set_xlim([15680, 15740])


        # axes.set_xlim([16160,16280])
        axes.set_ylim([0.5, 1 + 0.5 * N])
        # axes.set_yticks(np.arange(0.8,1.21,0.1))

        # plt.xlabel("Wave length $\AA$", fontsize=20)
        plt.ylabel("Flux", fontsize=20)
        plt.title("The fluxes of individual visits for one star from APOGGE team", fontsize=20)


        # plot inf_flux:

        plt.subplot(2, 1, 2)



        for i in range(N):
            plt.plot(wl, inf_flux[i] + 0.5 * i, label=name[i])


            # plt.errorbar(wl,flux[i] + 0.3 * i, ecolor='k', alpha=0.02, capthick=0.2, yerr=ivar[i]**(-0.5))

        # find peak:


        ### Add vertical lines

        peak_i = self.find_peak(inf_flux[0,:])

        n_p = len(peak_i)

        for i in range(0,n_p):
            index = peak_i[i]
            plt.plot((wl[index], wl[index]), (0.5, 1 + 0.5 * N), 'k-')




        axes = plt.gca()
        #axes.set_xlim([15660, 15780])
        axes.set_xlim([15680, 15740])


        # axes.set_xlim([16160,16280])
        axes.set_ylim([0.5, 1 + 0.5 * N])
        # axes.set_yticks(np.arange(0.8,1.21,0.1))

        plt.xlabel("Wave length $\AA$", fontsize=20)
        plt.ylabel("Flux", fontsize=20)
        plt.title("The fluxes of individual visits from the Cannon", fontsize=20)

        plt.show()

model = plot()


#model.MJD_ve_colorbar(MJD,velocity,meanivar)

#model.Fiber_ve_colorbar(fiber_id,velocity,meanivar)

#model.RA_ve_colorbar(RA,velocity,meanivar)

#model.DEC_ve_colorbar(DEC,velocity,meanivar)

#model.ve_subplot_4(velocity=velocity,MJD=MJD,Fiber=fiber_id,RA=RA,DEC=DEC,mean_ivar=meanivar)

# print(parameters.shape,MJD.shape)

#model.a_subplot_4(a=parameters,MJD=MJD,Fiber=fiber_id,RA=RA,DEC=DEC,mean_ivar=meanivar)


#model.c_subplot_4(c=parameters,MJD=MJD,Fiber=fiber_id,RA=RA,DEC=DEC,mean_ivar=meanivar)


"""
plt.plot(MJD,parameters,"ro")

axes = plt.gca()

axes.set_ylim([-5,5])
axes.set_yticks(np.arange(-5,5.1,1))

plt.show()

"""


# model.dchi_subplot_4(dchi=dchi,MJD=MJD,Fiber=fiber_id,RA=RA,DEC=DEC,mean_ivar=meanivar)


#model2 = plot_histogram()

#model2.histogram_abc_visits(a,b,c)

#model2.histogram_delta_chi_visits(dchi)
#model2.histogram_ve_visits(velocity)




p = 628

star_choose = fits.open(path_fits[p])

nor_i = star_choose[0].data
ivar_i = star_choose[1].data
inf_i = star_choose[2].data

print(nor_i.shape,ivar_i.shape)

label_i = star_choose[8].data
SNR = star_choose[0].header["SNR"]

print(SNR,label_i)


model3 = plot_continuum_pixel()

#model3.plot_single_star(flux = nor_i,ivar = ivar_i)

model3.plot_peak_single_star(flux = nor_i,inf_flux=inf_i)

