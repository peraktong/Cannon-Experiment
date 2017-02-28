
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


# here parameter_sim represent the results from simultaneously fitting


class plot():

    def read_data(self):

        N = len(path_fits)
        print(N)

        velocity = []
        velocity_new = []
        fiber_id = []
        mean_ivar = []
        parameters = np.array([0,1,0])
        parameters_new = np.array([0,1,0])
        parameters_sim = np.array([0,1,0])



        inf_label = []
        dchi = []

        MJD = []
        HJD = []
        meanivar = []

        RA = []
        DEC = []
        SNR = []

        airmass = []

        # star name and the number of visit
        # dimension N*2
        star_visit = []
        star_name = []



        for i in range(0, N):
            print("loading star %d" % (i + 1))

            star_name_i = path_fits[i]

            star_i = fits.open(path_fits[i])

            ni = len(star_i[4].data[:, 0])



            # mean ivar
            one = np.ones(ni - 2)
            for si in range(0,ni-2):
                star_name = np.append(star_name, star_name_i)
                star_visit.append(si)

            meanivar = np.append(meanivar, one * mi[i])



            dchi = np.append(dchi, star_i[6].data[2:ni])
            # SNR RA DEC
            SNR = np.append(SNR, (star_i[0].header["SNR"] * one))
            RA = np.append(RA, (star_i[0].header["RA"] * one))
            DEC = np.append(DEC, (star_i[0].header["DEC"] * one))

            velocity = np.append(velocity, star_i[10].data[2:ni, 0])
            velocity_new = np.append(velocity_new, star_i[15].data[2:ni, 0])
            fiber_id = np.append(fiber_id, star_i[7].data)
            mean_ivar.append(np.mean(star_i[1].data[0]))
            parameters = np.vstack((parameters,star_i[4].data[2:ni,0:3]))
            parameters_new = np.vstack((parameters_new, star_i[14].data[2:ni, 0:3]))
            parameters_sim = np.vstack((parameters_sim,star_i[18].data[2:ni, 3:6]))


            MJD = np.append(MJD, star_i[11].data)
            HJD = np.append(HJD,star_i[16].data)
            airmass = np.append(airmass,star_i[17].data)

            print(star_i[4].data[:, 0].shape, star_i[0].data.shape, star_i[11].data.shape, star_i[12].data.shape)
            print(star_i[11].data)

        self.path_fits = path_fits
        velocity = np.array(velocity)
        self.velocity = velocity

        velocity_new = np.array(velocity_new)
        self.velocity_new = velocity_new

        fiber_id = np.array(fiber_id)
        self.fiber_id = fiber_id

        na = len(parameters[:,0])
        parameters = parameters[1:na,:]
        self.parameters = parameters

        parameters_new = parameters_new[1:na,:]
        self.parameters_new = parameters_new

        parameters_sim = parameters_sim[1:na,:]
        self.parameters_sim = parameters_sim



        dchi = np.array(dchi)
        self.dchi = dchi


        MJD = np.array(MJD)
        self.MJD = MJD

        HJD = np.array(HJD)
        self.HJD = HJD

        airmass = np.array(airmass)
        self.airmass = airmass

        SNR = np.array(SNR)
        self.SNR = SNR

        RA = np.array(RA)
        self.RA = RA

        DEC = np.array(DEC)
        self.DEC = DEC

        meanivar = np.array(meanivar)
        self.meanivar = meanivar

        star_name =np.array(star_name)
        self.star_name = star_name
        star_visit = np.array(star_visit)
        self.star_visit = star_visit

        print("star name shape")
        print(star_name.shape,star_visit.shape)



        # check shape

        print(MJD.shape, SNR.shape, RA.shape, DEC.shape)

        # give values:

    ## Choose single stars:


    def histogram_2_2_rv_abc(self):

        a = self.parameters_sim[:,0]
        b = self.parameters_sim[:,1]
        c = self.parameters_sim[:,2]
        RV = (c-a)/(a+b+c)*4144.68





        font = {'weight': 'bold', 'size': 15}
        matplotlib.rc('font', **font)

        fig = plt.figure()



        f, ((ax1, ax2), (ax3, ax4)) = \
            plt.subplots(2, 2)

        colors = ["cyan",'b', 'g', 'r']
        name = ["RV","a", "b", "c"]

        # histogram of rv
        #ax1

        rms_RV = (np.sum(RV*RV)/len(RV))**0.5
        rms_a = (np.sum(a * a) / len(a)) ** 0.5
        rms_b = (np.sum(b*b) / len(b)) ** 0.5
        rms_c = (np.sum(c * c) / len(c)) ** 0.5

        ax1.hist(RV, bins=40, color=colors[0], label="%s RMS = %.2f"%(name[0],rms_RV))


        #ax1.set_title('Histogram of Radial velocity shifts', fontsize=30)
        ax1.set_xlabel('values of radial velocity shifts $m/s$', fontsize=15)
        ax1.set_ylabel('Number', fontsize=15)
        ax1.legend(prop={'size': 15})

        # add vertical grey line
        # ax1.plot((wl[index], wl[index]), (0.5, 1 + 0.5 * N), 'k-', linewidth=1.5)


        # histogram of a
        #ax2

        ax2.hist(a, bins=40, color=colors[1], label="%s RMS = %.2f"%(name[1],rms_a))


        #ax2.set_title('Histogram of parameter a', fontsize=30)
        ax2.set_xlabel('values of parameter a', fontsize=15)
        ax2.set_ylabel('Number', fontsize=15)
        ax2.legend(prop={'size': 15})

        # add vertical grey line
        # ax1.plot((wl[index], wl[index]), (0.5, 1 + 0.5 * N), 'k-', linewidth=1.5)


        # histogram of b
        #ax3

        ax3.hist(b, bins=40, color=colors[2], label="%s RMS = %.2f"%(name[2],rms_b))
        ax3.legend(prop={'size': 15})


        #ax3.set_title('Histogram of paramete b', fontsize=30)
        ax3.set_xlabel("values of parameter b", fontsize=15)
        ax3.set_ylabel('Number', fontsize=15)

        # add vertical grey line
        # ax1.plot((wl[index], wl[index]), (0.5, 1 + 0.5 * N), 'k-', linewidth=1.5)



        # histogram of c
        #ax4

        ax4.hist(c, bins=40, color=colors[3], label="%s RMS = %.2f"%(name[3],rms_c))
        ax4.legend(prop={'size': 15})


        #ax4.set_title('Histogram of parameter c', fontsize=30)
        ax4.set_xlabel("values of parameter c", fontsize=15)
        ax4.set_ylabel('Number', fontsize=15)

        # add vertical grey line
        # ax1.plot((wl[index], wl[index]), (0.5, 1 + 0.5 * N), 'k-', linewidth=1.5)
        f.suptitle("Histogram of RV shifts, a, b and c by simultaneous fitting labels and abc")
        #f.suptitle("Histogram of RV shifts, a, b and c by using the absorption lines")




        plt.show()




    # RV vs HJD RA DEC Fiber Airmass new
    def ve_sim_subplot_5(self,HJD,Fiber,RA,DEC,airmass,mean_ivar,SNR):

        a = self.parameters_sim[:,0]
        b = self.parameters_sim[:,1]
        c = self.parameters_sim[:,2]


        velocity = (c-a)/(a+b+c)*4144.68

        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)

        fig = plt.figure()


        f, ((ax1, ax2,ax3), (ax4, ax5,ax6)) = \
            plt.subplots(2, 3)


        alpha = 0.3
        #ax1

        ax1.scatter(HJD, velocity, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax1.set_title('RV shifts vs HJD', fontsize=24,y=0.85)

        ax1.set_xlabel('HJD', fontsize=20)
        ax1.set_ylabel('RV shifts $m/s$', fontsize=20)

        # add vertical line:
        #ax1.plot((np.min(HJD),np.max(HJD)), (0,0), 'k-', linewidth=1)
        ax1.axhline(y=0, linewidth=1, color="k", alpha=0.5)



        ax1.set_ylim([-6000, 8000])
        ax1.set_yticks(np.arange(-6000, 8001, 3500))


        #ax1.set_position([0,0.6,0.4,0.4])


        #ax2

        ax2.scatter(RA, velocity, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax2.set_title('RV shifts vs RA', fontsize=24,y=0.85)
        # add vertical line:
        #ax2.plot((np.min(RA),np.max(RA)), (0,0), 'k-', linewidth=1)
        ax2.axhline(y=0, linewidth=1, color="k", alpha=0.5)



        ax2.set_xlabel('RA', fontsize=20)
        #ax2.set_ylabel('RV shifts $m/s$', fontsize=20)

        #ax2.set_position([0.5, 0.6, 0.4, 0.4])

        ax2.set_ylim([-6000, 8000])
        ax2.set_yticks(np.arange(-6000, 8001, 3500))



        #ax3

        ax3.scatter(DEC, velocity, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax3.set_title('RV shifts vs DEC', fontsize=24,y=0.85)

        # add vertical line:
        #ax3.plot((np.min(DEC),np.max(DEC)), (0,0), 'k-', linewidth=1)
        ax3.axhline(y=0, linewidth=1, color="k", alpha=0.5)


        ax3.set_xlabel('DEC', fontsize=20)
        #ax3.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax3.set_ylim([-6000, 8000])
        ax3.set_yticks(np.arange(-6000, 8001, 3500))



        #ax4

        ax4.scatter(Fiber, velocity, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax4.set_title('RV shifts vs FiberID', fontsize=24,y=0.85)

        # add vertical line:
        ax4.axhline(y=0, linewidth=1, color="k", alpha=0.5)
        #ax4.plot((np.min(Fiber),np.max(Fiber)), (0,0), 'k-', linewidth=1)



        ax4.set_xlabel('FIberID', fontsize=20)
        ax4.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax4.set_ylim([-6000, 8000])
        ax4.set_yticks(np.arange(-6000, 8001, 3500))


        #ax5

        ax5.scatter(airmass, velocity, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax5.set_title('RV shifts vs air mass', fontsize=24,y=0.85)

        # add vertical line:
        ax5.axhline(y=0, linewidth=1, color="k", alpha=0.5)
        #ax5.plot((np.min(airmass),np.max(airmass)), (0,0), 'k-', linewidth=1)



        ax5.set_xlabel('FIberID', fontsize=20)
        #ax5.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax5.set_ylim([-6000, 8000])
        ax5.set_yticks(np.arange(-6000, 8001, 3500))

        #ax6

        ax6.scatter(SNR, velocity, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax6.set_title('RV shifts vs SNR', fontsize=24,y=0.85)

        # add vertical line:
        ax6.axhline(y=0, linewidth=1,color="k",alpha=0.5)
        #ax6.plot((np.min(SNR),np.max(SNR)), (0,0), )



        ax6.set_xlabel('SNR', fontsize=20)

        ax6.set_ylim([-6000, 8000])
        ax6.set_yticks(np.arange(-6000, 8001, 3500))

        f.subplots_adjust(right=0.8)

        pl = ax1.scatter(HJD, velocity, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)


        cb.set_label("Mean inverse variance", fontsize=20)
        f.suptitle("RV shifts from simultaneous fitting vs HJD, RA, DEC, FiberID, Airmass and SNR", fontsize=30)

        plt.show()

## a vs them

    # RV vs HJD RA DEC Fiber Airmass
    def a_subplot_5(self,a,HJD,Fiber,RA,DEC,airmass,mean_ivar,SNR):

        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)

        fig = plt.figure()


        f, ((ax1, ax2,ax3), (ax4, ax5,ax6)) = \
            plt.subplots(2, 3)


        alpha = 0.3
        #ax1

        ax1.scatter(HJD, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax1.set_title('a vs HJD', fontsize=24,y=0.85)

        ax1.set_xlabel('HJD', fontsize=20)
        ax1.set_ylabel('a', fontsize=20)

        # add vertical line:
        #ax1.plot((np.min(HJD),np.max(HJD)), (0,0), 'k-', linewidth=1)
        ax1.axhline(y=0, linewidth=1, color="k", alpha=0.5)


        ax1.set_ylim([-2,3])
        ax1.set_yticks(np.arange(-2,3.1,1))

        #ax1.set_position([0,0.6,0.4,0.4])


        #ax2

        ax2.scatter(RA, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax2.set_title('a vs RA', fontsize=24,y=0.85)
        # add vertical line:
        #ax2.plot((np.min(RA),np.max(RA)), (0,0), 'k-', linewidth=1)
        ax2.axhline(y=0, linewidth=1, color="k", alpha=0.5)



        ax2.set_xlabel('RA', fontsize=20)
        #ax2.set_ylabel('a', fontsize=20)

        #ax2.set_position([0.5, 0.6, 0.4, 0.4])

        ax2.set_ylim([-2,3])
        ax2.set_yticks(np.arange(-2,3.1,1))


        #ax3

        ax3.scatter(DEC, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax3.set_title('a vs DEC', fontsize=24,y=0.85)

        # add vertical line:
        #ax3.plot((np.min(DEC),np.max(DEC)), (0,0), 'k-', linewidth=1)
        ax3.axhline(y=0, linewidth=1, color="k", alpha=0.5)


        ax3.set_xlabel('DEC', fontsize=20)
        #ax3.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax3.set_ylim([-2,3])
        ax3.set_yticks(np.arange(-2,3.1,1))

        #ax3.set_position([0, 0, 0.4, 0.4])



        #ax4

        ax4.scatter(Fiber, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax4.set_title('a vs FiberID', fontsize=24,y=0.85)

        # add vertical line:
        ax4.axhline(y=0, linewidth=1, color="k", alpha=0.5)
        #ax4.plot((np.min(Fiber),np.max(Fiber)), (0,0), 'k-', linewidth=1)



        ax4.set_xlabel('FIberID', fontsize=20)
        ax4.set_ylabel('a', fontsize=20)

        ax4.set_ylim([-2,3])
        ax4.set_yticks(np.arange(-2,3.1,1))

        #ax4.set_position([0.5,0, 0.4, 0.4])


        #ax5

        ax5.scatter(airmass, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax5.set_title('a vs air mass', fontsize=24,y=0.85)

        # add vertical line:
        ax5.axhline(y=0, linewidth=1, color="k", alpha=0.5)
        #ax5.plot((np.min(airmass),np.max(airmass)), (0,0), 'k-', linewidth=1)



        ax5.set_xlabel('FIberID', fontsize=20)
        #ax5.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax5.set_ylim([-2,3])
        ax5.set_yticks(np.arange(-2,3.1,1))


        #ax6

        ax6.scatter(SNR, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax6.set_title('a vs SNR', fontsize=24,y=0.85)

        # add vertical line:
        ax6.axhline(y=0, linewidth=1,color="k",alpha=0.5)
        #ax6.plot((np.min(SNR),np.max(SNR)), (0,0), )



        ax6.set_xlabel('SNR', fontsize=20)
        #ax6.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax6.set_ylim([-2,3])
        ax6.set_yticks(np.arange(-2,3.1,1))


        f.subplots_adjust(right=0.8)

        pl = ax1.scatter(HJD, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)


        cb.set_label("Mean inverse variance", fontsize=20)
        f.suptitle("Parameter a from simultaneous fitting vs HJD, RA, DEC, FiberID, Airmass and SNR", fontsize=30)

        plt.show()





    # RV vs HJD RA DEC Fiber Airmass new

    # RV vs HJD RA DEC Fiber Airmass
    def a_new_subplot_5(self,a,HJD,Fiber,RA,DEC,airmass,mean_ivar,SNR):

        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)

        fig = plt.figure()


        f, ((ax1, ax2,ax3), (ax4, ax5,ax6)) = \
            plt.subplots(2, 3)


        alpha = 0.3
        #ax1

        ax1.scatter(HJD, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax1.set_title('a vs HJD', fontsize=24,y=0.85)

        ax1.set_xlabel('HJD', fontsize=20)
        ax1.set_ylabel('a', fontsize=20)

        # add vertical line:
        #ax1.plot((np.min(HJD),np.max(HJD)), (0,0), 'k-', linewidth=1)
        ax1.axhline(y=0, linewidth=1, color="k", alpha=0.5)


        ax1.set_ylim([-2,3])
        ax1.set_yticks(np.arange(-2,3.1,1))

        #ax1.set_position([0,0.6,0.4,0.4])


        #ax2

        ax2.scatter(RA, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax2.set_title('a vs RA', fontsize=24,y=0.85)
        # add vertical line:
        #ax2.plot((np.min(RA),np.max(RA)), (0,0), 'k-', linewidth=1)
        ax2.axhline(y=0, linewidth=1, color="k", alpha=0.5)



        ax2.set_xlabel('RA', fontsize=20)
        #ax2.set_ylabel('a', fontsize=20)

        #ax2.set_position([0.5, 0.6, 0.4, 0.4])

        ax2.set_ylim([-2,3])
        ax2.set_yticks(np.arange(-2,3.1,1))


        #ax3

        ax3.scatter(DEC, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax3.set_title('a vs DEC', fontsize=24,y=0.85)

        # add vertical line:
        #ax3.plot((np.min(DEC),np.max(DEC)), (0,0), 'k-', linewidth=1)
        ax3.axhline(y=0, linewidth=1, color="k", alpha=0.5)


        ax3.set_xlabel('DEC', fontsize=20)
        #ax3.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax3.set_ylim([-2,3])
        ax3.set_yticks(np.arange(-2,3.1,1))

        #ax3.set_position([0, 0, 0.4, 0.4])



        #ax4

        ax4.scatter(Fiber, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax4.set_title('a vs FiberID', fontsize=24,y=0.85)

        # add vertical line:
        ax4.axhline(y=0, linewidth=1, color="k", alpha=0.5)
        #ax4.plot((np.min(Fiber),np.max(Fiber)), (0,0), 'k-', linewidth=1)



        ax4.set_xlabel('FIberID', fontsize=20)
        ax4.set_ylabel('a', fontsize=20)

        ax4.set_ylim([-2,3])
        ax4.set_yticks(np.arange(-2,3.1,1))

        #ax4.set_position([0.5,0, 0.4, 0.4])


        #ax5

        ax5.scatter(airmass, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax5.set_title('a vs air mass', fontsize=24,y=0.85)

        # add vertical line:
        ax5.axhline(y=0, linewidth=1, color="k", alpha=0.5)
        #ax5.plot((np.min(airmass),np.max(airmass)), (0,0), 'k-', linewidth=1)



        ax5.set_xlabel('FIberID', fontsize=20)
        #ax5.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax5.set_ylim([-2,3])
        ax5.set_yticks(np.arange(-2,3.1,1))


        #ax6

        ax6.scatter(SNR, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax6.set_title('a vs SNR', fontsize=24,y=0.85)

        # add vertical line:
        ax6.axhline(y=0, linewidth=1,color="k",alpha=0.5)
        #ax6.plot((np.min(SNR),np.max(SNR)), (0,0), )



        ax6.set_xlabel('SNR', fontsize=20)
        #ax6.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax6.set_ylim([-2,3])
        ax6.set_yticks(np.arange(-2,3.1,1))


        f.subplots_adjust(right=0.8)

        pl = ax1.scatter(HJD, a, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)


        cb.set_label("Mean inverse variance", fontsize=20)
        f.suptitle("Parameter a from the absorption line vs HJD, RA, DEC, FiberID, Airmass and SNR", fontsize=30)

        plt.show()

        ## a vs them

        # RV vs HJD RA DEC Fiber Airmass
        def a_subplot_5(self, a, HJD, Fiber, RA, DEC, airmass, mean_ivar, SNR):
            font = {'weight': 'bold', 'size': 13}
            matplotlib.rc('font', **font)

            fig = plt.figure()

            f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = \
                plt.subplots(2, 3)

            alpha = 0.3
            # ax1

            ax1.scatter(HJD, a, marker='x', c=mean_ivar,
                        vmin=10000, vmax=40000, alpha=alpha)

            ax1.set_title('a vs HJD', fontsize=24, y=0.85)

            ax1.set_xlabel('HJD', fontsize=20)
            ax1.set_ylabel('a', fontsize=20)

            # add vertical line:
            # ax1.plot((np.min(HJD),np.max(HJD)), (0,0), 'k-', linewidth=1)
            ax1.axhline(y=0, linewidth=1, color="k", alpha=0.5)

            ax1.set_ylim([-2, 3])
            ax1.set_yticks(np.arange(-2, 3.1, 1))

            # ax1.set_position([0,0.6,0.4,0.4])


            # ax2

            ax2.scatter(RA, a, marker='x', c=mean_ivar,
                        vmin=10000, vmax=40000, alpha=alpha)

            ax2.set_title('a vs RA', fontsize=24, y=0.85)
            # add vertical line:
            # ax2.plot((np.min(RA),np.max(RA)), (0,0), 'k-', linewidth=1)
            ax2.axhline(y=0, linewidth=1, color="k", alpha=0.5)

            ax2.set_xlabel('RA', fontsize=20)
            # ax2.set_ylabel('a', fontsize=20)

            # ax2.set_position([0.5, 0.6, 0.4, 0.4])

            ax2.set_ylim([-2, 3])
            ax2.set_yticks(np.arange(-2, 3.1, 1))

            # ax3

            ax3.scatter(DEC, a, marker='x', c=mean_ivar,
                        vmin=10000, vmax=40000, alpha=alpha)

            ax3.set_title('a vs DEC', fontsize=24, y=0.85)

            # add vertical line:
            # ax3.plot((np.min(DEC),np.max(DEC)), (0,0), 'k-', linewidth=1)
            ax3.axhline(y=0, linewidth=1, color="k", alpha=0.5)

            ax3.set_xlabel('DEC', fontsize=20)
            # ax3.set_ylabel('RV shifts $m/s$', fontsize=20)

            ax3.set_ylim([-2, 3])
            ax3.set_yticks(np.arange(-2, 3.1, 1))

            # ax3.set_position([0, 0, 0.4, 0.4])



            # ax4

            ax4.scatter(Fiber, a, marker='x', c=mean_ivar,
                        vmin=10000, vmax=40000, alpha=alpha)

            ax4.set_title('a vs FiberID', fontsize=24, y=0.85)

            # add vertical line:
            ax4.axhline(y=0, linewidth=1, color="k", alpha=0.5)
            # ax4.plot((np.min(Fiber),np.max(Fiber)), (0,0), 'k-', linewidth=1)



            ax4.set_xlabel('FIberID', fontsize=20)
            ax4.set_ylabel('a', fontsize=20)

            ax4.set_ylim([-2, 3])
            ax4.set_yticks(np.arange(-2, 3.1, 1))

            # ax4.set_position([0.5,0, 0.4, 0.4])


            # ax5

            ax5.scatter(airmass, a, marker='x', c=mean_ivar,
                        vmin=10000, vmax=40000, alpha=alpha)

            ax5.set_title('a vs air mass', fontsize=24, y=0.85)

            # add vertical line:
            ax5.axhline(y=0, linewidth=1, color="k", alpha=0.5)
            # ax5.plot((np.min(airmass),np.max(airmass)), (0,0), 'k-', linewidth=1)



            ax5.set_xlabel('FIberID', fontsize=20)
            # ax5.set_ylabel('RV shifts $m/s$', fontsize=20)

            ax5.set_ylim([-2, 3])
            ax5.set_yticks(np.arange(-2, 3.1, 1))

            # ax6

            ax6.scatter(SNR, a, marker='x', c=mean_ivar,
                        vmin=10000, vmax=40000, alpha=alpha)

            ax6.set_title('a vs SNR', fontsize=24, y=0.85)

            # add vertical line:
            ax6.axhline(y=0, linewidth=1, color="k", alpha=0.5)
            # ax6.plot((np.min(SNR),np.max(SNR)), (0,0), )



            ax6.set_xlabel('SNR', fontsize=20)
            # ax6.set_ylabel('RV shifts $m/s$', fontsize=20)

            ax6.set_ylim([-2, 3])
            ax6.set_yticks(np.arange(-2, 3.1, 1))

            f.subplots_adjust(right=0.8)

            pl = ax1.scatter(HJD, a, marker='x', c=mean_ivar,
                             vmin=10000, vmax=40000, alpha=alpha)

            cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
            cb = f.colorbar(pl, cax=cbar_ax)

            cb.set_label("Mean inverse variance", fontsize=20)
            f.suptitle("Parameter a from the whole spectrum vs HJD, RA, DEC, FiberID, Airmass and SNR", fontsize=30)

            plt.show()

        # RV vs HJD RA DEC Fiber Airmass new

        # RV vs HJD RA DEC Fiber Airmass
        def a_new_subplot_5(self, a, HJD, Fiber, RA, DEC, airmass, mean_ivar, SNR):
            font = {'weight': 'bold', 'size': 13}
            matplotlib.rc('font', **font)

            fig = plt.figure()

            f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = \
                plt.subplots(2, 3)

            alpha = 0.3
            # ax1

            ax1.scatter(HJD, a, marker='x', c=mean_ivar,
                        vmin=10000, vmax=40000, alpha=alpha)

            ax1.set_title('a vs HJD', fontsize=24, y=0.85)

            ax1.set_xlabel('HJD', fontsize=20)
            ax1.set_ylabel('a', fontsize=20)

            # add vertical line:
            # ax1.plot((np.min(HJD),np.max(HJD)), (0,0), 'k-', linewidth=1)
            ax1.axhline(y=0, linewidth=1, color="k", alpha=0.5)

            ax1.set_ylim([-2, 3])
            ax1.set_yticks(np.arange(-2, 3.1, 1))

            # ax1.set_position([0,0.6,0.4,0.4])


            # ax2

            ax2.scatter(RA, a, marker='x', c=mean_ivar,
                        vmin=10000, vmax=40000, alpha=alpha)

            ax2.set_title('a vs RA', fontsize=24, y=0.85)
            # add vertical line:
            # ax2.plot((np.min(RA),np.max(RA)), (0,0), 'k-', linewidth=1)
            ax2.axhline(y=0, linewidth=1, color="k", alpha=0.5)

            ax2.set_xlabel('RA', fontsize=20)
            # ax2.set_ylabel('a', fontsize=20)

            # ax2.set_position([0.5, 0.6, 0.4, 0.4])

            ax2.set_ylim([-2, 3])
            ax2.set_yticks(np.arange(-2, 3.1, 1))

            # ax3

            ax3.scatter(DEC, a, marker='x', c=mean_ivar,
                        vmin=10000, vmax=40000, alpha=alpha)

            ax3.set_title('a vs DEC', fontsize=24, y=0.85)

            # add vertical line:
            # ax3.plot((np.min(DEC),np.max(DEC)), (0,0), 'k-', linewidth=1)
            ax3.axhline(y=0, linewidth=1, color="k", alpha=0.5)

            ax3.set_xlabel('DEC', fontsize=20)
            # ax3.set_ylabel('RV shifts $m/s$', fontsize=20)

            ax3.set_ylim([-2, 3])
            ax3.set_yticks(np.arange(-2, 3.1, 1))

            # ax3.set_position([0, 0, 0.4, 0.4])



            # ax4

            ax4.scatter(Fiber, a, marker='x', c=mean_ivar,
                        vmin=10000, vmax=40000, alpha=alpha)

            ax4.set_title('a vs FiberID', fontsize=24, y=0.85)

            # add vertical line:
            ax4.axhline(y=0, linewidth=1, color="k", alpha=0.5)
            # ax4.plot((np.min(Fiber),np.max(Fiber)), (0,0), 'k-', linewidth=1)



            ax4.set_xlabel('FIberID', fontsize=20)
            ax4.set_ylabel('a', fontsize=20)

            ax4.set_ylim([-2, 3])
            ax4.set_yticks(np.arange(-2, 3.1, 1))

            # ax4.set_position([0.5,0, 0.4, 0.4])


            # ax5

            ax5.scatter(airmass, a, marker='x', c=mean_ivar,
                        vmin=10000, vmax=40000, alpha=alpha)

            ax5.set_title('a vs air mass', fontsize=24, y=0.85)

            # add vertical line:
            ax5.axhline(y=0, linewidth=1, color="k", alpha=0.5)
            # ax5.plot((np.min(airmass),np.max(airmass)), (0,0), 'k-', linewidth=1)



            ax5.set_xlabel('FIberID', fontsize=20)
            # ax5.set_ylabel('RV shifts $m/s$', fontsize=20)

            ax5.set_ylim([-2, 3])
            ax5.set_yticks(np.arange(-2, 3.1, 1))

            # ax6

            ax6.scatter(SNR, a, marker='x', c=mean_ivar,
                        vmin=10000, vmax=40000, alpha=alpha)

            ax6.set_title('a vs SNR', fontsize=24, y=0.85)

            # add vertical line:
            ax6.axhline(y=0, linewidth=1, color="k", alpha=0.5)
            # ax6.plot((np.min(SNR),np.max(SNR)), (0,0), )



            ax6.set_xlabel('SNR', fontsize=20)
            # ax6.set_ylabel('RV shifts $m/s$', fontsize=20)

            ax6.set_ylim([-2, 3])
            ax6.set_yticks(np.arange(-2, 3.1, 1))

            f.subplots_adjust(right=0.8)

            pl = ax1.scatter(HJD, a, marker='x', c=mean_ivar,
                             vmin=10000, vmax=40000, alpha=alpha)

            cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
            cb = f.colorbar(pl, cax=cbar_ax)

            cb.set_label("Mean inverse variance", fontsize=20)
            f.suptitle("Parameter a from the absorption line vs HJD, RA, DEC, FiberID, Airmass and SNR", fontsize=30)

            plt.show()


## b vs them

    # RV vs HJD RA DEC Fiber Airmass
    def b_subplot_5(self,b,HJD,Fiber,RA,DEC,airmass,mean_ivar,SNR):

        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)

        fig = plt.figure()


        f, ((ax1, ax2,ax3), (ax4, ax5,ax6)) = \
            plt.subplots(2, 3)


        alpha = 0.3
        #ax1

        ax1.scatter(HJD, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax1.set_title('b vs HJD', fontsize=24,y=0.85)

        ax1.set_xlabel('HJD', fontsize=20)
        ax1.set_ylabel('b', fontsize=20)

        # add vertical line:
        #ax1.plot((np.min(HJD),np.max(HJD)), (0,0), 'k-', linewidth=1)
        ax1.axhline(y=0, linewidth=1, color="k", alpha=0.5)


        ax1.set_ylim([-3,4])
        ax1.set_yticks(np.arange(-3,4.1,1))

        #ax1.set_position([0,0.6,0.4,0.4])


        #ax2

        ax2.scatter(RA, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax2.set_title('b vs RA', fontsize=24,y=0.85)
        # add vertical line:
        #ax2.plot((np.min(RA),np.max(RA)), (0,0), 'k-', linewidth=1)
        ax2.axhline(y=0, linewidth=1, color="k", alpha=0.5)



        ax2.set_xlabel('RA', fontsize=20)
        #ax2.set_ylabel('a', fontsize=20)

        #ax2.set_position([0.5, 0.6, 0.4, 0.4])

        ax2.set_ylim([-3,4])
        ax2.set_yticks(np.arange(-3,4.1,1))


        #ax3

        ax3.scatter(DEC, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax3.set_title('b vs DEC', fontsize=24,y=0.85)

        # add vertical line:
        #ax3.plot((np.min(DEC),np.max(DEC)), (0,0), 'k-', linewidth=1)
        ax3.axhline(y=0, linewidth=1, color="k", alpha=0.5)


        ax3.set_xlabel('DEC', fontsize=20)
        #ax3.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax3.set_ylim([-3,4])
        ax3.set_yticks(np.arange(-3,4.1,1))

        #ax3.set_position([0, 0, 0.4, 0.4])



        #ax4

        ax4.scatter(Fiber, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax4.set_title('b vs FiberID', fontsize=24,y=0.85)

        # add vertical line:
        ax4.axhline(y=0, linewidth=1, color="k", alpha=0.5)
        #ax4.plot((np.min(Fiber),np.max(Fiber)), (0,0), 'k-', linewidth=1)



        ax4.set_xlabel('FIberID', fontsize=20)
        ax4.set_ylabel('b', fontsize=20)

        ax4.set_ylim([-3,4])
        ax4.set_yticks(np.arange(-3,4.1,1))

        #ax4.set_position([0.5,0, 0.4, 0.4])


        #ax5

        ax5.scatter(airmass, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax5.set_title('b vs air mass', fontsize=24,y=0.85)

        # add vertical line:
        ax5.axhline(y=0, linewidth=1, color="k", alpha=0.5)
        #ax5.plot((np.min(airmass),np.max(airmass)), (0,0), 'k-', linewidth=1)



        ax5.set_xlabel('FIberID', fontsize=20)
        #ax5.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax5.set_ylim([-3,4])
        ax5.set_yticks(np.arange(-3,4.1,1))


        #ax6

        ax6.scatter(SNR, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax6.set_title('b vs SNR', fontsize=24,y=0.85)

        # add vertical line:
        ax6.axhline(y=0, linewidth=1,color="k",alpha=0.5)
        #ax6.plot((np.min(SNR),np.max(SNR)), (0,0), )



        ax6.set_xlabel('SNR', fontsize=20)
        #ax6.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax6.set_ylim([-3,4])
        ax6.set_yticks(np.arange(-3,4.1,1))


        f.subplots_adjust(right=0.8)

        pl = ax1.scatter(HJD, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)


        cb.set_label("Mean inverse variance", fontsize=20)
        f.suptitle("Parameter b from the whole spectrum vs HJD, RA, DEC, FiberID, Airmass and SNR", fontsize=30)

        plt.show()





    # RV vs HJD RA DEC Fiber Airmass
    def b_new_subplot_5(self,b,HJD,Fiber,RA,DEC,airmass,mean_ivar,SNR):

        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)

        fig = plt.figure()


        f, ((ax1, ax2,ax3), (ax4, ax5,ax6)) = \
            plt.subplots(2, 3)


        alpha = 0.3
        #ax1

        ax1.scatter(HJD, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax1.set_title('b vs HJD', fontsize=24,y=0.85)

        ax1.set_xlabel('HJD', fontsize=20)
        ax1.set_ylabel('b', fontsize=20)

        # add vertical line:
        #ax1.plot((np.min(HJD),np.max(HJD)), (0,0), 'k-', linewidth=1)
        ax1.axhline(y=0, linewidth=1, color="k", alpha=0.5)


        ax1.set_ylim([-3,4])
        ax1.set_yticks(np.arange(-3,4.1,1))

        #ax1.set_position([0,0.6,0.4,0.4])


        #ax2

        ax2.scatter(RA, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax2.set_title('b vs RA', fontsize=24,y=0.85)
        # add vertical line:
        #ax2.plot((np.min(RA),np.max(RA)), (0,0), 'k-', linewidth=1)
        ax2.axhline(y=0, linewidth=1, color="k", alpha=0.5)



        ax2.set_xlabel('RA', fontsize=20)
        #ax2.set_ylabel('a', fontsize=20)

        #ax2.set_position([0.5, 0.6, 0.4, 0.4])

        ax2.set_ylim([-3,4])
        ax2.set_yticks(np.arange(-3,4.1,1))


        #ax3

        ax3.scatter(DEC, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax3.set_title('b vs DEC', fontsize=24,y=0.85)

        # add vertical line:
        #ax3.plot((np.min(DEC),np.max(DEC)), (0,0), 'k-', linewidth=1)
        ax3.axhline(y=0, linewidth=1, color="k", alpha=0.5)


        ax3.set_xlabel('DEC', fontsize=20)
        #ax3.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax3.set_ylim([-3,4])
        ax3.set_yticks(np.arange(-3,4.1,1))

        #ax3.set_position([0, 0, 0.4, 0.4])



        #ax4

        ax4.scatter(Fiber, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax4.set_title('b vs FiberID', fontsize=24,y=0.85)

        # add vertical line:
        ax4.axhline(y=0, linewidth=1, color="k", alpha=0.5)
        #ax4.plot((np.min(Fiber),np.max(Fiber)), (0,0), 'k-', linewidth=1)



        ax4.set_xlabel('FIberID', fontsize=20)
        ax4.set_ylabel('b', fontsize=20)

        ax4.set_ylim([-3,4])
        ax4.set_yticks(np.arange(-3,4.1,1))

        #ax4.set_position([0.5,0, 0.4, 0.4])


        #ax5

        ax5.scatter(airmass, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax5.set_title('b vs air mass', fontsize=24,y=0.85)

        # add vertical line:
        ax5.axhline(y=0, linewidth=1, color="k", alpha=0.5)
        #ax5.plot((np.min(airmass),np.max(airmass)), (0,0), 'k-', linewidth=1)



        ax5.set_xlabel('FIberID', fontsize=20)
        #ax5.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax5.set_ylim([-3,4])
        ax5.set_yticks(np.arange(-3,4.1,1))


        #ax6

        ax6.scatter(SNR, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax6.set_title('b vs SNR', fontsize=24,y=0.85)

        # add vertical line:
        ax6.axhline(y=0, linewidth=1,color="k",alpha=0.5)
        #ax6.plot((np.min(SNR),np.max(SNR)), (0,0), )



        ax6.set_xlabel('SNR', fontsize=20)
        #ax6.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax6.set_ylim([-3,4])
        ax6.set_yticks(np.arange(-3,4.1,1))


        f.subplots_adjust(right=0.8)

        pl = ax1.scatter(HJD, b, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)


        cb.set_label("Mean inverse variance", fontsize=20)
        f.suptitle("Parameter b from the absorption line vs HJD, RA, DEC, FiberID, Airmass and SNR", fontsize=30)

        plt.show()


## c vs them

    # RV vs HJD RA DEC Fiber Airmass
    def c_subplot_5(self,c,HJD,Fiber,RA,DEC,airmass,mean_ivar,SNR):

        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)

        fig = plt.figure()


        f, ((ax1, ax2,ax3), (ax4, ax5,ax6)) = \
            plt.subplots(2, 3)


        alpha = 0.3
        #ax1

        ax1.scatter(HJD, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax1.set_title('c vs HJD', fontsize=24,y=0.85)

        ax1.set_xlabel('HJD', fontsize=20)
        ax1.set_ylabel('c', fontsize=20)

        # add vertical line:
        #ax1.plot((np.min(HJD),np.max(HJD)), (0,0), 'k-', linewidth=1)
        ax1.axhline(y=0, linewidth=1, color="k", alpha=0.5)


        ax1.set_ylim([-2,3])
        ax1.set_yticks(np.arange(-2,3.1,1))

        #ax1.set_position([0,0.6,0.4,0.4])


        #ax2

        ax2.scatter(RA, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax2.set_title('c vs RA', fontsize=24,y=0.85)
        # add vertical line:
        #ax2.plot((np.min(RA),np.max(RA)), (0,0), 'k-', linewidth=1)
        ax2.axhline(y=0, linewidth=1, color="k", alpha=0.5)



        ax2.set_xlabel('RA', fontsize=20)
        #ax2.set_ylabel('a', fontsize=20)

        #ax2.set_position([0.5, 0.6, 0.4, 0.4])

        ax2.set_ylim([-2,3])
        ax2.set_yticks(np.arange(-2,3.1,1))


        #ax3

        ax3.scatter(DEC, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax3.set_title('c vs DEC', fontsize=24,y=0.85)

        # add vertical line:
        #ax3.plot((np.min(DEC),np.max(DEC)), (0,0), 'k-', linewidth=1)
        ax3.axhline(y=0, linewidth=1, color="k", alpha=0.5)


        ax3.set_xlabel('DEC', fontsize=20)
        #ax3.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax3.set_ylim([-2,3])
        ax3.set_yticks(np.arange(-2,3.1,1))

        #ax3.set_position([0, 0, 0.4, 0.4])



        #ax4

        ax4.scatter(Fiber, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax4.set_title('c vs FiberID', fontsize=24,y=0.85)

        # add vertical line:
        ax4.axhline(y=0, linewidth=1, color="k", alpha=0.5)
        #ax4.plot((np.min(Fiber),np.max(Fiber)), (0,0), 'k-', linewidth=1)



        ax4.set_xlabel('FIberID', fontsize=20)
        ax4.set_ylabel('c', fontsize=20)
        ax4.set_ylim([-2,3])
        ax4.set_yticks(np.arange(-2,3.1,1))

        #ax4.set_position([0.5,0, 0.4, 0.4])


        #ax5

        ax5.scatter(airmass, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax5.set_title('c vs air mass', fontsize=24,y=0.85)

        # add vertical line:
        ax5.axhline(y=0, linewidth=1, color="k", alpha=0.5)
        #ax5.plot((np.min(airmass),np.max(airmass)), (0,0), 'k-', linewidth=1)



        ax5.set_xlabel('FIberID', fontsize=20)
        #ax5.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax5.set_ylim([-2,3])
        ax5.set_yticks(np.arange(-2,3.1,1))


        #ax6

        ax6.scatter(SNR, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax6.set_title('c vs SNR', fontsize=24,y=0.85)

        # add vertical line:
        ax6.axhline(y=0, linewidth=1,color="k",alpha=0.5)
        #ax6.plot((np.min(SNR),np.max(SNR)), (0,0), )



        ax6.set_xlabel('SNR', fontsize=20)
        #ax6.set_ylabel('RV shifts $m/s$', fontsize=20)
        ax6.set_ylim([-2,3])
        ax6.set_yticks(np.arange(-2,3.1,1))


        f.subplots_adjust(right=0.8)

        pl = ax1.scatter(HJD, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)


        cb.set_label("Mean inverse variance", fontsize=20)
        f.suptitle("Parameter c from the whole spectrum vs HJD, RA, DEC, FiberID, Airmass and SNR", fontsize=30)

        plt.show()



    # RV vs HJD RA DEC Fiber Airmass
    def c_new_subplot_5(self,c,HJD,Fiber,RA,DEC,airmass,mean_ivar,SNR):

        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)

        fig = plt.figure()


        f, ((ax1, ax2,ax3), (ax4, ax5,ax6)) = \
            plt.subplots(2, 3)


        alpha = 0.3
        #ax1

        ax1.scatter(HJD, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax1.set_title('c vs HJD', fontsize=24,y=0.85)

        ax1.set_xlabel('HJD', fontsize=20)
        ax1.set_ylabel('c', fontsize=20)

        # add vertical line:
        #ax1.plot((np.min(HJD),np.max(HJD)), (0,0), 'k-', linewidth=1)
        ax1.axhline(y=0, linewidth=1, color="k", alpha=0.5)


        ax1.set_ylim([-2,3])
        ax1.set_yticks(np.arange(-2,3.1,1))

        #ax1.set_position([0,0.6,0.4,0.4])


        #ax2

        ax2.scatter(RA, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax2.set_title('c vs RA', fontsize=24,y=0.85)
        # add vertical line:
        #ax2.plot((np.min(RA),np.max(RA)), (0,0), 'k-', linewidth=1)
        ax2.axhline(y=0, linewidth=1, color="k", alpha=0.5)



        ax2.set_xlabel('RA', fontsize=20)
        #ax2.set_ylabel('a', fontsize=20)

        #ax2.set_position([0.5, 0.6, 0.4, 0.4])

        ax2.set_ylim([-2,3])
        ax2.set_yticks(np.arange(-2,3.1,1))


        #ax3

        ax3.scatter(DEC, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax3.set_title('c vs DEC', fontsize=24,y=0.85)

        # add vertical line:
        #ax3.plot((np.min(DEC),np.max(DEC)), (0,0), 'k-', linewidth=1)
        ax3.axhline(y=0, linewidth=1, color="k", alpha=0.5)


        ax3.set_xlabel('DEC', fontsize=20)
        #ax3.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax3.set_ylim([-2,3])
        ax3.set_yticks(np.arange(-2,3.1,1))

        #ax3.set_position([0, 0, 0.4, 0.4])



        #ax4

        ax4.scatter(Fiber, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax4.set_title('c vs FiberID', fontsize=24,y=0.85)

        # add vertical line:
        ax4.axhline(y=0, linewidth=1, color="k", alpha=0.5)
        #ax4.plot((np.min(Fiber),np.max(Fiber)), (0,0), 'k-', linewidth=1)



        ax4.set_xlabel('FIberID', fontsize=20)
        ax4.set_ylabel('c', fontsize=20)
        ax4.set_ylim([-2,3])
        ax4.set_yticks(np.arange(-2,3.1,1))

        #ax4.set_position([0.5,0, 0.4, 0.4])


        #ax5

        ax5.scatter(airmass, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax5.set_title('c vs air mass', fontsize=24,y=0.85)

        # add vertical line:
        ax5.axhline(y=0, linewidth=1, color="k", alpha=0.5)
        #ax5.plot((np.min(airmass),np.max(airmass)), (0,0), 'k-', linewidth=1)



        ax5.set_xlabel('FIberID', fontsize=20)
        #ax5.set_ylabel('RV shifts $m/s$', fontsize=20)

        ax5.set_ylim([-2,3])
        ax5.set_yticks(np.arange(-2,3.1,1))


        #ax6

        ax6.scatter(SNR, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        ax6.set_title('c vs SNR', fontsize=24,y=0.85)

        # add vertical line:
        ax6.axhline(y=0, linewidth=1,color="k",alpha=0.5)
        #ax6.plot((np.min(SNR),np.max(SNR)), (0,0), )



        ax6.set_xlabel('SNR', fontsize=20)
        #ax6.set_ylabel('RV shifts $m/s$', fontsize=20)
        ax6.set_ylim([-2,3])
        ax6.set_yticks(np.arange(-2,3.1,1))


        f.subplots_adjust(right=0.8)

        pl = ax1.scatter(HJD, c, marker='x', c=mean_ivar,
                    vmin=10000, vmax=40000, alpha=alpha)

        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)


        cb.set_label("Mean inverse variance", fontsize=20)
        f.suptitle("Parameter c from the absorption line vs HJD, RA, DEC, FiberID, Airmass and SNR", fontsize=30)

        plt.show()

    def old_vs_new_sim(self):

        parameters = self.parameters

        parameters_sim = self.parameters_sim

        a = parameters[:,0]
        b = parameters[:,1]
        c = parameters[:,2]

        RV = (c-a)/(a+b+c)*4144.68

        a = parameters_sim[:, 0]
        b = parameters_sim[:, 1]
        c = parameters_sim[:, 2]

        RV_sim = (c - a) / (a + b + c) * 4144.68

        font = {'weight': 'bold', 'size': 12}
        matplotlib.rc('font', **font)

        f, ((ax1, ax2), (ax3, ax4)) = \
            plt.subplots(2, 2)

        """

        a = parameters_new[:,0]
        b = parameters_new[:,1]
        c = parameters_new[:,2]
        RV_n = (c-a)/(a+b+c)*4144.68

        """

        alpha = 0.3

#rv
        ax1.plot(RV,RV_sim,"ro",label="RV shifts", markersize=3,alpha = alpha)
        ax1.plot(RV,RV,"k-")
        ax1.set_xlabel("RV shifts by fitting the whole spectrum $m/s$", fontsize=14)
        ax1.set_ylabel("RV shifts by fitting simultaneously $m/s$", fontsize=14)

#a

        ax2.plot(parameters[:,0],parameters_sim[:,0],"ro",label="Parameter a",markersize=3,alpha = alpha)
        ax2.plot(parameters[:,0],parameters[:,0], "k-")
        ax2.set_xlabel("a by fitting the whole spectrum", fontsize=14)
        ax2.set_ylabel("a by fitting simultaneously", fontsize=14)



#b

        ax3.plot(parameters[:,1],parameters_sim[:,1],"ro",label="Parameter b", markersize=3,alpha = alpha)
        ax3.plot(parameters[:,1],parameters[:,1], "k-")
        ax3.set_xlabel("b by fitting the whole spectrum", fontsize=14)
        ax3.set_ylabel("b by fitting simultaneously", fontsize=14)


#c

        ax4.plot(parameters[:,2],parameters_sim[:,2],"ro",label="Parameter c", markersize=3,alpha = alpha)
        ax4.plot(parameters[:, 2], parameters[:, 2], "k-")
        ax4.set_xlabel("c by fitting the whole spectrum", fontsize=14)
        ax4.set_ylabel("c by fitting simultaneously", fontsize=14)


        f.suptitle("Comparison of the old and the new method", fontsize=20)


        plt.show()

    def choose_four_biggest_RV_for_new_method(self):
        # return index
        N = len(self.velocity_new)
        index = self.velocity_new.argsort()[N-8:N-4]
        print(self.velocity_new[index])

        # from small to big
        print(index)

        return index

    def choose_four_biggest_delta_rv(self):
        N = len(self.velocity_new)
        index = abs(self.velocity-self.velocity_new).argsort()[N-8:N-4]

        print(index)
        return index



    def plot_visit_sim_old(self,index):


        # only choose individual visits:
        N = len(index)
        for i in range(0,N):

            # mask
            plt.subplot(N,2,2*i+1)
            star = fits.open(self.star_name[index[i]])

            name = str(self.star_name[index[i]]).replace(".fits", "")
            name = name.replace("/Users/caojunzhi/Desktop/Data/n_900/", "")

            ind = self.star_visit[index[i]]
            flux = star[0].data[ind+2,:]
            inf = star[2].data[ind + 2, :]
            mask = star[13].data[ind+2,:]

            plt.step(wl, flux, "k", label="One visit of star %s"%name, linewidth=0.7, alpha=1)
            plt.plot(wl, flux*mask, "ro", label="The absorption line", markersize=1, alpha=0.5)
            plt.ylabel("Flux", fontsize=20)
            axes = plt.gca()
            axes.set_xlim([15660, 15780])
            # axes.set_xlim([16160,16280])
            axes.set_ylim([0.5, 1.5])
            plt.legend()

            # inf

            plt.subplot(N,2,2*i+2)

            plt.step(wl, flux, "k", label="From the whole spectrum RV=%.2f $m/s$ a=%.2f b=%.2f c=%.2f" % (
            self.velocity[index[i]], self.parameters[index[i],0], self.parameters[index[i],1], self.parameters[index[i],2]), linewidth=0.7, alpha=1)
            plt.plot(wl, inf, "b", label="From the absorption line RV=%.2f $m/s$ a=%.2f b=%.2f c=%.2f" % (
            self.velocity_new[index[i]], self.parameters_new[index[i],0], self.parameters_new[index[i],1], self.parameters_new[index[i],2]), linewidth=0.7,
                     alpha=0.5)

            # plt.errorbar(wl,flux[i], ecolor='k', alpha=0.02, capthick=0.2, yerr=ivar[i]**(-0.5))


            axes = plt.gca()
            axes.set_xlim([15660, 15780])
            # axes.set_xlim([16160,16280])
            axes.set_ylim([0.5, 1.5])
            # axes.set_yticks(np.arange(0.8,1.21,0.1))

            # plt.xlabel("Wave length $\AA$", fontsize=20)
            #plt.ylabel("Flux", fontsize=20)

            plt.legend()
        plt.suptitle("The fitting result of visits with the biggest delta RV shifts", fontsize=20)

        plt.show()




"""
p = 103
star_choose = fits.open(path_fits[p])

nor_i = star_choose[0].data
ivar_i = star_choose[1].data
inf_i = star_choose[2].data

print(nor_i.shape,ivar_i.shape)


"""
model = plot()
# plot continuum pixel
#model.plot_continuum_pixel_single_star(flux=nor_i,ivar=ivar_i)

# read data
model.read_data()

#model.histogram_2_2_rv_abc()

# plot 2*3 histogram of RV

#model.ve_sim_subplot_5(HJD=model.HJD,Fiber=model.fiber_id,RA=model.RA,DEC=model.DEC,airmass=model.airmass,mean_ivar=model.meanivar,SNR=model.SNR)

# a
#model.a_subplot_5(a=model.parameters[:,0],HJD=model.HJD,Fiber=model.fiber_id,RA=model.RA,DEC=model.DEC,airmass=model.airmass,mean_ivar=model.meanivar,SNR=model.SNR)
#model.a_new_subplot_5(a=model.parameters_new[:,0],HJD=model.HJD,Fiber=model.fiber_id,RA=model.RA,DEC=model.DEC,airmass=model.airmass,mean_ivar=model.meanivar,SNR=model.SNR)

# b vs all

#model.b_subplot_5(b=model.parameters[:,1],HJD=model.HJD,Fiber=model.fiber_id,RA=model.RA,DEC=model.DEC,airmass=model.airmass,mean_ivar=model.meanivar,SNR=model.SNR)
#model.b_new_subplot_5(b=model.parameters_new[:,1],HJD=model.HJD,Fiber=model.fiber_id,RA=model.RA,DEC=model.DEC,airmass=model.airmass,mean_ivar=model.meanivar,SNR=model.SNR)



# c vs all

#model.c_subplot_5(c=model.parameters[:,2],HJD=model.HJD,Fiber=model.fiber_id,RA=model.RA,DEC=model.DEC,airmass=model.airmass,mean_ivar=model.meanivar,SNR=model.SNR)
#model.c_new_subplot_5(c=model.parameters_new[:,2],HJD=model.HJD,Fiber=model.fiber_id,RA=model.RA,DEC=model.DEC,airmass=model.airmass,mean_ivar=model.meanivar,SNR=model.SNR)

# old vs new


model.old_vs_new_sim()
# choose five stars and plot:
# individual visit:
# save in peak 1,2,3,4,5


#name = [628,577,461,321,77]
#path = path_fits[name[3]]

#model.plot_single_star_mask_result(path=path)





#big_4 = model.choose_four_biggest_RV_for_new_method()
#big_delta_4 = model.choose_four_biggest_delta_rv()

#model.plot_visit_mask_result(big_4)