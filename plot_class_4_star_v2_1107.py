import matplotlib.pyplot as plt
import random
import matplotlib
import numpy as np
import pickle
from astropy.table import Table
import TheCannon_2
from TheCannon_2 import dataset,apogee
from astropy.table import Table
import numpy as np
import os
from astropy.io import fits
import pickle
import AnniesLasso_2 as tc


class plot_4:
    def __init__(self):
        n = 1

    def train(self):

        training_set = Table.read("reference_labels.csv")

        pkl_file = open('training_set_flux.pkl', 'rb')
        training_set_flux = pickle.load(pkl_file)
        pkl_file.close()

        # read table
        pkl_file = open('training_set_ivar.pkl', 'rb')
        training_set_ivar = pickle.load(pkl_file)
        pkl_file.close()



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
        self.model =model

    def plot_v2(self,flux,ivar,error,star_name):

        #import data

        pkl_file = open('wl.pkl', 'rb')
        wl = pickle.load(pkl_file)
        pkl_file.close()

        all_set_flux = np.array(flux)
        all_set_ivar = np.array(ivar)
        all_set_error = np.array(error)
        star_name = np.array(star_name)

        inf_label = self.model.fit(all_set_flux, all_set_ivar)
        v = self.model.vectorizer.get_label_vector(inf_label)
        # v_500 = model.vectorizer.get_label_vector(test_label)

        inf_flux = np.dot(v, self.model.theta.T)
        # opt four star





        # plot
        nor = all_set_flux
        inf_old = inf_flux

        inf_opt, parameters_4 = self.model.fitting_spectrum_parameters_single(
            all_set_flux, all_set_ivar, inf_old)

        error = all_set_error

        delta_chi = self.model.delta_chi_squared(all_set_flux, all_set_ivar, inf_flux)

        # plot-opt
        trans = 0.5
        p1 = 0
        p2 = 1
        p3 = 2
        p4 = 3

        n1 = star_name[0]
        n2 = star_name[1]
        n3 = star_name[2]
        n4 = star_name[3]
        print(n1,n2,n3,n4)

        # test label

        l1 = "Teff=%.3f Log g = %.3f" % (inf_label[0,0],inf_label[0,1])
        l2 = "Teff=%.3f Log g = %.3f" % (inf_label[1, 0], inf_label[1, 1])
        l3 = "Teff=%.3f Log g = %.3f" % (inf_label[2, 0], inf_label[2, 1])
        l4 = "Teff=%.3f Log g = %.3f" % (inf_label[3, 0], inf_label[3, 1])

        # delta-chi and parameters a, b and c.
        delta_para_1 = "delta-chi-squared=%f a=%.3f b=%.3f c=%.3f" % (
        delta_chi[0], parameters_4[0, 0], parameters_4[0, 1], parameters_4[0, 2])

        delta_para_2 = "delta-chi-squared=%f a=%.3f b=%.3f c=%.3f" % (
        delta_chi[1], parameters_4[1, 0], parameters_4[1, 1], parameters_4[1, 2])

        delta_para_3 = "delta-chi-squared=%f a=%.3f b=%.3f c=%.3f" % (
        delta_chi[2], parameters_4[2, 0], parameters_4[2, 1], parameters_4[2, 2])

        delta_para_4 = "delta-chi-squared=%f a=%.3f b=%.3f c=%.3f" % (
        delta_chi[3], parameters_4[3, 0], parameters_4[3, 1], parameters_4[3, 2])

        print(parameters_4)

        ## Let's plot

        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)
        fig = plt.figure()

        f, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = \
            plt.subplots(4, 2, sharex='col', sharey='row')

        # ax1

        ax1.step(wl, nor[p1, :], "k",label =l1, alpha=1, linewidth=1.5)
        ax1.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax1.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax1.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax1.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax1.set_xlim([15660, 15780])
        ax1.set_ylim([0.8, 1.21])
        ax1.set_yticks(np.arange(0.8, 1.21, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax2

        ax2.step(wl, nor[p1, :], "k",label = delta_para_1, alpha=1, linewidth=1.5)
        ax2.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax2.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax2.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax2.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax2.set_xlim([16160, 16280])
        ax2.set_ylim([0.8, 1.21])
        ax2.set_yticks(np.arange(0.8, 1.21, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax3

        ax3.step(wl, nor[p1, :], "k", label =l2,alpha=1, linewidth=1.5)
        ax3.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax3.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax3.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax3.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax3.set_xlim([15660, 15780])
        ax3.set_ylim([0.8, 1.2])
        ax3.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax4

        ax4.step(wl, nor[p1, :], "k",label =delta_para_2, alpha=1, linewidth=1.5)
        ax4.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax4.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax4.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax4.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax4.set_xlim([16160, 16280])
        ax4.set_ylim([0.8, 1.2])
        ax4.set_yticks(np.arange(0.8, 1.2, 0.1))


        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax5

        ax5.step(wl, nor[p1, :], "k",label = l3, alpha=1, linewidth=1.5)
        ax5.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax5.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax5.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax5.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax5.set_xlim([15660, 15780])
        ax5.set_ylim([0.8, 1.2])
        ax5.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax6

        ax6.step(wl, nor[p1, :], "k",label =delta_para_3, alpha=1, linewidth=1.5)
        ax6.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax6.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax6.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax6.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax6.set_xlim([16160, 16280])
        ax6.set_ylim([0.8, 1.2])
        ax6.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax7

        ax7.step(wl, nor[p1, :], "k", label =l4,alpha=1, linewidth=1.5)
        ax7.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax7.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax7.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax7.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax7.set_xlim([15660, 15780])
        ax7.set_ylim([0.8, 1.2])
        ax7.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax8

        ax8.step(wl, nor[p1, :], "k",label=delta_para_4, alpha=1, linewidth=1.5)
        ax8.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax8.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax8.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax8.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax8.set_xlim([16160, 16280])
        ax8.set_ylim([0.8, 1.2])
        ax8.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        ax3.set_ylabel('Normalized flux', fontsize=30)
        ax3.yaxis.set_label_coords(-0.05, 0)

        ax7.set_xlabel('Wave length $\AA$ ', fontsize=25)
        ax8.set_xlabel('Wave length $\AA$ ', fontsize=25)

        f.subplots_adjust(hspace=0)

        # Don't use this
        # plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        
        plt.show()
        #f.savefig(delta_para_1+".png",dpi=400,papertype = "a4")

    def plot(self, flux, ivar, error, star_name, test_labels,parameters_4):
        # import data

        pkl_file = open('wl.pkl', 'rb')
        wl = pickle.load(pkl_file)
        pkl_file.close()

        all_set_flux = np.array(flux)
        all_set_ivar = np.array(ivar)
        all_set_error = np.array(error)
        star_name = np.array(star_name)

        #inf_label = self.model.fit(all_set_flux, all_set_ivar)

###################################
        # bugs
        inf_label = test_labels


        v = self.model.vectorizer.get_label_vector(inf_label)
        # v_500 = model.vectorizer.get_label_vector(test_label)

        inf_flux = np.dot(v, self.model.theta.T)
        # opt four star





        # plot
        nor = all_set_flux
        inf_old = inf_flux

        inf_opt, p = self.model.fitting_spectrum_parameters_single(
            all_set_flux, all_set_ivar, inf_old)

        print(parameters_4)

        #error = all_set_error

        delta_chi = self.model.delta_chi_squared(all_set_flux, all_set_ivar, inf_flux)

        # plot-opt
        trans = 0.5
        p1 = 0
        p2 = 1
        p3 = 2
        p4 = 3

        n1 = star_name[0]
        n2 = star_name[1]
        n3 = star_name[2]
        n4 = star_name[3]
        print(n1, n2, n3, n4)

        # test label

        l1 = "Teff=%.3f Log g = %.3f" % (inf_label[0, 0], inf_label[0, 1])
        l2 = "Teff=%.3f Log g = %.3f" % (inf_label[1, 0], inf_label[1, 1])
        l3 = "Teff=%.3f Log g = %.3f" % (inf_label[2, 0], inf_label[2, 1])
        l4 = "Teff=%.3f Log g = %.3f" % (inf_label[3, 0], inf_label[3, 1])

        #bug
        # para_4 is obtained by the optimization process using APOGEE team labels.


        # delta-chi and parameters a, b and c.
        delta_para_1 = "delta-chi-squared=%f a=%.3f b=%.3f c=%.3f" % (
            delta_chi[0], parameters_4[0, 0], parameters_4[0, 1], parameters_4[0, 2])

        delta_para_2 = "delta-chi-squared=%f a=%.3f b=%.3f c=%.3f" % (
            delta_chi[1], parameters_4[1, 0], parameters_4[1, 1], parameters_4[1, 2])

        delta_para_3 = "delta-chi-squared=%f a=%.3f b=%.3f c=%.3f" % (
            delta_chi[2], parameters_4[2, 0], parameters_4[2, 1], parameters_4[2, 2])

        delta_para_4 = "delta-chi-squared=%f a=%.3f b=%.3f c=%.3f" % (
            delta_chi[3], parameters_4[3, 0], parameters_4[3, 1], parameters_4[3, 2])

        print(parameters_4)

        # star-1

        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)
        fig = plt.figure()

        plt.subplot(4, 2, 1)

        plt.step(wl, nor[p1, :], "k", label=delta_para_1, alpha=1, linewidth=1.5)
        # plt.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 4, capthick=0.2, yerr=error[p1, :])
        plt.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        plt.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        plt.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        fig.suptitle('Comparison of the spectrum', fontsize=24)

        plt.xlabel('wave length', fontsize=20)
        plt.ylabel('Spectrum', fontsize=20)

        axes = plt.gca()
        axes.set_xlim([15660, 15780])

        # axes.set_xlim([16160,16280])
        axes.set_ylim([0.8, 1.21])
        axes.set_yticks(np.arange(0.8, 1.21, 0.1))

        plt.subplot(4, 2, 2)
        plt.step(wl, nor[p1, :], "k", label=l1, alpha=1, linewidth=1.5)
        # plt.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 4, capthick=0.2, yerr=error[p1, :])
        plt.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        plt.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        plt.legend(bbox_to_anchor=(0.1, 0.7), loc=3,
                   ncol=2, mode="expand", borderaxespad=0.)
        fig.suptitle('Comparison of the spectrum', fontsize=24)
        plt.xlabel('wave length', fontsize=20)
        plt.ylabel('Spectrum', fontsize=20)

        axes = plt.gca()
        # axes.set_xlim([15660,15780])
        axes.set_xlim([16160, 16280])
        axes.set_ylim([0.8, 1.21])
        axes.set_yticks(np.arange(0.8, 1.21, 0.1))

        # star-2

        plt.subplot(4, 2, 3)

        plt.step(wl, nor[p2, :], "k", label=delta_para_2, alpha=1, linewidth=1.5)
        # plt.errorbar(wl, nor[p2, :], ecolor='k', alpha=trans / 4, capthick=0.2, yerr=error[p1, :])
        plt.plot(wl, inf_opt[p2, :], "r", alpha=trans, linewidth=1.5)
        plt.plot(wl, inf_old[p2, :], "g", alpha=trans, linewidth=1.5)

        plt.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        # fig.suptitle('Comparison of the spectrum', fontsize=24)
        plt.xlabel('wave length', fontsize=20)
        plt.ylabel('Spectrum', fontsize=20)

        axes = plt.gca()
        axes.set_xlim([15660, 15780])

        # axes.set_xlim([16160,16280])
        axes.set_ylim([0.8, 1.21])
        axes.set_yticks(np.arange(0.8, 1.21, 0.1))

        plt.subplot(4, 2, 4)
        plt.step(wl, nor[p2, :], "k", label=l2, alpha=1, linewidth=1.5)
        # plt.errorbar(wl, nor[p2, :], ecolor='k', alpha=trans / 4, capthick=0.2, yerr=error[p1, :])
        plt.plot(wl, inf_opt[p2, :], "r", alpha=trans, linewidth=1.5)
        plt.plot(wl, inf_old[p2, :], "g", alpha=trans, linewidth=1.5)

        plt.legend(bbox_to_anchor=(0.1, 0.7), loc=3,
                   ncol=2, mode="expand", borderaxespad=0.)
        # fig.suptitle('Comparison of the spectrum', fontsize=24)
        plt.xlabel('wave length', fontsize=20)
        plt.ylabel('Spectrum', fontsize=20)

        axes = plt.gca()
        # axes.set_xlim([15660,15780])
        axes.set_xlim([16160, 16280])
        axes.set_ylim([0.8, 1.21])
        axes.set_yticks(np.arange(0.8, 1.21, 0.1))
        # star-3

        plt.subplot(4, 2, 5)

        plt.step(wl, nor[p3, :], "k", label=delta_para_3, alpha=1, linewidth=1.5)
        # plt.errorbar(wl, nor[p3, :], ecolor='k', alpha=trans / 4, capthick=0.2, yerr=error[p1, :])
        plt.plot(wl, inf_opt[p3, :], "r", alpha=trans, linewidth=1.5)
        plt.plot(wl, inf_old[p3, :], "g", alpha=trans, linewidth=1.5)

        plt.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        # fig.suptitle('Comparison of the spectrum', fontsize=24)
        plt.xlabel('wave length', fontsize=20)
        plt.ylabel('Spectrum', fontsize=20)

        axes = plt.gca()
        axes.set_xlim([15660, 15780])

        # axes.set_xlim([16160,16280])
        axes.set_ylim([0.8, 1.21])
        axes.set_yticks(np.arange(0.8, 1.21, 0.1))

        plt.subplot(4, 2, 6)
        plt.step(wl, nor[p3, :], "k", label=l3, alpha=1, linewidth=1.5)
        # plt.errorbar(wl, nor[p3, :], ecolor='k', alpha=trans / 4, capthick=0.2, yerr=error[p1, :])
        plt.plot(wl, inf_opt[p3, :], "r", alpha=trans, linewidth=1.5)
        plt.plot(wl, inf_old[p3, :], "g", alpha=trans, linewidth=1.5)

        plt.legend(bbox_to_anchor=(0.1, 0.7), loc=3,
                   ncol=2, mode="expand", borderaxespad=0.)
        # fig.suptitle('Comparison of the spectrum', fontsize=24)
        plt.xlabel('wave length', fontsize=20)
        plt.ylabel('Spectrum', fontsize=20)

        axes = plt.gca()
        # axes.set_xlim([15660,15780])
        axes.set_xlim([16160, 16280])
        axes.set_ylim([0.8, 1.21])
        axes.set_yticks(np.arange(0.8, 1.21, 0.1))
        # star-4

        plt.subplot(4, 2, 7)

        plt.step(wl, nor[p4, :], "k", label=delta_para_4, alpha=1, linewidth=1.5)
        # plt.errorbar(wl, nor[p4, :], ecolor='k', alpha=trans / 4, capthick=0.2, yerr=error[p1, :])
        plt.plot(wl, inf_opt[p4, :], "r", alpha=trans, linewidth=1.5)
        plt.plot(wl, inf_old[p4, :], "g", alpha=trans, linewidth=1.5)

        plt.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        # fig.suptitle('Comparison of the spectrum', fontsize=24)
        plt.xlabel('wave length', fontsize=20)
        plt.ylabel('Spectrum', fontsize=20)

        axes = plt.gca()
        axes.set_xlim([15660, 15780])

        # axes.set_xlim([16160,16280])
        axes.set_ylim([0.8, 1.21])
        axes.set_yticks(np.arange(0.8, 1.21, 0.1))

        plt.subplot(4, 2, 8)
        plt.step(wl, nor[p4, :], "k", label=l4, alpha=1, linewidth=1.5)
        # plt.errorbar(wl, nor[p4, :], ecolor='k', alpha=trans / 4, capthick=0.2, yerr=error[p1, :])
        plt.plot(wl, inf_opt[p4, :], "r", alpha=trans, linewidth=1.5)
        plt.plot(wl, inf_old[p4, :], "g", alpha=trans, linewidth=1.5)

        plt.legend(bbox_to_anchor=(0.1, 0.7), loc=3,
                   ncol=2, mode="expand", borderaxespad=0.)
        # fig.suptitle('Comparison of the spectrum', fontsize=24)
        plt.xlabel('wave length', fontsize=20)
        plt.ylabel('Spectrum', fontsize=20)

        axes = plt.gca()
        # axes.set_xlim([15660,15780])
        axes.set_xlim([16160, 16280])
        axes.set_ylim([0.8, 1.21])
        axes.set_yticks(np.arange(0.8, 1.21, 0.1))

        plt.show()

        # calculate chi-squared
        ivar = all_set_ivar

        for p in range(0, 4):
            ivar_r = ivar[p, :]
            ni = len(ivar_r)

            c = np.zeros((ni, ni))

            for i in range(0, ni):
                c[i, i] = ivar_r[i]

            # correct chi-squared
            a_old = np.dot(np.dot((nor[p, :] - inf_old[p, :]).T, c), (nor[p, :] - inf_old[p, :]))
            a_opt = np.dot(np.dot((nor[p, :] - inf_opt[p, :]).T, c), (nor[p, :] - inf_opt[p, :]))

            print(a_old, a_opt)

    def plot_v3(self, flux, ivar, error, star_name,inf_flux,inf_label):

        # import data

        pkl_file = open('wl.pkl', 'rb')
        wl = pickle.load(pkl_file)
        pkl_file.close()

        all_set_flux = np.array(flux)
        all_set_ivar = np.array(ivar)
        all_set_error = np.array(error)
        star_name = np.array(star_name)

        # plot
        nor = all_set_flux
        inf_old = inf_flux

        inf_opt, parameters_4 = self.model.fitting_spectrum_parameters_single(
            all_set_flux, all_set_ivar, inf_old)

        error = all_set_error

        delta_chi = self.model.delta_chi_squared(all_set_flux, all_set_ivar, inf_flux)

        # plot-opt
        trans = 0.5
        p1 = 0
        p2 = 1
        p3 = 2
        p4 = 3

        n1 = star_name[0]
        n2 = star_name[1]
        n3 = star_name[2]
        n4 = star_name[3]
        print(n1, n2, n3, n4)

        # test label

        l1 = "Teff=%.3f Log g = %.3f" % (inf_label[0, 0], inf_label[0, 1])
        l2 = "Teff=%.3f Log g = %.3f" % (inf_label[1, 0], inf_label[1, 1])
        l3 = "Teff=%.3f Log g = %.3f" % (inf_label[2, 0], inf_label[2, 1])
        l4 = "Teff=%.3f Log g = %.3f" % (inf_label[3, 0], inf_label[3, 1])

        # delta-chi and parameters a, b and c.
        delta_para_1 = "delta-chi-squared=%f a=%.3f b=%.3f c=%.3f" % (
            delta_chi[0], parameters_4[0, 0], parameters_4[0, 1], parameters_4[0, 2])

        delta_para_2 = "delta-chi-squared=%f a=%.3f b=%.3f c=%.3f" % (
            delta_chi[1], parameters_4[1, 0], parameters_4[1, 1], parameters_4[1, 2])

        delta_para_3 = "delta-chi-squared=%f a=%.3f b=%.3f c=%.3f" % (
            delta_chi[2], parameters_4[2, 0], parameters_4[2, 1], parameters_4[2, 2])

        delta_para_4 = "delta-chi-squared=%f a=%.3f b=%.3f c=%.3f" % (
            delta_chi[3], parameters_4[3, 0], parameters_4[3, 1], parameters_4[3, 2])

        print(parameters_4)

        ## Let's plot

        font = {'weight': 'bold', 'size': 13}
        matplotlib.rc('font', **font)
        fig = plt.figure()

        f, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = \
            plt.subplots(4, 2, sharex='col', sharey='row')

        # ax1

        ax1.step(wl, nor[p1, :], "k", label=l1, alpha=1, linewidth=1.5)
        ax1.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax1.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax1.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax1.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax1.set_xlim([15660, 15780])
        ax1.set_ylim([0.8, 1.21])
        ax1.set_yticks(np.arange(0.8, 1.21, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax2

        ax2.step(wl, nor[p1, :], "k", label=delta_para_1, alpha=1, linewidth=1.5)
        ax2.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax2.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax2.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax2.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax2.set_xlim([16160, 16280])
        ax2.set_ylim([0.8, 1.21])
        ax2.set_yticks(np.arange(0.8, 1.21, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax3

        ax3.step(wl, nor[p1, :], "k", label=l2, alpha=1, linewidth=1.5)
        ax3.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax3.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax3.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax3.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax3.set_xlim([15660, 15780])
        ax3.set_ylim([0.8, 1.2])
        ax3.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax4

        ax4.step(wl, nor[p1, :], "k", label=delta_para_2, alpha=1, linewidth=1.5)
        ax4.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax4.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax4.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax4.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax4.set_xlim([16160, 16280])
        ax4.set_ylim([0.8, 1.2])
        ax4.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax5

        ax5.step(wl, nor[p1, :], "k", label=l3, alpha=1, linewidth=1.5)
        ax5.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax5.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax5.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax5.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax5.set_xlim([15660, 15780])
        ax5.set_ylim([0.8, 1.2])
        ax5.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax6

        ax6.step(wl, nor[p1, :], "k", label=delta_para_3, alpha=1, linewidth=1.5)
        ax6.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax6.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax6.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax6.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax6.set_xlim([16160, 16280])
        ax6.set_ylim([0.8, 1.2])
        ax6.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax7

        ax7.step(wl, nor[p1, :], "k", label=l4, alpha=1, linewidth=1.5)
        ax7.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax7.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax7.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax7.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax7.set_xlim([15660, 15780])
        ax7.set_ylim([0.8, 1.2])
        ax7.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax8

        ax8.step(wl, nor[p1, :], "k", label=delta_para_4, alpha=1, linewidth=1.5)
        ax8.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax8.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
        ax8.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

        ax8.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax8.set_xlim([16160, 16280])
        ax8.set_ylim([0.8, 1.2])
        ax8.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        ax3.set_ylabel('Normalized flux', fontsize=30)
        ax3.yaxis.set_label_coords(-0.05, 0)

        ax7.set_xlabel('Wave length $\AA$ ', fontsize=25)
        ax8.set_xlabel('Wave length $\AA$ ', fontsize=25)

        f.subplots_adjust(hspace=0)

        # Don't use this
        # plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

        plt.show()
        # f.savefig(delta_para_1+".png",dpi=400,papertype = "a4")









