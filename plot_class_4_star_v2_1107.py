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
import math
import re



def log10(x):
    return math.log10(x)

# sinc interpolation

# Distance between s is equal and they are your raw data
# u is what you want. Then length of u is not equal to s.



class plot_4:
    def __init__(self):
        n = 1

    def plot_name_4(self, path_array):

        # import data

        pkl_file = open('wl.pkl', 'rb')
        wl = pickle.load(pkl_file)
        pkl_file.close()

        all_set_flux = []
        all_set_ivar = []
        star_name = []
        inf_old = []
        inf_opt = []
        inf_label = []
        parameters_4 = []
        delta_chi = []

        n_star = len(path_array)
        for i in range(0,n_star):

            star = fits.open(path_array[i])
            all_set_flux.append(star[0].data[0])
            all_set_ivar.append(star[1].data[0])
            star_name_i = path_array[i]
            star_name_i = re.sub('/Users/caojunzhi/Desktop/Data/n_900/', '', star_name_i)

            star_name_i = re.sub('\.fits', '',star_name_i)
            star_name.append(star_name_i)

            inf_old.append(star[2].data[0])
            inf_opt.append(star[3].data[0])
            inf_label.append(star[9].data[0])
            parameters_4.append(star[4].data[0])
            delta_chi.append(star[6].data[0])

        all_set_flux = np.array(all_set_flux)
        all_set_ivar = np.array(all_set_ivar)
        star_name = np.array(star_name)
        inf_old = np.array(inf_old)
        inf_opt = np.array(inf_opt)
        inf_label = np.array(inf_label)
        parameters_4 = np.array(parameters_4)
        delta_chi = np.array(delta_chi)

        error = all_set_ivar ** (-0.5)

        # plot
        nor = all_set_flux


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

        ax3.step(wl, nor[p2, :], "k", label=l2, alpha=1, linewidth=1.5)
        ax3.errorbar(wl, nor[p2, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax3.plot(wl, inf_opt[p2, :], "r", alpha=trans, linewidth=1.5)
        ax3.plot(wl, inf_old[p2, :], "g", alpha=trans, linewidth=1.5)

        ax3.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax3.set_xlim([15660, 15780])
        ax3.set_ylim([0.8, 1.2])
        ax3.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax4

        ax4.step(wl, nor[p2, :], "k", label=delta_para_2, alpha=1, linewidth=1.5)
        ax4.errorbar(wl, nor[p2, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax4.plot(wl, inf_opt[p2, :], "r", alpha=trans, linewidth=1.5)
        ax4.plot(wl, inf_old[p2, :], "g", alpha=trans, linewidth=1.5)

        ax4.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax4.set_xlim([16160, 16280])
        ax4.set_ylim([0.8, 1.2])
        ax4.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax5

        ax5.step(wl, nor[p3, :], "k", label=l3, alpha=1, linewidth=1.5)
        ax5.errorbar(wl, nor[p3, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax5.plot(wl, inf_opt[p3, :], "r", alpha=trans, linewidth=1.5)
        ax5.plot(wl, inf_old[p3, :], "g", alpha=trans, linewidth=1.5)

        ax5.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax5.set_xlim([15660, 15780])
        ax5.set_ylim([0.8, 1.2])
        ax5.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax6

        ax6.step(wl, nor[p3, :], "k", label=delta_para_3, alpha=1, linewidth=1.5)
        ax6.errorbar(wl, nor[p3, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax6.plot(wl, inf_opt[p3, :], "r", alpha=trans, linewidth=1.5)
        ax6.plot(wl, inf_old[p3, :], "g", alpha=trans, linewidth=1.5)

        ax6.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax6.set_xlim([16160, 16280])
        ax6.set_ylim([0.8, 1.2])
        ax6.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax7

        ax7.step(wl, nor[p4, :], "k", label=l4, alpha=1, linewidth=1.5)
        ax7.errorbar(wl, nor[p4, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax7.plot(wl, inf_opt[p4, :], "r", alpha=trans, linewidth=1.5)
        ax7.plot(wl, inf_old[p4, :], "g", alpha=trans, linewidth=1.5)

        ax7.legend(bbox_to_anchor=(0, 0.65), loc=3,
                   ncol=1)
        ax7.set_xlim([15660, 15780])
        ax7.set_ylim([0.8, 1.2])
        ax7.set_yticks(np.arange(0.8, 1.2, 0.1))

        fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

        # ax8

        ax8.step(wl, nor[p4, :], "k", label=delta_para_4, alpha=1, linewidth=1.5)
        ax8.errorbar(wl, nor[p4, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
        ax8.plot(wl, inf_opt[p4, :], "r", alpha=trans, linewidth=1.5)
        ax8.plot(wl, inf_old[p4, :], "g", alpha=trans, linewidth=1.5)

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






