import os
import numpy as np
from numpy import genfromtxt
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle
import matplotlib.pyplot as plt
import matplotlib
import matplotlib as mpl

import AnniesLasso_2 as tc
import time


# load data

pkl_file = open('parameters_tr_inf.pkl', 'rb')
parameters_tr_inf = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('parameters_tr.pkl', 'rb')
parameters_tr = pickle.load(pkl_file)
pkl_file.close()


# labels
tr_label = genfromtxt("reference-3.csv",delimiter=",")
tr_label = np.array(tr_label)


pkl_file = open('inf_labels_tr_opt.pkl', 'rb')
inf_labels_tr_opt = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('inf_labels_tr.pkl', 'rb')
inf_labels_tr = pickle.load(pkl_file)
pkl_file.close()

#load n_900 data

pkl_file = open('n_testing_labels_900.pkl', 'rb')
ref_labels_900 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n_inf_labels_900_opt.pkl', 'rb')
inf_labels_900_opt = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n_inf_label_900.pkl', 'rb')
inf_labels_900 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n_testing_ivar_900.pkl', 'rb')
testing_ivar_900 = pickle.load(pkl_file)
pkl_file.close()

ivar= np.mean(testing_ivar_900,axis=1)
ivar = np.array(ivar)


# chi-squared:

# chi squared

xp=2

old = np.sum((ref_labels_900[:,xp]-inf_labels_900[:,xp])**2)
opt = np.sum((ref_labels_900[:,xp]-inf_labels_900_opt[:,xp])**2)

print(old,opt)






font = {'weight': 'bold','size': 20}
mpl.rc('font', **font)
fig = plt.figure()


# ax1 = fig.add_axes([0.10, 0.1, 0.85, 1])

plt.scatter(ref_labels_900[:,xp],inf_labels_900_opt[:,xp], marker='x',  c=ivar,
                    vmin=10000, vmax=40000, alpha=0.5,label = "inf_label_opt chi-squared = %.f"%opt)

plt.plot(ref_labels_900[:,xp],inf_labels_900[:,xp],"ro",alpha = 0.3,label = "inf_label_old chi-squared = %.f"%old)

plt.plot(ref_labels_900[:,xp],ref_labels_900[:,xp],"k")

plt.legend(bbox_to_anchor=(0.48, 0.35), loc=2,
                   ncol=1)

plt.xlabel('reference Fe/H', fontsize=30)
plt.ylabel('inferred Fe/H', fontsize=30)
fig.suptitle('1-1 plot of label $Fe/H$', fontsize=24 , ha='center')

#axes = plt.gca()
#axes.set_xlim([-4,4])
#axes.set_ylim([2000,6000])


#cb = plt.colorbar(pl, ax=ax1, orientation='horizontal')
#cb.set_label("Mean inverse variance", fontsize=30)

plt.show()

