import os
import numpy as np
from numpy import genfromtxt
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle
import matplotlib.pyplot as plt
import matplotlib

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



"""
# chi squared

xp=0

old = np.sum((tr_label[:,xp]-inf_labels_tr[:,xp])**2)
opt = np.sum((tr_label[:,xp]-inf_labels_tr_opt[:,xp])**2)

print(old,opt)


plt.plot(tr_label[:,xp],inf_labels_tr[:,xp],"ro")
plt.plot(tr_label[:,xp],inf_labels_tr_opt[:,xp],"go")
plt.show()



"""

# chi-squared:

xp=0

old = np.sum((ref_labels_900[:,xp]-inf_labels_900[:,xp])**2)
opt = np.sum((ref_labels_900[:,xp]-inf_labels_900_opt[:,xp])**2)

print(old,opt)


plt.plot(ref_labels_900[:,xp],inf_labels_900[:,xp],"ro")
plt.plot(ref_labels_900[:,xp],inf_labels_900_opt[:,xp],"go")
plt.show()

