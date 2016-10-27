from numpy import genfromtxt
import matplotlib.pyplot as plt
import random
import matplotlib
import numpy as np
import pickle

# import data
pkl_file = open('delta_chi_500.pkl', 'rb')
delta_chi_500 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('FiberID.pkl', 'rb')
FiberID = pickle.load(pkl_file)
pkl_file.close()

#Plot
font = {'weight': 'bold','size': 25}
matplotlib.rc('font', **font)
fig = plt.figure()

axes = plt.gca()
axes.set_ylim([0,50000])

plt.plot(FiberID,delta_chi_500,"ro")

fig.suptitle('Delta-Chi-Squared', fontsize=30)
plt.xlabel('Fiber ID', fontsize=30)
plt.ylabel('Delta chi-squared', fontsize=30)
plt.show()