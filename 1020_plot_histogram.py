import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from copy import deepcopy

import os
import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle
import matplotlib.pyplot as plt
import matplotlib
from astropy.io.fits import getdata

import AnniesLasso_2 as tc


# open data

pkl_file = open('parameters_500.pkl', 'rb')
parameters_500 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('FiberID.pkl', 'rb')
FiberID = pickle.load(pkl_file)
pkl_file.close()

# delete 0,0,0, which are the zero order of
print(parameters_500.shape,FiberID.shape)

colors = ['b', 'g', 'r']
name = ["a","b","c"]

# Plot histogram
fig = plt.figure()
plt.hist(parameters_500,bins=15,color=colors,label=name)
plt.legend(prop={'size': 30})
fig.suptitle('Histogram of parameters a,b and c',fontsize =30)
plt.xlabel('values of a, b and c', fontsize=30)
plt.ylabel('Number of stars', fontsize=30)

plt.show()


print("Plot histogram")

# plot the parameters vs metadata


font = {'weight': 'bold','size': 30}
matplotlib.rc('font', **font)
fig = plt.figure()

plt.plot(FiberID,parameters_500[:,0], "bo", label="a", linewidth=3.0)
plt.plot(FiberID,parameters_500[:,1], "go", label="b", linewidth=3.0)
plt.plot(FiberID,parameters_500[:,2], "ro", label="c", linewidth=3.0)

plt.legend(loc="best")

fig.suptitle('a, b, c against mean fiber number', fontsize=30)
plt.xlabel('fiber number', fontsize=28)
plt.ylabel('values of a, b, c', fontsize=26)
plt.show()

