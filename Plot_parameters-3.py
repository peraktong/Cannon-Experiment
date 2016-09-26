import matplotlib.pyplot as plt
from numpy import genfromtxt
import matplotlib

# data
para = genfromtxt("fitting parameters-3.csv",delimiter=",")
wl = genfromtxt("wl")

p = 118
font = {'weight': 'bold',
            'size': 30}

matplotlib.rc('font', **font)
fig = plt.figure()
plt.plot(wl, para[0,:]+para[1,:]+para[2,:],linewidth =5.0)
fig.suptitle('fitting parameters', fontsize=40)
plt.xlabel('Wave Length', fontsize=38)
plt.ylabel('a+b+c', fontsize=36)
axes = plt.gca()
axes.set_xlim([15660,15780])
#axes.set_xlim([16160,16280])
axes.set_ylim([0.8,1.1])
plt.show()
