from __future__ import print_function, division
import math
import pickle
import PyAstronomy
import matplotlib

from PyAstronomy import pyasl
import numpy as np
import matplotlib.pylab as plt


import plot_class_4_star_v2_1107

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()

def log(x):
    return math.log10(x)

def power10(x):
    return 10**x

wl_log = list(map(log,wl))
wl_log = np.array(wl_log)

wl_log_long = np.linspace(wl_log[0],wl_log[-1],8575*50)
#print(wl_log_long)

# make a new wave length scale ws: 8575*10 points
ws = 10**(wl_log_long)
print("ws shape check")
print(ws.shape)

# divide

# May use this
"""
wl_log_long = np.linspace(wl_log[0],wl_log[-1],8575*10)
wl_log_long_s =np.linspace(wl_log[4000],wl_log[4575],8575*2)
wl_long = list(map(power10,wl_log_long))
wl_long = np.array(wl_log)

wl_long_s = list(map(power10,wl_log_long_s))
wl_long_s = np.array(wl_long_s)
print(wl_long_s[0],wl_long_s[-1])
"""



# Create the template

# gamma = 0.2086
gamma = 0.199

sigma =10

#tw = wl
tw =ws
#tw = np.linspace(wl[0],wl[-1],8575*10)
tf = np.exp(-(tw-16000.0)**2/(2.*sigma**2))+1

# choose 20%
choose = np.arange(0,1716)
choose = choose*5-1

# 0.07 -- 4.2


# Create data, which are not that well sampled
#dw = wl[choose]

#dw = np.linspace(wl[0],wl[-1],8575*2)

#dw = wl
dw =ws

# df = np.exp(-(dw-5004.007)**2/(2.*0.1**2))
x = np.exp(-(dw-16000.0+1*gamma)**2/(2.*sigma**2))+1
y = np.exp(-(dw-16000.0)**2/(2.*sigma**2))+1
z = np.exp(-(dw-16000.0-1*gamma)**2/(2.*sigma**2))+1

# remake-abc

a = 0.13
b = 0.8
c = 0.07

df = a*x+b*y+c*z
print(a,b,c)

simple = 4.144*(c-a)/(a+b+c)
jason = 4.144*(c-a)/((a**2+b**2+c**2)**0.5)
david = 4.144*(a-c)/((a+c-2*b)*2)

print(simple)
print(jason)
print(david)

N = 8575*25

mean = np.mean((df-tf)*10**8)
print("mean = %f *10^-8 "%mean)






# Carry out the cross-correlation.
# The RV-range is -30 - +30 km/s in steps of 0.6 km/s.
# The first and last 20 points of the data are skipped.
print("Doing cross-correlation")
rv, cc = pyasl.crosscorrRV(dw, df, tw, tf, -1.5, 1.5, 5./5000., skipedge=20)

# Find the index of maximum cross-correlation function
maxind = np.argmax(cc)
print(maxind)

print("Cross-correlation function is maximized at dRV = ", rv[maxind], " km/s")
if rv[maxind] > 0.0:
  print("  A red-shift with respect to the template")
else:
  print("  A blue-shift with respect to the template")


# Plot template and data

velocity = "v_David = %.3f km/s, v_Jason = %.3f km/s, v_ccf = %.3f"%(david,jason,rv[maxind])
parameters = "a=%.3f b=%.3f c=%.3f moments = %.3f *10^-8 "%(a,b,c,mean)

font = {'weight': 'bold','size': 15}
matplotlib.rc('font', **font)
fig = plt.figure()
fig.suptitle("Template (blue) and data (red)", fontsize=24)
plt.plot(tw, tf, 'b',label=velocity,alpha=0.5)
plt.plot(dw, df, 'r',label=parameters,alpha=0.5)

plt.xlabel("Wave length $\AA$")
plt.ylabel("flux")
plt.legend(loc="best")

axes = plt.gca()
axes.set_xlim([15990,16010])

axes.set_ylim([0.8,2.21])
axes.set_yticks(np.arange(0.8,2.21,0.35))

plt.show()


"""
plt.plot(rv, cc, 'bp-')
plt.plot(rv[maxind], cc[maxind], 'ro')
plt.show()

"""