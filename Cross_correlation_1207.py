from __future__ import print_function, division
import math
import pickle
import PyAstronomy
from PyAstronomy import pyasl
import numpy as np
import matplotlib.pylab as plt


def log(x):
    return math.log10(x)


def power10(x):
    return 10 ** x


def cross_correlation_velocity(a,b,c):

    pkl_file = open('wl.pkl', 'rb')
    wl = pickle.load(pkl_file)
    pkl_file.close()

    # Create the template
    # tw = wl
    tw = np.linspace(wl[0], wl[-1], 8575*10)
    tf = np.exp(-(tw - 16000.0) ** 2 / (2. * 20 ** 2))

    # choose 20%
    choose = np.arange(0, 1716)
    choose = choose * 5 - 1

    # 0.07 -- 4.2

    # gamma = 0.2086
    gamma = 0.199
    # Create data, which are not that well sampled
    # dw = wl[choose]
    dw = np.linspace(wl[0], wl[-1], 8575*2)

    # df = np.exp(-(dw-5004.007)**2/(2.*0.1**2))
    x = np.exp(-(dw - 16000.0 + gamma) ** 2 / (2. * 20 ** 2))
    y = np.exp(-(dw - 16000.0) ** 2 / (2. * 20 ** 2))
    z = np.exp(-(dw - 16000.0 - gamma) ** 2 / (2. * 20 ** 2))

    df = a * x + b * y + c * z
    print(a, b, c)
    print(4.144 * (c - a) / ((a ** 2 + b ** 2 + c ** 2) ** 0.5))
    print(4.144 * (c - a) / ((a + c - 2 * b) * 4))

    """
    # Plot template and data
    plt.title("Template (blue) and data (red)")
    plt.plot(tw, tf, 'b.-')
    plt.plot(dw, df, 'r.-')
    plt.show()
    """

    # Carry out the cross-correlation.
    # The RV-range is -30 - +30 km/s in steps of 0.6 km/s.
    # The first and last 20 points of the data are skipped.
    print("Doing cross-correlation")
    rv, cc = pyasl.crosscorrRV(dw, df, tw, tf, -5, 5, 5. / 5000., skipedge=20)

    # Find the index of maximum cross-correlation function
    maxind = np.argmax(cc)

    print("Cross-correlation function is maximized at dRV = ", rv[maxind], " km/s")
    if rv[maxind] > 0.0:
        print("  A red-shift with respect to the template")
    else:
        print("  A blue-shift with respect to the template")

    # CCF, David, Me
    return 1000*rv[maxind],4144 * (c - a) / ((a + c - 2 * b) * 4),4.144 * (c - a) / ((a ** 2 + b ** 2 + c ** 2) ** 0.5)



def cross_correlation_velocity_v2(a,b,c):

    pkl_file = open('wl.pkl', 'rb')
    wl = pickle.load(pkl_file)
    pkl_file.close()

    wl_log = list(map(log, wl))
    wl_log = np.array(wl_log)

    # use more pixels

    wl_log_long = np.linspace(wl_log[0], wl_log[-1], 8575 * 50)
    # print(wl_log_long)

    # make a new wave length scale ws: 8575*10 points
    ws = 10 ** (wl_log_long)
    #print("ws shape check")
    #print(ws.shape)

    # divide


    # Create the template


    # The shape of the Gaussian
    # gamma = 0.2086
    gamma = 0.199


    miu = 16000.

    sigma = 10

    # tw = wl
    tw = ws
    # tw = np.linspace(wl[0],wl[-1],8575*10)
    tf = np.exp(-(tw - miu) ** 2 / (2. * sigma ** 2)) + 1

    # choose 20%
    choose = np.arange(0, 1716)
    choose = choose * 5 - 1

    # 0.07 -- 4.2


    # Create data, which are not that well sampled
    # dw = wl[choose]

    # dw = np.linspace(wl[0],wl[-1],8575*2)

    # dw = wl
    dw = ws

    # df = np.exp(-(dw-5004.007)**2/(2.*0.1**2))
    x = np.exp(-(dw - miu + 1 * gamma) ** 2 / (2. * sigma ** 2)) + 1
    y = np.exp(-(dw - miu) ** 2 / (2. * sigma ** 2)) + 1
    z = np.exp(-(dw - miu - 1 * gamma) ** 2 / (2. * sigma ** 2)) + 1

    # remake-abc

    df = a * x + b * y + c * z
    print(a, b, c)
    print(4.144 * (c - a) / ((a ** 2 + b ** 2 + c ** 2) ** 0.5))
    print(4.144 * (c - a) / ((a + c - 2 * b) * 4))

    """
    # Plot template and data
    plt.title("Template (blue) and data (red)")
    plt.plot(tw, tf, 'b.-')
    plt.plot(dw, df, 'r.-')
    plt.show()

    """

    # Carry out the cross-correlation.
    # The RV-range is -30 - +30 km/s in steps of 0.6 km/s.
    # The first and last 20 points of the data are skipped.
    print("Doing cross-correlation")
    rv, cc = pyasl.crosscorrRV(dw, df, tw, tf, -1.5, 1.5, 5. / 5000., skipedge=20)

    # Find the index of maximum cross-correlation function
    maxind = np.argmax(cc)
    print(maxind)

    print("Cross-correlation function is maximized at dRV = ", rv[maxind], " km/s")
    if rv[maxind] > 0.0:
        print("  A red-shift with respect to the template")
    else:
        print("  A blue-shift with respect to the template")

    # CCF, David, Me
    return 1000*rv[maxind],4140 * (c - a) / ((a + c - 2 * b) * 4),4.140 * (c - a) / ((a ** 2 + b ** 2 + c ** 2) ** 0.5)

