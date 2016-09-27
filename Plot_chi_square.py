import numpy as np
from scipy.stats import chisquare
import matplotlib.pyplot as plt
import pickle
import matplotlib

# read data

pkl_file = open('normalized_flux.pkl', 'rb')
nor = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('inferred_spectrum.pkl', 'rb')
inf_old = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('optimized_spectrum.pkl', 'rb')
inf_opt = pickle.load(pkl_file)
pkl_file.close()


"""
nor = nor.ravel()
inf_old = inf_old.ravel()
inf_opt = inf_opt.ravel()

"""
chi_old = np.array([])
p_old = np.array([])

chi_opt = np.array([])
p_opt = np.array([])

for p in range(0,nor[:,0].size):
    nor = nor[p, :]
    inf_old = inf_old[p, :]
    inf_opt = inf_opt[p, :]
    a_old, b_old = chisquare(nor, inf_old)
    a_opt, b_opt = chisquare((nor, inf_opt))
    np.append(chi_old,a_old)
    #p_old.append(b_old)
    np.append(chi_opt,a_opt)
    #p_old.append(b_opt)

# plot
wl = np.ones(nor[:,0].size)

font = {'weight': 'bold','size': 30}
matplotlib.rc('font', **font)
fig = plt.figure()

plt.plot(wl, chi_old,wl,chi_opt,linewidth =5.0)

fig.suptitle('compare chi-square', fontsize=40)
plt.xlabel('stellar', fontsize=38)
plt.ylabel('chi-square', fontsize=36)

#axes = plt.gca()
#axes.set_xlim([15660,15780])
#axes.set_xlim([16160,16280])
#axes.set_ylim([0.8,1.1])
#axes.set_ylim([0.8,1.1])
plt.show()










