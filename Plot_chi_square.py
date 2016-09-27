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


chi_old = np.array([])
p_old = np.array([])

chi_opt = np.array([])
p_opt = np.array([])

for p in range(0,nor[:,0].size):
    nor_i = nor[p,:]
    inf_old_i = inf_old[p,:]
    inf_opt_i = inf_opt[p,:]
    a_old, b_old = chisquare(nor_i, inf_old_i)
    a_opt, b_opt = chisquare(nor_i, inf_opt_i)
    chi_old = np.append(chi_old,a_old)
    #p_old.append(b_old)
    chi_opt = np.append(chi_opt,a_opt)
    #p_old.append(b_opt)

print(chi_opt.shape)
# plot


font = {'weight': 'bold','size': 30}
matplotlib.rc('font', **font)
fig = plt.figure()

plt.plot(chi_old,"r",chi_opt,"g",linewidth =5.0)

fig.suptitle('compare chi-square', fontsize=40)
plt.xlabel('stellar', fontsize=38)
plt.ylabel('chi-square', fontsize=36)

#axes = plt.gca()
#axes.set_xlim([15660,15780])
#axes.set_xlim([16160,16280])
#axes.set_ylim([0.8,1.1])
#axes.set_ylim([0.8,1.1])
plt.show()










