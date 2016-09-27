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
