from numpy import genfromtxt
import matplotlib.pyplot as plt
import random
import matplotlib
import numpy as np
import pickle
from astropy.table import Table
import TheCannon_2
from TheCannon_2 import dataset,apogee
from astropy.table import Table
import numpy as np
import os
from astropy.io import fits
import pickle
import AnniesLasso_2 as tc

import plot_class_4_star_v2_1107








# import data
pkl_file = open('n2_delta_chi_900.pkl', 'rb')
delta_chi_500 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()


pkl_file = open('error_4.pkl', 'rb')
error = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n2_parameters_900.pkl', 'rb')
parameters_500 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n2_FiberID_900.pkl', 'rb')
FiberID = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n2_inf_flux_900.pkl', 'rb')
inf_flux = pickle.load(pkl_file)
pkl_file.close()




pkl_file = open('n2_tr_ID_900.pkl', 'rb')
tr_ID_500 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n2_testing_labels_900.pkl', 'rb')
test_labels = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open("n2_testing_set_900.pkl", 'rb')
testing_set = pickle.load(pkl_file)
pkl_file.close()


# flux,ivar and error

pkl_file = open("n_testing_flux_900.pkl", 'rb')
testing_flux = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open("n_testing_ivar_900.pkl", 'rb')
testing_ivar = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open("n_testing_error_900.pkl", 'rb')
testing_error = pickle.load(pkl_file)
pkl_file.close()


## check error
"""
plt.plot(testing_error[1,:],"r")

axes = plt.gca()

axes.set_ylim([-0.01,0.01])


plt.show()

"""





#print(choose_800_final.shape)
#print(choose_800_final.shape,FiberID.shape)


# check shape
#print(parameters_500.shape,delta_chi_500.shape)

#### find the biggest four delta-chi and smallest four
N_star = len(delta_chi_500)


delta_chi_500 = np.array(delta_chi_500)

big_four = delta_chi_500.argsort()[N_star-4:N_star]
small_four = delta_chi_500.argsort()[0:4]

#print(big_four,small_four)
#print(delta_chi_500[big_four],delta_chi_500[small_four])

big_four_set = testing_set[big_four]
small_four_set = testing_set[small_four]

# print(delta_chi_500[big_four],delta_chi_500[small_four])

# find parameters closest to b =-1


a_500 = parameters_500[:,0]
a_500 = np.array(a_500)

b_500 = parameters_500[:,1]
b_500 = np.array(b_500)

c_500 = parameters_500[:,2]
c_500 = np.array(c_500)

plus = np.ones(N_star)
minus = np.ones(N_star)*(-1)

#print(plus,minus)


close_4_1 =((b_500-plus)**2).argsort()[0:4]
close_4_minus_1 =((b_500-minus)**2).argsort()[0:4]
biggest_a_c_b = (a_500+c_500-b_500).argsort()[N_star-4:N_star]

# find flux and save them
"""
print(close_4_1,close_4_minus_1)
print(b_500[close_4_1],b_500[close_4_minus_1])
print(tr_ID_500[close_4_1])
print(tr_ID_500[close_4_minus_1])
print(parameters_500[close_4_1],parameters_500[close_4_minus_1])


print(close_4_1,close_4_minus_1)
print(b_500[close_4_1],b_500[close_4_minus_1])
print(tr_ID_500[close_4_1])
print(tr_ID_500[close_4_minus_1])
print(parameters_500[close_4_1],parameters_500[close_4_minus_1])

"""




#print(testing_set[close_4_1]["APOGEE_ID"],testing_set[close_4_minus_1]["APOGEE_ID"])

close_1_set = testing_set[close_4_1]
close_minus_1_set = testing_set[close_4_minus_1]
biggest_a_c_b_set = testing_set[biggest_a_c_b]
## name

"""
print(big_four_set,small_four_set)

print(close_1_set,close_minus_1_set)

"""

# also choose four stars with biggest |a|-|c| and |c|-|a|

# a,c

"""
big_a_c = (np.absolute(parameters_500[:,0])-np.absolute(parameters_500[:,2])).argsort()[N_star-4:N_star]
big_a_c_set = testing_set[big_a_c]

big_c_a = (np.absolute(parameters_500[:,2])-np.absolute(parameters_500[:,0])).argsort()[N_star-4:N_star]
big_c_a_set = testing_set[big_c_a]

"""

big_a_c = (parameters_500[:,0]-parameters_500[:,2]).argsort()[N_star-4:N_star]
big_a_c_set = testing_set[big_a_c]

big_c_a = (parameters_500[:,2]-parameters_500[:,0]).argsort()[N_star-4:N_star]
big_c_a_set = testing_set[big_c_a]



print(max(parameters_500[:,0]-parameters_500[:,2]))
print(test_labels[big_a_c])






"""
print(big_a_c,big_c_a)
print(testing_set[big_a_c],testing_set[big_c_a])
print(big_a_c_set,big_c_a_set)

print(parameters_500[big_a_c],parameters_500[big_c_a])

print(testing_set[big_a_c]["APOGEE_ID"],testing_set[big_c_a]["APOGEE_ID"])


"""


"""
plt.plot(test_labels[:,1],"ro")
plt.show()

"""





"""
#save

output = open('big_four_delta_set.pkl', 'wb')
pickle.dump(big_four_set, output)
output.close()

output = open('small_four_delta_set.pkl', 'wb')
pickle.dump(small_four_set, output)
output.close()

output = open('close_1_set.pkl', 'wb')
pickle.dump(close_1_set, output)
output.close()

output = open('close_minus_1_set.pkl', 'wb')
pickle.dump(close_minus_1_set, output)
output.close()

output = open('biggest_a_c_b_set.pkl', 'wb')
pickle.dump(biggest_a_c_b_set, output)
output.close()

output = open('big_a_c_set.pkl', 'wb')
pickle.dump(big_a_c_set, output)
output.close()

output = open('big_c_a_set.pkl', 'wb')
pickle.dump(big_c_a_set, output)
output.close()



"""

##################
# plot:
###################



#model = plot_class_4_star_v2_1107.plot_4()
#model.train()

print("checking")
print(tr_ID_500[big_a_c])
print(parameters_500[big_a_c])




"""
# plot a>c

model.plot_diag(flux=testing_flux[big_a_c],
           ivar=testing_ivar[big_a_c],
           error=testing_error[big_a_c],
           star_name=testing_set[big_a_c]["APOGEE_ID"],test_labels=test_labels[big_a_c],
           para_4=parameters_500[big_a_c])
"""




"""

model = plot_class_4_star_v2_1107.plot_4()
model.train()
# plot big delta-chi

model.plot(flux=testing_flux[big_four],
           ivar=testing_ivar[big_four],
           error=testing_error[big_four],
           star_name=testing_set[big_four]["APOGEE_ID"],
           test_label=test_labels[big_four])

# plot small delta-chi

model.plot(flux=testing_flux[small_four],
           ivar=testing_ivar[small_four],
           error=testing_error[small_four],
           star_name=testing_set[small_four]["APOGEE_ID"],
           test_label=test_labels[small_four])

# plot b = 1

model.plot(flux=testing_flux[close_4_1],
           ivar=testing_ivar[close_4_1],
           error=testing_error[close_4_1],
           star_name=testing_set[close_4_1]["APOGEE_ID"],
           test_label=test_labels[close_4_1])

# plot biggest a+c-b

model.plot(flux=testing_flux[biggest_a_c_b],
           ivar=testing_ivar[biggest_a_c_b],
           error=testing_error[biggest_a_c_b],
           star_name=testing_set[biggest_a_c_b]["APOGEE_ID"],
           test_label=test_labels[biggest_a_c_b])

# plot a>c

model.plot(flux=testing_flux[big_a_c],
           ivar=testing_ivar[big_a_c],
           error=testing_error[big_a_c],
           star_name=testing_set[big_a_c]["APOGEE_ID"],
           test_label=test_labels[big_a_c])

# plot c>a

model.plot(flux=testing_flux[big_c_a],
           ivar=testing_ivar[big_c_a],
           error=testing_error[big_c_a],
           star_name=testing_set[big_c_a]["APOGEE_ID"],
           test_label=test_labels[big_c_a])

"""

# plot





font = {'weight': 'bold', 'size': 13}
matplotlib.rc('font', **font)
fig = plt.figure()

p1 = 0
trans = 0.5
l1 = "star A"
l2 = "star B"
l3 = "star C"
l4 = "star D"

nor = testing_flux[big_a_c]
inf_old = inf_flux[big_a_c]

n_pixel = nor[0, :].size

one = np.ones(4)

# new method for building matrix
z_data = np.c_[one,inf_old]
z_data = z_data[:,0:n_pixel]

y_data =inf_old

x_data = np.c_[inf_old,one]
x_data = x_data[:,1:n_pixel+1]


# set a,b and c manually
a = 2
b = -2
c = 1


"""
parameters_4 = np.array([[0.28,0.815,-0.1],[0.287,0.717,-0.009],
                         [0.152,0.974,-0.127],[0.216,0.811,-0.031]])
"""


parameters_4 = [[ 0.053 , 0.885  , 0.058],
 [ 0.045 , 1.051, -0.097],
 [-0.127 , 1.039 , 0.088],
 [ 2.00 , -1.50, 0.50]]

parameters_4 = np.array(parameters_4)


inf_opt = np.ones(8575)
for i in range(0,4):
    inf_opt_i = parameters_4[i,0]*x_data[i,:]+parameters_4[i,1]*y_data[i,:]+parameters_4[i,2]*z_data[i,:]
    inf_opt = np.vstack((inf_opt,inf_opt_i))

inf_opt = np.array(inf_opt)
inf_opt = inf_opt[1:5,:]
print(inf_opt.shape)







"""
# ax1

plt.step(wl, nor[p1, :], "k", label=l1, alpha=1, linewidth=1.5)
plt.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
plt.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
plt.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

plt.legend(bbox_to_anchor=(0, 0.65), loc=3,
           ncol=1)

axes = plt.gca()
axes.set_xlim([15660,15780])

#axes.set_xlim([16160,16280])
axes.set_ylim([0.8,1.21])
axes.set_yticks(np.arange(0.8,1.21,0.1))

fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")


plt.show()
# f.savefig(delta_para_1+".png",dpi=400,papertype = "a4")

"""


# plot-opt
trans = 0.5
p1 = 0
p2 = 1
p3 = 2
p4 = 3



# test label



# delta-chi and parameters a, b and c.
delta_para_1 = "a=%.3f b=%.3f c=%.3f" % (
    parameters_4[0, 0], parameters_4[0, 1], parameters_4[0, 2])

delta_para_2 = "a=%.3f b=%.3f c=%.3f" % (
    parameters_4[1, 0], parameters_4[1, 1], parameters_4[1, 2])

delta_para_3 = "a=%.3f b=%.3f c=%.3f" % (
    parameters_4[2, 0], parameters_4[2, 1], parameters_4[2, 2])

delta_para_4 = "a=%.3f b=%.3f c=%.3f" % (
    parameters_4[3, 0], parameters_4[3, 1], parameters_4[3, 2])

print(parameters_4)

## Let's plot

font = {'weight': 'bold', 'size': 13}
matplotlib.rc('font', **font)
fig = plt.figure()

f, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = \
    plt.subplots(4, 2, sharex='col', sharey='row')

# ax1

ax1.step(wl, nor[p1, :], "k", label=l1, alpha=1, linewidth=1.5)
ax1.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
ax1.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
ax1.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

ax1.legend(bbox_to_anchor=(0, 0.65), loc=3,
           ncol=1)
ax1.set_xlim([15660, 15780])
ax1.set_ylim([0.8, 1.21])
ax1.set_yticks(np.arange(0.8, 1.21, 0.1))

fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

# ax2

ax2.step(wl, nor[p1, :], "k", label=delta_para_1, alpha=1, linewidth=1.5)
ax2.errorbar(wl, nor[p1, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
ax2.plot(wl, inf_opt[p1, :], "r", alpha=trans, linewidth=1.5)
ax2.plot(wl, inf_old[p1, :], "g", alpha=trans, linewidth=1.5)

ax2.legend(bbox_to_anchor=(0, 0.65), loc=3,
           ncol=1)
ax2.set_xlim([16160, 16280])
ax2.set_ylim([0.8, 1.21])
ax2.set_yticks(np.arange(0.8, 1.21, 0.1))

fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

# ax3

ax3.step(wl, nor[p2, :], "k", label=l2, alpha=1, linewidth=1.5)
ax3.errorbar(wl, nor[p2, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
ax3.plot(wl, inf_opt[p2, :], "r", alpha=trans, linewidth=1.5)
ax3.plot(wl, inf_old[p2, :], "g", alpha=trans, linewidth=1.5)

ax3.legend(bbox_to_anchor=(0, 0.65), loc=3,
           ncol=1)
ax3.set_xlim([15660, 15780])
ax3.set_ylim([0.8, 1.2])
ax3.set_yticks(np.arange(0.8, 1.2, 0.1))

fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

# ax4

ax4.step(wl, nor[p2, :], "k", label=delta_para_2, alpha=1, linewidth=1.5)
ax4.errorbar(wl, nor[p2, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
ax4.plot(wl, inf_opt[p2, :], "r", alpha=trans, linewidth=1.5)
ax4.plot(wl, inf_old[p2, :], "g", alpha=trans, linewidth=1.5)

ax4.legend(bbox_to_anchor=(0, 0.65), loc=3,
           ncol=1)
ax4.set_xlim([16160, 16280])
ax4.set_ylim([0.8, 1.2])
ax4.set_yticks(np.arange(0.8, 1.2, 0.1))

fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

# ax5

ax5.step(wl, nor[p3, :], "k", label=l3, alpha=1, linewidth=1.5)
ax5.errorbar(wl, nor[p3, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
ax5.plot(wl, inf_opt[p3, :], "r", alpha=trans, linewidth=1.5)
ax5.plot(wl, inf_old[p3, :], "g", alpha=trans, linewidth=1.5)

ax5.legend(bbox_to_anchor=(0, 0.65), loc=3,
           ncol=1)
ax5.set_xlim([15660, 15780])
ax5.set_ylim([0.8, 1.2])
ax5.set_yticks(np.arange(0.8, 1.2, 0.1))

fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

# ax6

ax6.step(wl, nor[p3, :], "k", label=delta_para_3, alpha=1, linewidth=1.5)
ax6.errorbar(wl, nor[p3, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
ax6.plot(wl, inf_opt[p3, :], "r", alpha=trans, linewidth=1.5)
ax6.plot(wl, inf_old[p3, :], "g", alpha=trans, linewidth=1.5)

ax6.legend(bbox_to_anchor=(0, 0.65), loc=3,
           ncol=1)
ax6.set_xlim([16160, 16280])
ax6.set_ylim([0.8, 1.2])
ax6.set_yticks(np.arange(0.8, 1.2, 0.1))

fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

# ax7

ax7.step(wl, nor[p4, :], "k", label=l4, alpha=1, linewidth=1.5)
ax7.errorbar(wl, nor[p4, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
ax7.plot(wl, inf_opt[p4, :], "r", alpha=trans, linewidth=1.5)
ax7.plot(wl, inf_old[p4, :], "g", alpha=trans, linewidth=1.5)

ax7.legend(bbox_to_anchor=(0, 0.65), loc=3,
           ncol=1)
ax7.set_xlim([15660, 15780])
ax7.set_ylim([0.8, 1.2])
ax7.set_yticks(np.arange(0.8, 1.2, 0.1))

fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

# ax8

ax8.step(wl, nor[p4, :], "k", label=delta_para_4, alpha=1, linewidth=1.5)
ax8.errorbar(wl, nor[p4, :], ecolor='k', alpha=trans / 10, capthick=0.2, yerr=error[p1, :])
ax8.plot(wl, inf_opt[p4, :], "r", alpha=trans, linewidth=1.5)
ax8.plot(wl, inf_old[p4, :], "g", alpha=trans, linewidth=1.5)

ax8.legend(bbox_to_anchor=(0, 0.65), loc=3,
           ncol=1)
ax8.set_xlim([16160, 16280])
ax8.set_ylim([0.8, 1.2])
ax8.set_yticks(np.arange(0.8, 1.2, 0.1))

fig.suptitle('Comparison of the spectrum', fontsize=30, horizontalalignment="center")

ax3.set_ylabel('Normalized flux', fontsize=30)
ax3.yaxis.set_label_coords(-0.05, 0)

ax7.set_xlabel('Wave length $\AA$ ', fontsize=25)
ax8.set_xlabel('Wave length $\AA$ ', fontsize=25)

f.subplots_adjust(hspace=0)

# Don't use this
# plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

plt.show()
# f.savefig(delta_para_1+".png",dpi=400,papertype = "a4")





















