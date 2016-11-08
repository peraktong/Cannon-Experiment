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
pkl_file = open('n_delta_chi_900.pkl', 'rb')
delta_chi_500 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n_parameters_900.pkl', 'rb')
parameters_500 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n_FiberID_900.pkl', 'rb')
FiberID = pickle.load(pkl_file)
pkl_file.close()


pkl_file = open('n_tr_ID_900.pkl', 'rb')
tr_ID_500 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('n_testing_labels_900.pkl', 'rb')
test_labels = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open("testing_set_900.pkl", 'rb')
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

"""

print(close_4_1,close_4_minus_1)
print(b_500[close_4_1],b_500[close_4_minus_1])
print(tr_ID_500[close_4_1])
print(tr_ID_500[close_4_minus_1])
print(parameters_500[close_4_1],parameters_500[close_4_minus_1])


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



##################
# plot:
###################
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



















