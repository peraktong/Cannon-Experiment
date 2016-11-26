import pickle
import matplotlib
import matplotlib.pyplot as plt


# read data

pkl_file = open('inf_labels_deviation_teff_900.pkl', 'rb')
inf_labels_deviation_teff_900 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('inf_labels_deviation_logg_900.pkl', 'rb')
inf_labels_deviation_logg_900 = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('inf_labels_deviation_fe_900.pkl', 'rb')
inf_labels_deviation_fe_900 = pickle.load(pkl_file)
pkl_file.close()


# plot histogram


"""
font = {'weight': 'bold','size': 30}
matplotlib.rc('font', **font)
fig = plt.figure()



colors = ['b', 'g', 'r']
name = ["Delta Teff","Delta logg","Delta Fe/H"]
print(inf_labels_deviation_teff_900.shape,inf_labels_deviation_logg_900.shape,inf_labels_deviation_fe_900.shape)

# Plot histogram
fig = plt.figure()
plt.hist(inf_labels_deviation_teff_900,bins=15,color=colors[0],label=name[0])
plt.legend(prop={'size': 30})
fig.suptitle('Histogram for the deviations of the inferred Teff',fontsize =30)
plt.xlabel('values of Delta Teff', fontsize=30)
plt.ylabel('Number of stars', fontsize=30)
plt.show()


colors = ['b', 'g', 'r']
name = ["Delta Teff","Delta logg","Delta Fe/H"]

# Plot histogram
fig = plt.figure()
plt.hist(inf_labels_deviation_logg_900,bins=15,color=colors[1],label=name[1])
plt.legend(prop={'size': 30})
fig.suptitle('Histogram for the deviations of the inferred logg',fontsize =30)
plt.xlabel('values of Delta logg', fontsize=30)
plt.ylabel('Number of stars', fontsize=30)
plt.show()
"""


colors = ['b', 'g', 'r']
name = ["Delta Teff","Delta logg","Delta Fe/H"]

# Plot histogram
fig = plt.figure()
plt.hist(inf_labels_deviation_fe_900,bins=15,color=colors[2],label=name[2])
plt.legend(prop={'size': 30})
fig.suptitle('Histogram for the deviations of the inferred Fe/H',fontsize =30)
plt.xlabel('values of Delta Fe/H', fontsize=30)
plt.ylabel('Number of stars', fontsize=30)
plt.show()



