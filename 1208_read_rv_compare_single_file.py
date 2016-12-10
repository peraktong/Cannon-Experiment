import pickle
import PyAstronomy
from PyAstronomy import pyasl
import numpy as np
from astropy.io import fits
import matplotlib.pylab as plt
import Cross_correlation_1207 as cc


# load path

pkl_file = open('n_900_path_fits.pkl', 'rb')
path_fits = pickle.load(pkl_file)
pkl_file.close()


pkl_file = open('n_900_path_flux.pkl', 'rb')
path_flux = pickle.load(pkl_file)
pkl_file.close()


N = len(path_fits)
print(N)

velocity = []
mean_fiber_id = []
mean_ivar = []
parameters = []
nor = []
ivar=[]
inf =[]
inf_label= []

for i in range(0,N):
    print("loading star %d"%(i+1))
    star_i = fits.open(path_fits[i])
    velocity.append(star_i[10].data[0])
    mean_fiber_id.append(np.mean(star_i[7].data))
    mean_ivar.append(np.mean(star_i[1].data[0]))
    parameters.append(star_i[4].data[0])

    nor.append(star_i[0].data[0])
    ivar.append(star_i[1].data[0])
    inf.append(star_i[2].data[0])
    inf_label.append(star_i[9].data[0])

    print(star_i[4].data[0].shape)

velocity =np.array(velocity)
mean_fiber_id = np.array(mean_fiber_id)
parameters = np.array(parameters)
nor=np.array(nor)
ivar = np.array(ivar)
inf = np.array(inf)
inf_label = np.array(inf_label)


velocity_compare = []

N = len(path_fits)

for i in range(0,N):
    name = "/Users/caojunzhi/Desktop/Data/velocity_900/" + "velocity" + str(i) + ".pkl"

    pkl_file = open(name, 'rb')
    v_i = pickle.load(pkl_file)
    pkl_file.close()

    velocity_compare.append(v_i)

velocity_compare = np.array(velocity_compare)
print(velocity_compare.shape)

plt.plot(velocity_compare[:,0],"ko",label="CCF")
plt.plot(velocity_compare[:,1],"go",label ="David")
plt.plot(velocity_compare[:,2],"ro",label="Jason")
axes = plt.gca()
#axes.set_xlim([15990,16010])

axes.set_ylim([-400,400])
#axes.set_yticks(np.arange(0.8,2.21,0.35))


plt.show()

output = open("velocity_compare_900.pkl", 'wb')
pickle.dump(velocity_compare, output)
output.close()


