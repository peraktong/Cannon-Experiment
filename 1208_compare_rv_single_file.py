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

# velocity compare  CFF,David,Jason
#N =605
a=0
b=100

for i in range(a,b):

    print('\033[1;31mdoing star %d of %d\033[1;m' % (i+1,N))
    # CC, David, Jason
    v1,v2,v3 = cc.cross_correlation_velocity_v2(parameters[i,0],parameters[i,1],parameters[i,2])

    velocity_i = np.array([v1,v2,v3])
    name = "/Users/caojunzhi/Desktop/Data/velocity_900/"+"velocity"+str(i)+".pkl"

    output = open(name, 'wb')
    pickle.dump(velocity_i, output)
    output.close()

