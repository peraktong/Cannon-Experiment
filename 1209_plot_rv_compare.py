import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math
from astropy.io import fits
import Plot_compare_velocity_1209 as pm

def log(x):
    return math.log10(x)

def power10(x):
    return 10**x

# load
pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('moments_900.pkl', 'rb')
moment = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('parameters_900.pkl', 'rb')
parameters = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('velocity_compare_900.pkl', 'rb')
velocity_compare = pickle.load(pkl_file)
pkl_file.close()

"""
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





#tw = wl



wl_log = list(map(log,wl))
wl_log = np.array(wl_log)

wl_log_long = np.linspace(wl_log[0],wl_log[-1],8575*50)
#print(wl_log_long)

# make a new wave length scale ws: 8575*10 points
ws = 10**(wl_log_long)

# Create the template

# gamma = 0.2086
gamma = 0.199

sigma =10

# Let's calculate moments of the distribution:

N = len(velocity_compare[:,0])

moment = []
for i in range(0,N):
    a = parameters[i,0]
    b = parameters[i,1]
    c = parameters[i,2]

    tw = ws
    # tw = np.linspace(wl[0],wl[-1],8575*10)
    tf = np.exp(-(tw - 16000.0) ** 2 / (2. * sigma ** 2)) + 1

    # 0.07 -- 4.2


    # Create data, which are not that well sampled
    # dw = wl[choose]

    # dw = np.linspace(wl[0],wl[-1],8575*2)

    # dw = wl
    dw = ws

    # df = np.exp(-(dw-5004.007)**2/(2.*0.1**2))
    x = np.exp(-(dw - 16000.0 + 1 * gamma) ** 2 / (2. * sigma ** 2)) + 1
    y = np.exp(-(dw - 16000.0) ** 2 / (2. * sigma ** 2)) + 1
    z = np.exp(-(dw - 16000.0 - 1 * gamma) ** 2 / (2. * sigma ** 2)) + 1

    # remake-abc

    df = a * x + b * y + c * z

    moment_i =np.mean(abs(df-tf)*10**3)
    print("mean = %f *10^-3 "%moment_i)
    moment.append(moment_i)

moment = np.array(moment)

output = open("moments_900.pkl", 'wb')
pickle.dump(moment, output)
output.close()

"""



# km/s
ccf = velocity_compare[:,0]/1000
david = velocity_compare[:,1]/1000
jason = velocity_compare[:,2]/1000



# save moment and parameters

"""
output = open("parameters_900.pkl", 'wb')
pickle.dump(parameters, output)
output.close()

output = open("moments_900.pkl", 'wb')
pickle.dump(moment, output)
output.close()

"""

########################################
# Plot and check!!
###########################

model = pm.plot_module()

# d-2 4 fake stars:
model.plot_4_fake()

# d-3  big c-a
#model.plot_4_ca()

# d-4  big a-c
#model.plot_4_ac()

# d-5 a+c-2b
#model.plot_4_acb()

# d-6 b=1
#model.plot_4_b_1()



# d-7
# model.plot_ccf_David_moment()

# d-figure 8
# model.plot_ccf_jason_moment()

#d-9
#model.plot_ccf_David_2bca()

#d-10
# model.plot_ccf_jason_2bca()













