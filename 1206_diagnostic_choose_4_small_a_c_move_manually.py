import os
import numpy as np
from astropy.table import Table
from astropy.io import fits
import AnniesLasso_2 as tc
import matplotlib.pyplot as plt
import pickle

import plot_class_4_star_v2_1107

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


v_jason = (parameters[:,2]-parameters[:,0])*(4144.68)



# find the four stars with biggest v


"""
# Jason

N =len(path_fits)
big_v_4 = (velocity*velocity).argsort()[N-4:N]
big_v_jason = (v_jason*v_jason).argsort()[N-4:N]

model = plot_class_4_star_v3_1201.plot_4()
# model.plot_name_4(path_fits[big_v_4])

model.plot_name_4(path_fits[big_v_jason])


"""
N =len(path_fits)
big_v_4 = (velocity*velocity).argsort()[N-4:N]
big_v_jason = (v_jason*v_jason).argsort()[N-4:N]

#model = plot_class_4_star_v3_1201.plot_4()
# model.plot_name_4(path_fits[big_v_4])
# model.plot_name_4(path_fits[big_v_4])



## Manually adjust the velocity shift

small_4_a_c = abs(parameters[:,2]-parameters[:,0]).argsort()[0:4]

print(parameters[small_4_a_c])
print(velocity[small_4_a_c])


# Let's move the spectrum one pixel left and see what happen to a,b,c

nor_4 =nor[small_4_a_c]
ivar_4 = nor[small_4_a_c]
inf_4 = inf[small_4_a_c]

#load the model

training_set_path = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/reference_labels.csv"
training_set_spectrum_dir = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/Data/"


training_set = Table.read("reference_labels.csv")



def get_pixmask(flux, err):
    bad_flux = ~np.isfinite(flux)
    bad_err = (~np.isfinite(err)) | (err <= 0)
    bad_pixels = bad_err | bad_flux
    return bad_pixels


N = len(training_set)
keep = np.ones(N, dtype=bool)

training_set_flux = []
training_set_ivar = []
training_set_error = []

for i, row in enumerate(training_set):
    image_path = os.path.join(training_set_spectrum_dir, row["ID"])
    if not os.path.exists(image_path):
        print("{}/{} could not be found: {}".format(i + 1, N, image_path))
        keep[i] = False
        continue
    print("{}/{}: {}".format(i + 1, N, image_path))
    image = fits.open(image_path)
    flux = image[1].data
    flux_err = image[2].data
    badpix = get_pixmask(flux, flux_err)
    ivar = 1.0/flux_err**2
    error = flux_err
    # badpix is a array and the length is 8575
    flux[badpix] = 1.0
    ivar[badpix] = 0.0
    training_set_flux.append(flux)
    training_set_ivar.append(ivar)
    training_set_error.append(error)

training_set_flux = np.array(training_set_flux)
training_set_ivar = np.array(training_set_ivar)
training_set_error = np.array(training_set_error)


assert all(keep)

# Construct model.
model = tc.L1RegularizedCannonModel(
    training_set, training_set_flux, training_set_ivar,threads=8)
model.s2 = 0
model.regularization = 0
model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(training_set,
    tc.vectorizer.polynomial.terminator(("Teff_{corr}", "logg_{corr}", "[M/H]_{corr}"), 2))

model.train()

# move the data spectrum


n_pixel = nor_4[0, :].size
n_star = nor_4[:, 0].size
one = np.ones(n_star)

# new method for building matrix
z_data = np.c_[one,nor_4]
z_data = z_data[:,0:n_pixel]

y_data =nor_4

x_data = np.c_[nor_4,one]
x_data = x_data[:,1:n_pixel+1]

# left -x

nor_4 = x_data

opt_4,para_4 = model.fitting_spectrum_parameters_single(nor_4,ivar_4,inf_4)

print(para_4)
print(parameters[small_4_a_c])






model_plot = plot_class_4_star_v2_1107.plot_4()
model_plot.train()

model_plot.plot_velocity(nor_4,ivar_4,ivar_4*0,path_fits[small_4_a_c],inf_4,inf_label=inf_label[small_4_a_c])



"""

#plt.plot(velocity,"ro")
#plt.plot(v_jason,"go")

plt.plot(v_jason/velocity,"ko")

axes = plt.gca()
#axes.set_ylim([-2000,2000])
axes.set_ylim([-5,5])

plt.show()

"""









