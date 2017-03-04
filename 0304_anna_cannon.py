from TheCannon import apogee
from TheCannon import dataset
import numpy as np
from TheCannon import model


tr_ID, wl, tr_flux, tr_ivar = apogee.load_spectra("/Users/caojunzhi/Downloads/example_DR10/Data")
tr_label = apogee.load_labels("/Users/caojunzhi/Downloads/example_DR10/reference_labels.csv")


test_ID = tr_ID
test_flux = tr_flux
test_ivar = tr_ivar

ds = dataset.Dataset(wl, tr_ID, tr_flux, tr_ivar, tr_label, test_ID, test_flux, test_ivar)
ds.set_label_names(['T_{eff}', '\log g', '[Fe/H]'])
#fig = ds.diagnostics_SNR()
#fig = ds.diagnostics_ref_labels()

ds.ranges = [[371,3192], [3697,5500], [5500,5997], [6461,8255]]
pseudo_tr_flux, pseudo_tr_ivar = ds.continuum_normalize_training_q(q=0.90, delta_lambda=50)

contmask = ds.make_contmask(pseudo_tr_flux, pseudo_tr_ivar, frac=0.07)

ds.set_continuum(contmask)
cont = ds.fit_continuum(3, "sinusoid")

norm_tr_flux, norm_tr_ivar, norm_test_flux, norm_test_ivar = ds.continuum_normalize(cont)

print(np.mean(norm_test_ivar))

ds.tr_flux = norm_tr_flux
ds.tr_ivar = norm_tr_ivar
ds.test_flux = norm_test_flux
ds.test_ivar = norm_test_ivar


md = model.CannonModel(2)
md.fit(ds)

# inf flux:



