import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
import pickle
from os.path import isfile, join
from os import listdir


from astropy.time import Time





pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('suspect_14.pkl', 'rb')
files = pickle.load(pkl_file)
pkl_file.close()

mypath = "/Users/caojunzhi/Desktop/Data/suspect_14"+"/"



def plot_visit_sim_old_customized(name):
    # set font size

    font = {'family': 'normal',
            'weight': 'bold',
            'size': 8}

    matplotlib.rc('font', **font)

    # only choose individual visits:
    star = fits.open(name)

    flux = star[8].data
    ivar = star[9].data
    inf_flux = star[10].data
    parameters =star[0].data
    inf_labels = star[12].data
    inf_flux_sim = star[11].data
    inf_labels_sim = star[6].data[:,0:3]
    parameters_sim = star[6].data[:,3:6]

    parameters_new = star[3].data

    # select visits with 2b<a+c

    a = parameters_new[:,0]
    b = parameters_new[:,1]
    c = parameters_new[:,2]

    mask = []

    for l in range(0,len(a)):
        if 2*b[l]>a[l]+c[l]:
            mask.append(1)

        else:
            mask.append(0)



    N = len(flux[:,0])

    # If you add combined spectra, change 0 to 2

    for i in range(0, N):

        # left
        plt.subplot(N, 2, 2 * i + 1)

        name = name.replace("/Users/caojunzhi/Desktop/Data/suspect_14/", "")
        name = name.replace(".fits", "")


        ## Add notifications for 2b<a+c
        if mask[i]==1:
            data_label = "Data flux"
            color = None
        else:
            #data_label = colored("Data flux of the visit with 2b < a+c","red")
            data_label = "Data flux"
            plt.plot(0,0,"ro",label="Data flux of the visit with 2b < a+c")
            color = "red"


        plt.step(wl, flux[i,:], "k", label=data_label, linewidth=0.7, alpha=0.8)
        plt.errorbar(wl,flux[i,:],yerr=ivar[i,:]**(-0.5),alpha = 0.2,ecolor='k')

        plt.plot(wl, inf_flux[i,:], "g", label="Inferred flux from fitting separately", linewidth=0.7, alpha=0.4)

        """

        plt.plot(wl, inf_flux_sim[i,:], "g", label="Inferred flux from fitting simultaneously Teff=%.2f K logg=%.2f Fe/H =%.2f" % (
        inf_labels_sim[i,0], inf_labels_sim[i,1], inf_labels_sim[i,2]), linewidth=0.7, alpha=0.4)


        """





        plt.xlabel("Wave length $\AA$")
        plt.ylabel("Flux", fontsize=20)
        axes = plt.gca()
        axes.set_xlim([15660,15780])
        # share x axis

        if i == N - 1:
            non = 1
        else:
            axes.set_xticklabels([])

        # axes.set_xlim([16160,16280])
        axes.set_ylim([0.5, 2])
        axes.set_yticks(np.arange(0.5, 1.99, 0.5))
        plt.legend()

        # inf

        plt.subplot(N, 2, 2 * i + 2)


        ## Add notifications for 2b<a+c
        if mask[i]==1:
            data_label = "Data flux"
            color = None
        else:
            #data_label = colored("Data flux of the visit with 2b < a+c","red")
            data_label = "Data flux"
            plt.plot(0,0,"ro",label="Data flux of the visit with 2b < a+c")
            color = "red"



        ## Add inf labels

        plt.step(wl, flux[i,:], "k", linewidth=0.7, alpha=0.8)
        plt.errorbar(wl,flux[i,:],yerr=ivar[i,:]**(-0.5),alpha = 0.2,ecolor='k')

        ai = parameters[i,0]
        bi = parameters[i,1]
        ci = parameters[i,2]

        plt.plot(wl, inf_flux[i,:], "g", label="Velocity from fitting separately a=%.2f b=%.2f c =%.2f $RV_{shift}$=%.2f $m/s$" % (
        ai,bi,ci,(ci-ai)/(ai+bi+ci)*4144.68), linewidth=0.7, alpha=0.4)

        ai = parameters_sim[i,0]
        bi = parameters_sim[i,1]
        ci = parameters_sim[i,2]

        """

        plt.plot(wl, inf_flux_sim[i,:], "g", label="Velocity from fitting simultaneously a=%.2f b=%.2f c =%.2f $RV_{shift}$=%.2f $m/s$" % (
        ai, bi, ci, (ci - ai) / (ai + bi + ci) * 4144.68), linewidth=0.7, alpha=0.4)


        """

        # add absorption lines:

        ai = parameters_new[i,0]
        bi = parameters_new[i,1]
        ci = parameters_new[i,2]

        plt.plot(0,0, "ko",markersize=0.5,
                 label="Velocity from fitting absorption lines a=%.2f b=%.2f c =%.2f $RV_{shift}$=%.2f $m/s$" % (
                     ai, bi, ci, (ci - ai) / (ai + bi + ci) * 4144.68), linewidth=0.7, alpha=0.4)

        plt.xlabel("Wave length $\AA$")
        plt.ylabel("Flux", fontsize=20)
        axes = plt.gca()
        axes.set_xlim([16160,16280])
        # share x axis

        if i == N - 1:
            non = 1
        else:
            axes.set_xticklabels([])

        # axes.set_xlim([16160,16280])
        axes.set_ylim([0.5, 2])
        axes.set_yticks(np.arange(0.5, 1.99, 0.5))
        plt.legend()


    plt.suptitle("The fitting result of star %s"%name, fontsize=20)

    # share x
    plt.subplots_adjust(hspace=.0)

    #plt.show()

    fig = matplotlib.pyplot.gcf()

    # adjust the size based on the number of visit
    visit = len(flux[:,0])

    height = 4+2*visit

    fig.set_size_inches(18.5, height)

    save_path = "/Users/caojunzhi/Downloads/upload_20170309"+"/"+str(name)+".png"
    fig.savefig(save_path, dpi=300)

    plt.close()


# input apstar
def Plot_RV_before_after(name):


    font = {'family': 'normal',

            'size': 24}

    matplotlib.rc('font', **font)


    dat = Table.read(name)
    image = fits.open(name)
    VHELIO = dat["VHELIO"]


    VHELIO = VHELIO.ravel()



    # read HJDs

    # Read BJD for individual visits
    HJD = []
    VHELIO_un = []


    for k in range(0, len(VHELIO)):
        name_k = "HJD" + str(k + 1)
        name_k2 = "VERR"+ str(k + 1)
        #print(name_k)

        HJD_i = image[0].header[name_k]
        HJD.append(HJD_i)
        VHELIO_un.append(image[0].header[name_k2])

    HJD = np.array(HJD)
    VHELIO_un = np.array(VHELIO_un)

    name = name.replace("apStar-r6-", "")

    # Your data from new
    image = fits.open(name)

    """

    a =image[6].data[:,3]
    b = image[6].data[:,4]
    c = image[6].data[:,5]

    ve_sim = (c-a)/(a+b+c)*4.14468

    ve_sim_un = image[7].data


    """

    ve_sim = image[4].data/1000

    a =image[3].data[:,0]
    b = image[3].data[:,1]
    c = image[3].data[:,2]

    mask = []
    mask_data = []
    for m in range(0,len(a)):
        if (2*b[m] > a[m]+c[m]):
            mask.append(0)
            mask_data.append(1)

        else:
            mask.append(1)
            mask_data.append(0)

    mask = np.array(mask,dtype=bool)
    mask_data = np.array(mask_data,dtype=bool)


    plt.subplot(1,1,1)

    std = np.std(VHELIO)
    std_sim = np.std(ve_sim+VHELIO)

    plt.plot(HJD[mask_data],VHELIO[mask_data],"ko",markersize=5,label="DR13 standard deviation = %.3f $km/s$"%std)
    plt.plot(HJD[mask], VHELIO[mask], "ro", markersize=5, label="Visits with 2b<a+c")




    plt.plot(HJD[mask_data], (VHELIO+ve_sim)[mask_data], "bo", markersize=5, label="After shift standard deviation = %.3f $km/s$"%std_sim)

    plt.legend(loc='upper center')

    plt.plot(HJD[mask], (VHELIO + ve_sim)[mask], "ro", markersize=5, label="Visits with 2b<a+c")



    plt.errorbar(HJD[mask_data], VHELIO[mask_data], yerr=VHELIO_un[mask_data], fmt='ko',alpha = 0.5)
    plt.errorbar(HJD[mask_data], (VHELIO + ve_sim)[mask_data], yerr=VHELIO_un[mask_data], fmt='bo',alpha = 0.5)

    plt.errorbar(HJD[mask], VHELIO[mask], yerr=VHELIO_un[mask], fmt='ro',alpha = 0.5)
    plt.errorbar(HJD[mask], (VHELIO + ve_sim)[mask], yerr=VHELIO_un[mask], fmt='ro',alpha = 0.5)



    name = name.replace("/Users/caojunzhi/Desktop/Data/suspect_14/","")
    name = name.replace(".fits","")

    plt.xlabel("HJD")
    plt.ylabel("Radial velocity $Km/s$")
    plt.suptitle("%s"%name)



    # plt.show()
    # save them:

    fig = matplotlib.pyplot.gcf()

    # adjust the size based on the number of visit

    fig.set_size_inches(14.5, 8.5)

    save_path = "/Users/caojunzhi/Downloads/upload_20170309/"+"velocity_"+name+".png"
    fig.savefig(save_path, dpi=300)

    plt.close()


# Write your result into a table.
# The input is Apstar.

def RV_table(name):

    length = len(name)

    id = []
    MJD = []
    HJD = []

    RV = []
    RV_un = []

    ve = []
    ve_un = []

    ve_abs = []
    ve_abs_un = []

    for i in range(0,length):

        dat = Table.read(name[i])
        image = fits.open(name[i])
        VHELIO = dat["VHELIO"]
        MJD_i = dat["MJD"]

        VHELIO = VHELIO.ravel()

        # read HJDs

        # Read BJD for individual visits
        HJD_i = []
        VHELIO_un = []


        for k in range(0, len(VHELIO)):
            name_k = "HJD" + str(k + 1)
            name_k2 = "VERR" + str(k + 1)
            # print(name_k)

            HJD_k = image[0].header[name_k]
            HJD_i.append(HJD_k)
            VHELIO_un.append(image[0].header[name_k2])

        HJD_i = np.array(HJD_i)
        VHELIO_un = np.array(VHELIO_un)

        name_i = name[i].replace("apStar-r6-", "")

        # Your data from new
        image = fits.open(name_i)

        """

        a =image[6].data[:,3]
        b = image[6].data[:,4]
        c = image[6].data[:,5]

        ve_sim = (c-a)/(a+b+c)*4.14468

        ve_sim_un = image[7].data


        """

        ve_i  = image[1].data.ravel()
        ve_un_i = image[2].data.ravel()

        ve_abs_i  = image[4].data.ravel()
        ve_abs_un_i = image[5].data.ravel()

        # add

        id_i = []
        for z in range(0,len(ve_i)):
            x = name_i.replace(".fits","")
            x = x.replace("/Users/caojunzhi/Desktop/Data/suspect_14/","")
            id_i.append(x)

        id_i = np.array(id_i)


        ## append

        id = np.append(id,id_i)
        id = np.array(id)
        id = id.ravel()

        MJD = np.append(MJD,MJD_i)
        MJD = np.array(MJD)
        MJD = MJD.ravel()

        HJD = np.append(HJD,HJD_i)
        HJD = np.array(HJD)
        HJD = HJD.ravel()

        # RV

        RV = np.append(RV,VHELIO)
        RV = np.array(RV)
        RV = RV.ravel()

        RV_un = np.append(RV_un,VHELIO_un)
        RV_un = np.array(RV_un)
        RV_un = RV_un.ravel()

        # VE

        ve = np.append(ve,ve_i)
        ve = np.array(ve)
        ve = ve.ravel()

        ve_un = np.append(ve_un,ve_un_i)
        ve_un = np.array(ve_un)
        ve_un = ve_un.ravel()

        # VE_new


        ve_abs = np.append(ve_abs,ve_abs_i)
        ve_abs = np.array(ve_abs)
        ve_abs = ve_abs.ravel()

        ve_abs_un = np.append(ve_abs_un,ve_abs_un_i)
        ve_abs_un = np.array(ve_abs_un)
        ve_abs_un = ve_abs_un.ravel()


    print(id.shape,HJD.shape,MJD.shape,RV.shape,RV_un.shape,ve.shape,ve_un.shape,ve_abs.shape,ve_abs_un.shape)


    # save them in the header


    path = "/Users/caojunzhi/Downloads/upload_20170309/table_14.fits"

    prihdr = fits.Header()
    prihdr['COMMENT'] = "id MJD HJD VELIO VELIOUN Shift ShiftUN Shiftabs ShiftabsUN"

    prihdu = fits.PrimaryHDU( header=prihdr)

    # Table list


    col1 = fits.Column(name='APOGEEID', format="25A", array=id)
    col2 = fits.Column(name='MJD', format='E', array=MJD)
    col3 = fits.Column(name='HJD', format="E", array=HJD)

    col4 = fits.Column(name='VHELIO', format='E', array=RV)
    col5 = fits.Column(name='VHELIOUN', format="E", array=RV_un)

    col6 = fits.Column(name='RV_shift', format='E', array=ve)
    col7 = fits.Column(name='RV_shift_UN', format="E", array=ve_un)

    col8 = fits.Column(name='RV_shift_abs', format="E", array=ve_abs)
    col9 = fits.Column(name='RV_shift_abs_UN', format="E", array=ve_abs_un)

    cols = fits.ColDefs(
        [col1, col2, col3, col4, col5, col6, col7, col8, col9])

    tbhdu = fits.BinTableHDU.from_columns(cols)

    thdulist = fits.HDUList([prihdu, tbhdu])

    thdulist.writeto(path, clobber=True)



















    #print(HJD)



## Plot fluxes




for i in range(0,len(files)):
    name_i = mypath + files[i]
    name_i = name_i.replace("apStar-r6-", "")
    plot_visit_sim_old_customized(name_i)








## Plot RVs


"""

mypath = "/Users/caojunzhi/Desktop/Data/suspect_14/original_files"

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
onlyfiles = onlyfiles[1:len(onlyfiles)]
print(onlyfiles)

mypath = "/Users/caojunzhi/Desktop/Data/suspect_14/"

for i in range(0,len(onlyfiles)):
    name_i = mypath+onlyfiles[i]
    print(name_i)

    Plot_RV_before_after(name_i)


"""


# construct table:

"""

mypath = "/Users/caojunzhi/Desktop/Data/suspect_14/original_files"

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
onlyfiles = onlyfiles[1:len(onlyfiles)]
print(onlyfiles)

mypath = "/Users/caojunzhi/Desktop/Data/suspect_14/"

index = []
for i in range(0,len(onlyfiles)):
    name_i = mypath+onlyfiles[i]
    print(name_i)

    index.append(name_i)

index = np.array(index)

RV_table(index)


"""








