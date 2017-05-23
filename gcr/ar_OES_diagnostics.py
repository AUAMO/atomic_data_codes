# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 10:06:27 2016

@author: ALEXIS
"""
## %%
import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir('/home/ivan/GDrive/Research/ar_OES')
spect_dir = '/home/ivan/GDrive/Research/ar_OES/DATA'

from OES_lib import retrieve_ALEXIS_data
filename, temp_dens_data, temp_dens_av, wavelengths, spect_array, spect_mean  = retrieve_ALEXIS_data(spect_dir)

temp_av = np.sqrt(np.mean(np.square(temp_dens_data[1,0:15])))
dens_av = np.sqrt(np.mean(np.square(temp_dens_data[2,0:15])))
temp_min = np.max(temp_dens_data[1,0:15])
temp_max = np.min(temp_dens_data[1,0:15])

header_array = np.array([['radius'],['t_e'],['n_e'],['isat'],['di_dv0']])
from OES_lib import get_spect_ratios
ratio_763_826, ratio_826_852, ratio_852_919 = get_spect_ratios(wavelengths, spect_mean)
ratio = ratio_763_826

limit = 16
temp_av = np.sqrt(np.mean(np.square(temp_dens_data[1,0:limit])))
dens_av = np.sqrt(np.mean(np.square(temp_dens_data[2,0:limit])))
temp_min = np.max(temp_dens_data[1,0:limit])
temp_max = np.min(temp_dens_data[1,0:limit])

from OES_lib import get_meta_scale
meta_indx, meta_scale, pec_temps, pec_array = get_meta_scale(temp_av, ratio, 763.5, 826.1)
print(meta_indx)
# print("Average Density = ", dens_av)
# print("Filename: ", filename[44:108] + '.png')
filetag = filename[44:108]

print(filename[44:108] + "|{:6.2f}".format(temp_av)  + "|{:6.2e}".format(dens_av) + "| {:6.2f}".format(ratio) + "|{:5.1f}".format(meta_scale))
## %%
def plot_ALEXIS_data():
    os.chdir('/home/ivan/GDrive/Research/ar_OES/meta_scale_images')
    fig = plt.subplots(2,2,figsize = (9,6))
    plt.subplot(221)
    label = "Average = {:5.2f} (eV)".format(temp_av)
    plt.plot(temp_dens_data[0,:], temp_dens_data[1,:],label = label)
    plt.legend(loc=3)
    plt.subplot(222)
    label = "Average = {:5.2e} m^-3".format(dens_av)
    plt.plot(temp_dens_data[0,:], temp_dens_data[2,:], label = label)
    plt.legend(loc=3)
    plt.subplot(212)
    plt.title
    plt.plot(pec_temps, pec_array, 'rx', label = 'mscale = ' + str(meta_scale))
    plt.axvline(temp_av, color = "orange")
    plt.axhline(ratio, color = 'red')
    plt.legend(loc=4)
    plt.suptitle('Temperature, Density and Metastable Ratio: ALEXIS')
    plt.savefig(filetag + ".png", dpi = 300)
    plt.show(fig)
plot_ALEXIS_data()


###################################    PLOTS  #############################################

# %%
from scipy import stats
os.chdir("/home/ivan/GDrive/Research/ar_OES/DATA")
BC1E1_dat = np.loadtxt('BC1E_1.dat', delimiter='|')
slope, intercept, r_value, p_value, std_err = stats.linregress(BC1E1_dat[:,2], BC1E1_dat[:,3])
lin_fit = BC1E1_dat[:,2] * slope + intercept
r_squared = r_value ** 2

meta_av_E = np.mean(BC1E1_dat[:,3])
label = ("Meta_Av = {:6.2f} ".format(meta_av_E))

fig1 = plt.figure(figsize = (9,6))
plt.title("BC1E_1")
plt.plot(BC1E1_dat[:,2], BC1E1_dat[:,3], 'r*', label = label)
plt.plot(BC1E1_dat[:,2], lin_fit, 'r', label = ("R^2= {:6.2f}".format(r_squared)))
plt.legend(loc=2)
plt.show(fig1)


# %%

BC1E2_dat = np.loadtxt('BC1E_2.dat', delimiter='|')
slope, intercept, r_value, p_value, std_err = stats.linregress(BC1E2_dat[:,2], BC1E2_dat[:,3])
lin_fit = BC1E2_dat[:,2] * slope + intercept
r_squared = r_value ** 2
meta_av_E = np.mean(BC1E2_dat[:,3])
label = ("Meta_Av = {:6.2f} ".format(meta_av_E))

fig2 = plt.figure(figsize = (9,6))
plt.title("BC1E_2")
plt.plot(BC1E2_dat[:,2], BC1E2_dat[:,3], 'r*', label = label)
plt.plot(BC1E2_dat[:,2], lin_fit, 'r', label = ("R^2= {:6.2f}".format(r_squared)))
plt.legend(loc=2)
plt.show(fig2)


# %%

BC1C_dat = np.loadtxt('BC1C.dat', delimiter='|')
slope, intercept, r_value, p_value, std_err = stats.linregress(BC1E2_dat[:,2], BC1E2_dat[:,3])
lin_fit = BC1C_dat[:,2] * slope + intercept
r_squared = r_value ** 2
meta_av_E = np.mean(BC1C_dat[:,3])

fig3 = plt.figure(figsize = (9,6))
plt.title("BC1C")
plt.plot(BC1C_dat[:,2], BC1C_dat[:,3], 'r*', label = "Viewport C")
plt.plot(BC1C_dat[:,2], lin_fit, 'r', label = ("R^2= {:6.2f}".format(r_squared)))
plt.legend(loc=2)
plt.show(fig3)

# %%

BC1A_dat = np.loadtxt('BC1A.dat', delimiter='|')
slope, intercept, r_value, p_value, std_err = stats.linregress(BC1A_dat[:,2], BC1A_dat[:,3])
lin_fit = BC1A_dat[:,2] * slope + intercept
r_squared = r_value ** 2
meta_av_E = np.mean(BC1A_dat[:,3])
label = ("Meta_Av = {:6.2f} ".format(meta_av_E))

fig3 = plt.figure(figsize = (9,6))
plt.title("BC1A")
plt.plot(BC1A_dat[:,2], BC1A_dat[:,3], 'r*', label = label)
plt.plot(BC1A_dat[:,2], lin_fit, 'r', label = ("R^2= {:6.2f}".format(r_squared)))
plt.legend(loc = 2)
plt.show(fig3)

# %%

all_dat = np.loadtxt('ALLPORTS.dat', delimiter='|')
slope, intercept, r_value, p_value, std_err = stats.linregress(all_dat[:,2], all_dat[:,3])
lin_fit = all_dat[:,2] * slope + intercept
r_squared = r_value ** 2
meta_av_E = np.mean(all_dat[:,2])

print(slope, intercept)
fig3 = plt.figure(figsize = (9,6))
plt.title("All Viewports")
plt.plot(all_dat[:,2], all_dat[:,3], 'r*', label = "All Viewports")
plt.plot(all_dat[:,2], lin_fit, 'r', label = ("R^2= {:6.2f}".format(r_squared)))
plt.legend(loc = 2)
plt.show(fig3)


# %%
line_ratio = 5.22
mscale = 3.66275938755 * line_ratio - 12.7236408394 ; print(mscale)

#####################################################################################################
# %%
def ar_temp_finder(line_ratio):
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    os.chdir('/home/ivan/GDrive/Research/ar_OES/DATA')

    def find_nearest(array, value):
        idx = (np.abs(array-value)).argmin()
        return idx, array[idx]

    # Import the .fits file which contains all the pec dat and
    # extract the information.
    hdulist = fits.open('pec_dat_650_950_lowdens.fits')
    pec_temps = hdulist[0].data       # The array of temps in eV
    n_meta = hdulist[1].data       # Array listing the metastable #
    pec_dens = hdulist[2].data       # Density array
    pec_wave = hdulist[3].data       # Wavelengths corresponiding to each PEC
    pec_pec = hdulist[4].data.T     # 3-D array containing all PEC's

    # Declare some variables
    size_meta = n_meta.size          # Number of metastables
    n_pec = pec_wave.size        # Number of PEC's for each metastable
    wl_min = int(pec_wave[0])           # smallest wavelength
    wl_max = pec_wave[n_pec - 1]          # largest wavelength

    # Create dictionaries to hold the invidual meta PEC's and Broadened flux
    pec_meta = dict()
    flux = dict()
    for i in range(size_meta):
        pec_meta[i] = pec_pec[:, :, n_pec*i:n_pec*(i+1)]

    gscale = 2000
    mscale = 3.66275938755 * line_ratio - 12.7236408394
    # mscale = 0.249341860455 * line_ratio + 3.64323291233

    lower, lval = find_nearest(pec_wave, 763.5)
    upper, uval = find_nearest(pec_wave, 826.1)

    pec_ratio_weighted = (gscale * pec_meta[0][:, :, lower] + mscale \
        * pec_meta[1][:, :, lower] + pec_meta[2][:, :, lower]) / (gscale * \
        pec_meta[0][:, :, upper] + mscale * pec_meta[1][:, :, upper] \
        + pec_meta[2][:, :, upper])

    pec_array = np.mean(pec_ratio_weighted, axis=1)
    from scipy.interpolate import interp1d
    pec_fit = interp1d(pec_temps, pec_array, kind = 'cubic')

    indx, val = find_nearest(pec_array, line_ratio)

    pec_x = np.linspace(0.5,20,1000)
    indx2, val2 = find_nearest(pec_fit(pec_x), line_ratio)

    indx, val = find_nearest(pec_array, line_ratio)
    pec_x = np.linspace(0.5,20,1000)
    indx2, val2 = find_nearest(pec_fit(pec_x), line_ratio)

    fig2 = plt.figure(figsize=(9,4), facecolor='white')
    plt.title('Temperature Ratio: mscale = {:6.2f}'.format(mscale))
    plt.xlabel('Temperature (eV)')
    plt.plot(pec_temps, pec_array, 'r*')
    plt.plot(pec_x, pec_fit(pec_x), 'r')
    plt.axhline(line_ratio)
    plt.axvline(pec_x[indx2], label = 'Temp = {:6.2f}'.format(pec_x[indx2]))
    plt.legend(loc = 4)
    plt.xlim(0, 20)
    plt.ylim(0,10)
    plt.savefig("temp_finder_" + str(line_ratio) + "_826_852.png", dpi=300)
    plt.show(fig2)

    print("line ratio = {:6.2f}".format(line_ratio))
    print("metastable ratio (1s5/1s3) = {:6.2f}".format(mscale))
    print("temp (from PEC's) = {:6.2f}".format(pec_temps[indx]))
    print("temp (from interpolated data) = {:6.2f}".format(pec_x[indx2]))

    return(pec_temps[indx], pec_x[indx2], mscale, pec_x, pec_array)
pec_temp, int_temp, mscale, pec_x, pec_array = ar_temp_finder(4.0)

# %%
import numpy as np
all_dat = np.loadtxt('ALLPORTS.dat', delimiter='|')
ratios = all_dat[:,2]
ratios.size
temp_pec = list() ; temp_int = list () ; mfit = list()
for i in range(np.size(ratios)):
    pec_temp, int_temp, mscale, pec_x, pec_array = ar_temp_finder(ratios[i])
    temp_pec.append(pec_temp)
    temp_int.append(int_temp)
    mfit.append(mscale)
temp_pec = np.array(temp_pec)
temp_int = np.array(temp_int)
mfit = np.array(mfit)
test_array = np.array([temp_pec, temp_int, mfit])

# %%
# plt.plot(all_dat[:,0] , mfit, 'g*')
plt.plot(all_dat[:,2], temp_int,  'g*')
# plt.plot(all_dat[:,2], mfit,  'g*')

plt.show()


np.size(all_dat[:,0])
np.size(test_array[:,0].T)








# %%

print(pec_temps)

from scipy.interpolate import interp1d
pec_fit = interp1d(pec_temps, pec_array, kind = 'cubic')
fig3 = plt.figure(figsize = (9,6))
plt.plot(pec_temps, pec_array, 'r*')
pec_x = np.linspace(0.5,20,1000)
plt.plot(pec_x, pec_fit(pec_x), 'r')
plt.show(fig3)
