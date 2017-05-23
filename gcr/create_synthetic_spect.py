from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

# Import the .fits file which contains all the pec dat and
# extract the information.
hdulist = fits.open('pec_dat_650_950_lowdens.fits')
pec_temps = hdulist[0].data       # The array of temps in eV
n_meta = hdulist[1].data       # Array listing the metastable #
pec_dens = hdulist[2].data       # Density array
pec_wave = hdulist[3].data       # Wavelengths corresponiding to each PEC
pec_pec = hdulist[4].data.T     # 3-D array containing all PEC's

data  = np.loadtxt('temp_5_dens2.5e16.dat')
ALEXIS_wavelengths = data[0]
ALEXIS_spect = data[1]

wavelist = np.zeros((np.size(pec_wave),2))
for i in range(0,np.size(pec_wave)):
    wavelist[i,0] = i
    wavelist[i,1] = pec_wave[i]

templist = np.zeros((np.size(pec_temps),2))
for i in range(0,np.size(pec_temps)):
    templist[i,0] = i
    templist[i,1] = pec_temps[i]

denslist = np.zeros((np.size(pec_dens),2))
for i in range(0,np.size(pec_dens)):
    denslist[i,0] = i
    denslist[i,1] = pec_dens[i]

# Declare some variables
size_meta = n_meta.size          # Number of metastables
n_pec = pec_wave.size        # Number of PEC's for each metastable
wl_min = int(pec_wave[0])           # smallest wavelength
wl_max = pec_wave[n_pec - 1]          # largest wavelength

# Spectrometer Resolution (nm)
lambda_d = 2.0
fwhm = 2.0*(np.log(2.0))**0.5*lambda_d
# These variables will be needed for Doppler Broadening
fwhm = 2.0*(np.log(2.0))**0.5*lambda_d
dwavel = np.amin(fwhm)/20.0
n_wavel = int(((wl_max-wl_min)/dwavel))
flux_lam = np.linspace(wl_min, wl_max, n_wavel)

# Create dictionaries to hold the invidual meta PEC's and Broadened flux
pec_meta = dict()
flux = dict()
for i in range(size_meta):
    pec_meta[i] = pec_pec[:, :, n_pec*i:n_pec*(i+1)]
    flux[i] = np.zeros(n_wavel)

# Choose a temp and density
t_indx = 8
d_indx = 18
print("t_e = {:6.2e} eV".format(templist[t_indx,1]))
print("n_e = {:6.2e} cm^-3 ".format(denslist[d_indx,1]))


# %% Doppler Broadening

# Doppler Broaden the array for the specified t_indx and d_indx

for i in range(n_wavel):
    for j in range(pec_wave.size):
        for k in range(size_meta):
            if flux_lam[i] > (pec_wave[j]-fwhm*5.0) and flux_lam[i] < \
               (pec_wave[j]+fwhm*5.0):
                flux[k][i] = flux[k][i] + pec_meta[k][t_indx, d_indx, j] * \
                (2 * np.pi) ** (-0.5)/lambda_d * np.exp(-abs((flux_lam[i] \
                - pec_wave[j]) / lambda_d) ** 2.0)



# %%


# %%


from astropy.io import fits
from scipy.optimize import curve_fit
from scipy import stats
import os
os.chdir(spect_dir)
import datetime
now = datetime.datetime.now()
import tkinter
from tkinter import filedialog
# we don't want a full GUI, so keep the root window from appearing
root = tkinter.Tk()
root.withdraw()
# show an "Open" dialog box and return the path to the selected file
filename = filedialog.askopenfilename()
hdulist = fits.open(filename)
# hdu_info = hdulist.info()

def spect_get():
    hdu_size = np.size(hdulist)
    scan_lambda = hdulist[hdu_size - 2].data.T
    wavelengths = np.array(scan_lambda[0]).T
    scan_spect = hdulist[hdu_size - 1].data.T

    spect = np.array(scan_spect[0])
    spec_av = np.mean(spect[:, 2000:2047])
    spect = spect - spec_av
    for i in range(1, 10):
        tmp_data = spect
        spect = np.array(scan_spect[i])
        spec_av = np.mean(spect[:, 2000:2047])
        spect = spect - spec_av
        spect_array = np.vstack((tmp_data, spect))
        spect = spect_array

    spect_av = spect.mean(axis=0)
    return (wavelengths, spect_array, spect_av)
wavelengths, spect_array, spect_av = spect_get()


# %% Plot the Doppler broadened flux
fig4 = plt.figure(figsize=(12, 6), facecolor='white')

mscale = 5

gscale = 2000
filetag = "mscale_" + str(mscale) + "_spec.png"

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]
spect_indx, spect_value = find_nearest(ALEXIS_wavelengths, 763.5)

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]
flux_indx, flux_value = find_nearest(flux_lam, 763.5)

lam = flux_lam
y0 = flux[0] * gscale
y1 = flux[1] * mscale
y2 = flux[2] * 1.
flux_sum = (y0 + y1 + y2)
flux_sum = flux_sum / (flux_sum[flux_indx])
plt.plot(flux_lam, flux_sum , 'm.', ms=3, label = 'Theory')

plt.plot(ALEXIS_wavelengths - 0.2, ALEXIS_spect /  ALEXIS_spect[spect_indx],\
         label = 'ALEXIS')
plt.xlabel('Wavelength (nm)', fontsize=14)
plt.ylabel('Line Intensity (Arbitrary Units)', fontsize=14)
plt.title(r'Synthetic Spectrum vs ALEXIS: t_e = 5 eV, n_e = 2.5e16 m^{-3}', fontsize=12 )
plt.ylim(0,3.0)
plt.xlim(760, 900)
plt.legend(loc=2)
plt.savefig(filetag, dpi=300)
plt.show(fig4)
