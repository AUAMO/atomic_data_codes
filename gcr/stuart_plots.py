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

# Choose upper and lower indices for the lines.

# %% TEMP RATIO PLOT 1
gscale = 2000
mscale = 5
upper = 35
lower = 30



pec_ratio_weighted = (gscale * pec_meta[0][:, :, lower] + mscale \
    * pec_meta[1][:, :, lower] + pec_meta[2][:, :, lower]) / (gscale * \
    pec_meta[0][:, :, upper] + mscale * pec_meta[1][:, :, upper] \
    + pec_meta[2][:, :, upper])

print("weighted ratio average = " , np.mean(pec_ratio_weighted[6:12]))

fig2 = plt.figure(figsize=(8, 6), facecolor='white')
plt.title('Temperature Ratio: mscale = ' + str(mscale) + "; gscale = " +str(gscale))
plt.xlabel('Temperature (eV)')
label1 = 'weighted ratio = '+str(pec_wave[lower]) + \
    '/' + str(pec_wave[upper])
plt.plot(pec_temps, pec_ratio_weighted[:,0], 'rx', label = label1 )
plt.plot(pec_temps, pec_ratio_weighted[:,1:19], 'rx')
plt.axhline(1.0)
plt.legend(loc=4)
plt.xlim(0, 20)
plt.ylim(0,4)

plt.savefig("temp_ratio_826_852.png", dpi=300)
plt.show(fig2)

# %% TEMP RATIO PLOT 1
gscale = 2000
mscale = 7
upper = 30
lower = 15

pec_ratio_weighted = (gscale * pec_meta[0][:, :, lower] + mscale \
    * pec_meta[1][:, :, lower] + pec_meta[2][:, :, lower]) / (gscale * \
    pec_meta[0][:, :, upper] + mscale * pec_meta[1][:, :, upper] \
    + pec_meta[2][:, :, upper])

print("weighted ratio average = " , np.mean(pec_ratio_weighted[3:]))

fig2 = plt.figure(figsize=(8, 6), facecolor='white')
plt.title('Temperature Ratio: mscale = ' + str(mscale) + "; gscale = " +str(gscale))
plt.xlabel('Temperature (eV)')

label1 = 'weighted ratio = '+str(pec_wave[lower]) + \
    '/' + str(pec_wave[upper])
plt.plot(pec_temps, pec_ratio_weighted[:,0], 'rx', label = label1 )
plt.plot(pec_temps, pec_ratio_weighted[:,1:19], 'rx')
plt.legend(loc=4)
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]
rindx, rval = find_nearest(pec_ratio_weighted[:,0], 5.49)
rindx, rval

plt.axvline(pec_temps[rindx], color = 'green')
plt.axvline(4.7, color = "cyan")
plt.axvline(7.1, color = "cyan")
plt.axvline(5.92, color = "orange")

plt.axhline(5.49)
plt.xlim(0, 10)
plt.ylim(0,10)

plt.savefig("temp_ratio_826_852.png", dpi=300)
plt.show(fig2)


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]
rindx, rval = find_nearest(pec_ratio_weighted[:,0], 5.49)
rindx, rval
# %%
fig5 = plt.figure(figsize=(12, 6), facecolor='white')

mscale = 7
gscale = 2000



for i in range(n_pec):
    fig = plt.figure(figsize=(12, 6), facecolor='white')
    pec_indx = i
    pec_ratio_1s3_1 = pec_meta[2][:,0,pec_indx] / pec_meta[2][:,0,pec_indx]
    pec_ratio_1s5_2 = pec_meta[1][:,0,pec_indx] / pec_meta[2][:,0,pec_indx]
    pec_ratio_grd_3 = pec_meta[0][:,0,pec_indx] / pec_meta[2][:,0,pec_indx]

    pec_ratio_1s3 = pec_meta[2][:,1:19,pec_indx] / pec_meta[2][:,1:19,pec_indx]
    pec_ratio_1s5 = pec_meta[1][:,1:19,pec_indx] / pec_meta[2][:,1:19,pec_indx]
    pec_ratio_grd = pec_meta[0][:,1:19,pec_indx] / pec_meta[2][:,1:19,pec_indx]

    plt.semilogy(pec_temps, pec_ratio_1s3_1, 'r', label = '1s3 limit' )
    plt.semilogy(pec_temps, pec_ratio_1s5_2, 'gx', label = 'ratio: 1s5 / 1s3')
    plt.semilogy(pec_temps, pec_ratio_grd_3, 'c+', label = 'ratio: ground / 1s3')

    plt.semilogy(pec_temps, pec_ratio_1s3, 'r' )
    plt.semilogy(pec_temps, pec_ratio_1s5, 'gx')
    plt.semilogy(pec_temps, pec_ratio_grd, 'c+')
    plt.ylim(1e-5, 1e2)
    plt.legend(loc=4)
    plt.title("PEC metastable ratio plots functions of temperature and density")
    filestring = str(wavelist[pec_indx,1])
    plt.title("PEC metastable ratio plots functions of temperature and density: " \
              + filestring + "nm. PEC index = " + str(pec_indx))

    filename = "pec_ratio" + str(i) + "_wl=" + filestring + ".png"
    plt.savefig(filename, dpi=150)
    plt.show(fig)

# %%


for i in range(n_pec):
    fig = plt.figure(figsize=(12, 6), facecolor='white')
    pec_indx = i
    pec_1s3_1 = pec_meta[2][:,0,pec_indx]
    pec_1s5_2 = pec_meta[1][:,0,pec_indx]
    pec_grd_3 = pec_meta[0][:,0,pec_indx]

    pec_1s3 = pec_meta[2][:,1:19,pec_indx]
    pec_1s5 = pec_meta[1][:,1:19,pec_indx]
    pec_grd = pec_meta[0][:,1:19,pec_indx]

    plt.semilogy(pec_temps, pec_1s3_1, 'rx', label = '1s3')
    plt.semilogy(pec_temps, pec_1s5_2, 'gx', label = '1s5')
    plt.semilogy(pec_temps, pec_grd_3, 'c+', label = 'ground')

    plt.semilogy(pec_temps, pec_1s3, 'rx' )
    plt.semilogy(pec_temps, pec_1s5, 'gx')
    plt.semilogy(pec_temps, pec_grd, 'c+')
    # plt.ylim(1e-5, 1e2)
    plt.legend(loc=4)
    plt.title("PEC metastable ratio plots functions of temperature and density")
    filestring = str(wavelist[pec_indx,1])
    plt.title("PEC metastable ratio plots functions of temperature and density: " \
              + filestring + "nm. PEC index = " + str(pec_indx))
    filename = "pec_" + str(i) + "_wl=" + filestring + ".png"
    plt.savefig(filename, dpi=150)
    plt.show(fig)

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

# %% Plot the Doppler broadened flux
fig4 = plt.figure(figsize=(12, 6), facecolor='white')

mscale = 7
gscale = 1400
filetag = "mscale_" + str(mscale) + "_spec.png"
wavelengths = np.array(wavelengths.reshape(2048))
wavelengths.shape
wavelengths.dtype
def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]
spect_indx, spect_value = find_nearest(wavelengths, 763.5)


def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]
flux_indx, flux_value = find_nearest(flux_lam, 763.5)

lam = flux_lam
y0 = flux[0] * gscale
y1 = flux[1] * mscale
y2 = flux[2] * 1.
flux_sum = (y0 + y1 + y2)
flux_sum = flux_sum / (1.0 * flux_sum[flux_indx])
plt.plot(flux_lam, flux_sum , 'm.', ms=3, label = 'summed')
plt.plot(wavelengths - 0.2, spect_mean/  spect_mean[spect_indx], label = 'ALEXIS')

plt.title('PEC spectrum: mscale =' +str(mscale) \
        + '; gscale = ' + str(gscale))
plt.xlabel("Wavelength (nm)", fontsize=16)
#pl.ylim(0,1.2)
plt.xlim(760, 900)
plt.legend(loc=2)
plt.savefig(filetag, dpi=300)
plt.show(fig4)
