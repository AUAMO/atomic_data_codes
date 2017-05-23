# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 10:48:59 2017

@author: nia00
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os

def ar_temp_finder(line_ratio, line):
    os.chdir('c:/Users/ivan/OneDrive/Research/ar_OES/DATA')

    def find_nearest(array, value):
        idx = (np.abs(array-value)).argmin()
        return idx, array[idx]

    # Import the .fits file which contains all the pec dat and
    # extract the information.
    hdulist = fits.open('pec_dat_650_950_lowdens.fits')
    pec_temps = hdulist[0].data       # The array of temps in eV
    n_meta = hdulist[1].data       # Array listing the metastable #
    pec_wave = hdulist[3].data       # Wavelengths corresponiding to each PEC
    pec_pec = hdulist[4].data.T     # 3-D array containing all PEC's

    # Declare some variables
    size_meta = n_meta.size          # Number of metastables
    n_pec = pec_wave.size        # Number of PEC's for each metastable

    # Create dictionaries to hold the invidual meta PEC's
    pec_meta = dict()
    for i in range(size_meta):
        pec_meta[i] = pec_pec[:, :, n_pec*i:n_pec*(i+1)]

    if line == 826:    
        upper, uval = find_nearest(pec_wave, 826.5)
        mscale =  3.58775946328 * line_ratio - 12.338615998
    if line == 852:
        upper, uval = find_nearest(pec_wave, 852.1)        
        mscale = 2.86026783707 * line_ratio - 7.82438767561

    lower, lval = find_nearest(pec_wave, 763.5)
    gscale = 2000

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

    return(pec_temps[indx], pec_x[indx2], mscale, pec_x, pec_array)


# %% Get info for 826 nm line

def get_all_temps_763_826():
    all_dat = np.loadtxt('C:/Users/ivan/OneDrive/Research/ar_OES/DATA/763_826/ALLPORTS.dat', delimiter='|')
    ratios = all_dat[:, 2]
    ratios.size
    temp_pec = list()
    temp_int = list()
    mfit = list()
    lrat = list()
    for i in range(np.size(ratios)):
        pec_temp, int_temp, mscale, pec_x, pec_array = \
                  ar_temp_finder(ratios[i], 826)
        temp_pec.append(pec_temp)
        temp_int.append(int_temp)
        mfit.append(mscale)
        lrat.append((mscale + 12.34)/3.59)

    temp_pec = np.array(temp_pec)
    temp_int = np.array(temp_int)
    mfit = np.array(mfit)
    temp_array = np.array([mfit, lrat, temp_pec, temp_int])
    return (temp_array, all_dat)

temp_array_826, all_dat = get_all_temps_763_826()
print(temp_array_826.T)


def square_fit(xx, aa, bb, cc):
    return (aa * xx**2 + bb * xx + cc)
popt, pcov = curve_fit(square_fit, temp_array_826[1], temp_array_826[3])

aa = popt[0]
bb = popt[1]
cc = popt[2]
xx = np.linspace(3.0, 8, 1000)
fit_826 = aa * xx ** 2 + bb * xx + cc
perr = np.sqrt(np.diag(pcov))





fig2 = plt.figure(figsize=(6,4), facecolor="white")
plt.plot(temp_array_826[1], temp_array_826[3], '*', color="orange")
plt.errorbar(temp_array_826[1], all_dat[:,0], yerr=(all_dat[:,0]*.25),fmt='o', fillstyle="none", color="blue", linestyle="None")
plt.ylim(3, 8.0)
plt.xlim(3, 8.0)

plt.show(fig2)

# %% Get info for 852 nm line

def get_all_temps_763_852():
    all_dat = np.loadtxt('C:/Users/ivan/OneDrive/Research/ar_OES/DATA/763_852/ALLPORTS.dat', \
                         delimiter='|')
    ratios = all_dat[:, 2]
    ratios.size
    temp_pec = list()
    temp_int = list()
    mfit = list()
    lrat = list()
    for i in range(np.size(ratios)):
        pec_temp, int_temp, mscale, pec_x, pec_array = \
                  ar_temp_finder(ratios[i], 852)
        temp_pec.append(pec_temp)
        temp_int.append(int_temp)
        mfit.append(mscale)
        lrat.append((mscale + 7.82)/2.86)

    temp_pec = np.array(temp_pec)
    temp_int = np.array(temp_int)
    mfit = np.array(mfit)
    temp_array = np.array([mfit, lrat, temp_pec, temp_int])
    return (temp_array, all_dat)

temp_array_852, all_dat = get_all_temps_763_852()
print(temp_array_852.T)


def square_fit(xx, aa, bb, cc):
    return (aa * xx**2 + bb * xx + cc)
popt, pcov = curve_fit(square_fit, temp_array_852[1], temp_array_852[3])

aa = popt[0]
bb = popt[1]
cc = popt[2]
xx = np.linspace(3.0, 8, 1000)
fit_852 = aa * xx ** 2 + bb * xx + cc
perr = np.sqrt(np.diag(pcov))

fit_string = r'T_e =  ' + '{:2.3f}'.format(aa) + r'x^2 ' + '{:2.3f}'.format(bb) + 'x +' + '{:2.3f}'.format(cc)
print(fit_string)
    #    print("n_e = {:6.2e} m^2 ".format(n_e))

fig2 = plt.figure(figsize=(8,6), facecolor="white")
plt.plot(temp_array_852[1], temp_array_852[3], '*', color="orange", ms=10, label="Predicted temps")
plt.plot(xx, fit_852, color="orange", label = "Quadratic fit")
plt.errorbar(temp_array_852[1], all_dat[:,0], yerr=(all_dat[:,0]*.25),fmt='o',\
             capsize=4, fillstyle="none", color="blue", label="Measured temps")
plt.legend(loc=1)
plt.text(4.7, 2.6, fit_string)
plt.title("ALEXIS: Predicted temperatures vs. measured line ratios.", weight="bold")
plt.xlabel("Line Ratio (763.5/852.1)")
plt.ylabel("Temperature (eV)")

plt.ylim(2.5, 8.5)
plt.xlim(4.5, 8.0)

plt.savefig("temp_finder.png", dpi=300)
plt.show(fig2)
# %%

fig1 = plt.figure(figsize=(6,4), facecolor="white")

plt.plot(temp_array_852[1], temp_array_852[3], '*', color="orange")
plt.plot(xx, fit_852, color="orange")

plt.ylim(4.5 , 6.0)
plt.xlim(4.5 , 8.0)
plt.show(fig1)