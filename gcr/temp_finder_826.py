# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 10:48:59 2017

@author: nia00
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os

def ar_temp_finder(line_ratio):

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

    # Create dictionaries to hold the invidual meta PEC's and Broadened flux
    pec_meta = dict()
    for i in range(size_meta):
        pec_meta[i] = pec_pec[:, :, n_pec*i:n_pec*(i+1)]

    gscale = 2000
    mscale =  3.58775946328 * line_ratio - 12.338615998
#    mscale = 0.249341860455 * line_ratio + 3.64323291233

    lower, lval = find_nearest(pec_wave, 763.5)
    upper, uval = find_nearest(pec_wave, 826.5)

    pec_ratio_weighted = (gscale * pec_meta[0][:, :, lower] + mscale \
        * pec_meta[1][:, :, lower] + pec_meta[2][:, :, lower]) / (gscale * \
        pec_meta[0][:, :, upper] + mscale * pec_meta[1][:, :, upper] \
        + pec_meta[2][:, :, upper])

    pec_array = np.mean(pec_ratio_weighted, axis=1)
    from scipy.interpolate import interp1d
    pec_fit = interp1d(pec_temps, pec_array, kind = 'quadratic')

    indx, val = find_nearest(pec_array, line_ratio)

    pec_x = np.linspace(0.5,20,1000)
    indx2, val2 = find_nearest(pec_fit(pec_x), line_ratio)

    indx, val = find_nearest(pec_array, line_ratio)
    pec_x = np.linspace(0.5,20,1000)
    indx2, val2 = find_nearest(pec_fit(pec_x), line_ratio)

    return(pec_temps[indx], pec_x[indx2], mscale, pec_x, pec_array)

def get_all_temps_763_826():
    os.chdir('C:/Users/ivan/OneDrive/Research/ar_OES/DATA/763_826')
    all_dat = np.loadtxt('ALLPORTS.dat', delimiter='|')
    ratios = all_dat[:, 2]
    ratios.size
    temp_pec = list()
    temp_int = list()
    mfit = list()
    lrat = list()
    for i in range(np.size(ratios)):
        pec_temp, int_temp, mscale, pec_x, pec_array = ar_temp_finder(ratios[i])
        temp_pec.append(pec_temp)
        temp_int.append(int_temp)
        mfit.append(mscale)
        lrat.append((mscale + 12.34)/3.59)
    temp_pec = np.array(temp_pec)
    temp_int = np.array(temp_int)
    mfit = np.array(mfit)
    temp_array = np.array([mfit, lrat, temp_pec, temp_int])
    return (temp_array, all_dat)


temp_array, all_dat = get_all_temps_763_826()
print(temp_array.T)

fig1 = plt.figure(figsize=(6,4), facecolor="white")
plt.plot(temp_array[1], temp_array[3], '*', color="orange")
#plt.ylim(4.5, 6.0)
plt.show(fig1)

fig2 = plt.figure(figsize=(6,4), facecolor="white")

plt.plot(temp_array[1], temp_array[3], '*', color="orange")
#plt.plot(temp_array[1], all_dat[:,0], 'bx', ms=3)

plt.errorbar(temp_array[1], all_dat[:,0], yerr=(all_dat[:,0]*.333),fmt='o', fillstyle="none", color="blue", linestyle="None")
#plt.ylim(3, 8.0)

plt.show(fig2)


####################################################################