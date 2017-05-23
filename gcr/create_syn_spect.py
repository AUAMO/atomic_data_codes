# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 12:56:16 2017

@author: nia00
"""
import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir('C:/Users/ivan/OneDrive/Research/ar_OES')
spect_dir = 'C:/Users/ivan/OneDrive/Research/ar_OES/DATA'
line_array=np.array([763.5, 794.8, 826.5, 852.1])
   
gscale = 2000
y0 = flux[0] * gscale
y1 = flux[1] * meta_scale
y2 = flux[2] * 1.
flux_sum = (y0 + y1 + y2)
def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]
flux_windx, flux_walue = find_nearest(flux_lam, \
                                     763.5)
spect_windx, spect_walue = find_nearest(wavelengths,\
                                       763.5)
fig1 = plt.figure(figsize=(8,6), facecolor="white")
for i in range(4):
    ylabel = str(line_array[i])
    def find_nearest(array, value):
        idx = (np.abs(array-value)).argmin()
        return idx, array[idx]
    flux_indx, flux_value = find_nearest(flux_lam, \
                                         line_array[i])   
    spect_indx, spect_value = find_nearest(wavelengths,\
                                           line_array[i])
    
    flux_sum = flux_sum / (flux_sum[flux_indx])
    
    plt.plot(wavelengths[spect_indx-6: spect_indx+10], \
             spect_mean[spect_indx-6:spect_indx+10] / \
             spect_mean[spect_windx], label = ylabel)
    plt.plot(flux_lam[flux_indx-6:flux_indx+10],\
              flux_sum[flux_indx-6:flux_indx+10] / \
              flux_sum[flux_windx], 'm.', ms=1)
plt.title("ALEXIS: Synthetic emission compared w/ experiment for strong lines.") 
plt.legend(loc=1)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Normalized Emission (arbitrary units)")
plt.xlim(750, 860)
plt.savefig("syn_spect.png", dpi=300)

plt.show(fig1)