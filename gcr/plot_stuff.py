# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 23:21:41 2017

@author: ivan
"""
import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir('C:/Users/ivan/GoogleDrive/Research/ar_OES/DATA/763_852')

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
data2 = np.loadtxt('pec_spect.dat')
data  = np.loadtxt('temp_5_dens2.5e16.dat')
indx2 = find_nearest(data2[0], 763.521)
indx = find_nearest(data[0], 763.521)

print(data[:,indx])
print(data2[:,indx2])

fig2 = plt.figure(figsize=(8, 6))
plt.plot(data[0], data[1] / data[1, indx])
plt.plot(data2[0], data2[1] / data2[1, indx2], 'g.', ms=3)
plt.xlabel('Wavelength (nm)', fontsize=14)
plt.ylabel('Line Intensity (Arbitrary Units)', fontsize=14)
plt.title(r'Synthetic Spectrum vs ALEXIS: t_e = 5 eV, n_e = 2.5e16 m^{-3}', fontsize=12 )
plt.xlim(760, 860)
plt.ylim(0, 2.5)
plt.savefig('syn_plot.png', dpi=300)
plt.show()