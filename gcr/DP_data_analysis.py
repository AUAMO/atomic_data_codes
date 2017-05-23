# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 10:06:27 2016

@author: ALEXIS
"""
import numpy as np
import matplotlib.pyplot as plt
import os

spect_dir = 'C:/Users/ivan/OneDrive/Research/ar_OES/DATA/BC1E'
os.chdir(spect_dir)

# %%
def retrieve_file_dir_data(spect_dir):
    from astropy.io import fits
    from scipy.optimize import curve_fit
    from scipy import stats
    print('filelist: ')
    filelist = list()
    spectlist = list()
    tempdenslist = list()
    temp_dens_av = list()
    yy_temps = list()
    yy_dens = list()
    
    for filename in os.listdir(spect_dir):
        if 'BC1G' in filename:
            continue
#        if 'CM100' in filename:
#            continue
#        if 'CM075' in filename:
#            continue
#        if 'CM050' in filename:
#            continue
        

        print(filename)
        hdulist = fits.open(filename)

        def spect_get():
            hdu_size = np.size(hdulist)
            scan_lambda = hdulist[hdu_size-2].data.T
            wavelengths = np.array(scan_lambda[0]).T
            scan_spect = hdulist[hdu_size-1].data.T

            spect = np.array(scan_spect[0])
            spec_av = np.mean(spect[:, 2000:2047])
            spect = spect - spec_av
            for i in range(1,10):
                tmp_data = spect
                spect = np.array(scan_spect[i])
                spec_av = np.mean(spect[:, 2000:2047])
                spect = spect - spec_av
                spect_array = np.vstack((tmp_data, spect))
                spect = spect_array

            spect_mean = spect.mean(axis=0)

            return (wavelengths, spect_array, spect_mean)

        def temp_dens_get(indx):
            hdu_size = np.size(hdulist)
            scan_size = hdu_size - 2
            scan_dat = {}
            for i in range(1,scan_size):
                scan_dat["scan{0}".format(i)] = hdulist[i].data

            DP_START = scan_dat['scan5'][indx][0]
            DP_STOP  = np.abs(DP_START)
            DP_STEP  = scan_dat['scan5'][indx][1]
            DP_RAW   = scan_dat['scan5'][indx][2]
            V_TRACE  = np.arange(DP_START, DP_STOP + DP_STEP, DP_STEP)
            isat_x = V_TRACE[0:20] ; isat_y = DP_RAW[0:20]
            isat_slope, isat_intercept, isat_rval, isat_pval, std_err = \
                stats.linregress(isat_x, isat_y)
            isat_trace = isat_slope * V_TRACE + isat_intercept
            zero_x = V_TRACE[75:85] ; zero_y = isat_trace[75:85]
            zero_slope, zero_intercept, zero_rval, zero_pval, std_err = \
                stats.linregress(zero_x, zero_y)
            def tanh_fit(xx, a0, a1, a2, a3, a4):
                return(a0 * (xx - a1) + a2 * np.tanh((xx - a1)/(a3)) + a4)
            ifit_params = curve_fit(tanh_fit, V_TRACE, DP_RAW)

            # Area of probe tip with 1 mm diameter and 2.5 mm length
            tipArea = 0.5 * (np.pi * 0.5e-3 ** 2) + (np.pi * 1e-3 *  2.5e-3)
            a0 = float(ifit_params[0][0])  # ion saturation slope
            a1 = float(ifit_params[0][1])  # voltage offset / V_f
            a2 = float(ifit_params[0][2])  # ion saturation current
            a3 = float(ifit_params[0][3])  # slope at zstatsero volts
            a4 = float(ifit_params[0][4])  # current offset
            ifit_trace = a0*(V_TRACE-a1) + a2*np.tanh((V_TRACE-a1)/(a3)) + a4
            zero_x = V_TRACE[75:85] ; zero_y = ifit_trace[75:85]
            zero_slope, zero_intercept, zero_rval, zero_pval, std_err = \
                stats.linregress(zero_x, zero_y)
            di_dv = a0 + a2 / (a3 * np.cosh(a1/a3)**2)
            t_e = np.abs(a2)/ (2. * di_dv)
        #    print("t_e = {:6.2e} eV ".format(t_e))
            n_e = (a2 / (1.6e-19 * tipArea)) * (6.6e-26 / \
                  (1.602e-19 * t_e))**0.5
        #    print("n_e = {:6.2e} m^2 ".format(n_e))
        #    print("isat = {:6.2e} Amps".format(a2))
        #    print("slope at V_0 = {:6.2e} Amps/Volt".format(di_dv))
            radial_position = indx * 2
        #    print("radial position = " , radial_position, "mm.")
            temp_dens_data = [[radial_position],[t_e],[n_e],[a2],[di_dv]]
            temp_dens_array = np.array(temp_dens_data)

            return temp_dens_array


        data = temp_dens_get(0)
        for i in range(1,26):
            tmp_data = data
            data = temp_dens_get(i)
            temp_dens_data = np.hstack((tmp_data, data))
            data = temp_dens_data
        wavelengths, spect_array, spect_mean = spect_get()
        tmp_av_data = temp_dens_data[0:39].mean(axis=1)
        temp_dens_av.append(tmp_av_data)

        filelist.append(filename)
        spectlist.append(spect_mean)
        tempdenslist.append(temp_dens_data)
        yy_temps.append(temp_dens_data[1])
        yy_dens.append(temp_dens_data[2])
    
    xx_rad = temp_dens_data[0]

    
    return (xx_rad, yy_temps, yy_dens, filelist)
            
xx, yy1, yy2, filelist = retrieve_file_dir_data(spect_dir)

# %%



for i in range(24):
#    tag = filelist[i][]
    plt.xlim(0,40)
    plt.ylim(3,8)
    plt.plot(xx, yy1[i])

# %%

for i in range(24):
#    tag = filelist[i][]
    plt.xlim(0,40)
    plt.ylim(0e16, 3.5e16)
    plt.plot(xx, yy2[i])







# %% Plot the results.
