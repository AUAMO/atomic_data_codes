# -*- coding: utf-8 -*-
# %%
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
import datetime; now = datetime.datetime.now()
import tkinter
from tkinter import filedialog

 # we don't want a full GUI, so keep the root window from appearing
root = tkinter.Tk()
root.withdraw()

# show an "Open" dialog box and return the path to the selected file
filename = filedialog.askopenfilename()

hdulist = fits.open(filename)
hdu_info = hdulist.info()
print(hdu_info)

# %%
scan_params = hdulist[0].data
hdu_size = np.size(hdulist)
scan_size = hdu_size - 2

scan_dat = {}
for i in range(1,scan_size):
    scan_dat["scan{0}".format(i)] = hdulist[i].data

scan_lambda = hdulist[hdu_size-2].data.T
wavelengths = np.array(scan_lambda[0]).T
scan_spect = hdulist[hdu_size-1].data.T
spect = np.array(scan_spect[2]).T


DP_START = scan_dat['scan1'][5][0]
DP_STOP  = np.abs(DP_START)
DP_STEP  = scan_dat['scan1'][5][1]
DP_RAW   = scan_dat['scan1'][5][2]
V_TRACE  = np.arange(DP_START, DP_STOP + DP_STEP, DP_STEP)

# Area of probe tip with 1 mm diameter and 2.5 mm length
tipArea = (np.pi * 0.5e-3 ** 2) + (np.pi * 1e-3 *  2.5e-3)
print('cell 1')
# %%
isat_x = V_TRACE[0:20] ; isat_y = DP_RAW[0:20]

isat_slope, isat_intercept, isat_rval, isat_pval, std_err = \
    stats.linregress(isat_x, isat_y)
isat_trace = isat_slope * V_TRACE + isat_intercept

zero_x = V_TRACE[75:85] ; zero_y = isat_trace[75:85]
zero_slope, zero_intercept, zero_rval, zero_pval, std_err = \
    stats.linregress(zero_x, zero_y)
zero_trace = zero_slope * V_TRACE[55:120] + zero_intercept

def tanh_fit(xx, a0, a1, a2, a3, a4):
    return(a0 * (xx - a1) + a2 * np.tanh((xx - a1)/(a3)) + a4)
ifit_params = curve_fit(tanh_fit, V_TRACE, DP_RAW)
a0 = float(ifit_params[0][0])  # ion saturation slope
a1 = float(ifit_params[0][1])  # voltage offset / V_f
a2 = float(ifit_params[0][2])  # ion saturation current
a3 = float(ifit_params[0][3])  # slope at zero volts
a4 = float(ifit_params[0][4])  # current offset
ifit_trace = a0*(V_TRACE-a1) + a2*np.tanh((V_TRACE-a1)/(a3)) + a4

zero_x = V_TRACE[75:85] ; zero_y = ifit_trace[75:85]
zero_slope, zero_intercept, zero_rval, zero_pval, std_err = \
    stats.linregress(zero_x, zero_y)
esatx = V_TRACE[145:155] ; esaty = ifit_trace[145:155]
esat_slope, esat_intercept, esat_rval, esat_pval, esat_err = \
    stats.linregress(esatx, esaty)
    
zero_trace = zero_slope * V_TRACE[55:120] + zero_intercept
esat_trace = esat_slope * V_TRACE[90:140] + esat_intercept

di_dv = a0 + a2 / (a3 * np.cosh(a1/a3)**2)

# %% Plot the results.
fig1 = plt.figure(figsize=(8,6), facecolor="white")
ax = fig1.add_subplot(111)
plt.grid()
plt.plot(V_TRACE, ifit_trace / 1e-3, 'r', lw=2)
plt.plot(V_TRACE, DP_RAW / 1e-3, 'b.', ms=3)
plt.plot(V_TRACE[0:81], isat_trace[0:81] / 1e-3, 'g-.', ms=1, label=r'$I_{si}$')
plt.plot(V_TRACE[55:120], zero_trace / 1e-3, 'm-.', ms=1, label = 'transition')
plt.plot(V_TRACE[90:140], esat_trace / 1e-3, 'o-.', ms=1, label=r'$I_{se}$')

plt.axhline(0.019, color="black", ls="dashdot")
ax.annotate(r'$V_p$', xy=(10, 0.025), xytext=(-40, 0.021), fontsize=16)

plt.axhline(0.0, color="black", ls="dashdot")
ax.annotate(r'$V_f$', xy=(0,0), xytext=(-40, 0.001), fontsize=16)


print()
print("isat = {:6.2e} Amps".format(a2))
print("slope at V_0 = {:6.2e} Amps/Volt".format(di_dv))

t_e = np.abs(a2)/ (2. * di_dv)
print("t_e = {:6.2e} eV ".format(t_e))

n_e = (a2 / (1.6e-19 * tipArea)) * (6.6e-26 / (1.602e-19 * t_e))**0.5
print("n_e = {:6.2e} m^2 ".format(n_e))
plt.title("Characteristic I-V trace (ALEXIS)", fontsize=16, weight="bold")
plt.legend(loc=4, prop={'size':14})
plt.xlabel("Bias Voltage", fontsize=16)
plt.ylabel("Current (mA)", fontsize=16)
plt.savefig("IV_trace.png", dpi=300)
plt.show(fig1)
# %% Plot the results.
plt.clf()
plt.grid()
plt.plot(wavelengths, spect - 2500)

plt.xlim(665, 860)
plt.ylim(0, 52000)
plt.show()
