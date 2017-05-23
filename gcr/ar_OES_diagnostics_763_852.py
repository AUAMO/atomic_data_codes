# %%  Open an appropriate ALEXIS data file and import temp, dens and spectral data.
import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir('C:/Users/ivan/Dropbox/Research/ar_OES')
spect_dir = 'C:/Users/ivan/Dropbox/Research/ar_OES/DATA'

# %%


def retrieve_ALEXIS_data(spect_dir):
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

    def temp_dens_get(indx):
        hdu_size = np.size(hdulist)
        scan_size = hdu_size - 2
        scan_dat = {}
        for i in range(1, scan_size):
            scan_dat["scan{0}".format(i)] = hdulist[i].data

        DP_START = scan_dat['scan5'][indx][0]
        DP_STOP = np.abs(DP_START)
        DP_STEP = scan_dat['scan5'][indx][1]
        DP_RAW = scan_dat['scan5'][indx][2]
        V_TRACE = np.arange(DP_START, DP_STOP + DP_STEP, DP_STEP)
        isat_x = V_TRACE[0:20]
        isat_y = DP_RAW[0:20]
        isat_slope, isat_intercept, isat_rval, isat_pval, std_err = \
            stats.linregress(isat_x, isat_y)
        isat_trace = isat_slope * V_TRACE + isat_intercept
        zero_x = V_TRACE[75:85]
        zero_y = isat_trace[75:85]
        zero_slope, zero_intercept, zero_rval, zero_pval, std_err = \
            stats.linregress(zero_x, zero_y)

        def tanh_fit(xx, a0, a1, a2, a3, a4):
            return(a0 * (xx - a1) + a2 * np.tanh((xx - a1) / (a3)) + a4)
        ifit_params = curve_fit(tanh_fit, V_TRACE, DP_RAW)

        # Area of probe tip with 1 mm diameter and 2.5 mm length
        tipArea = 0.5 * (np.pi * 0.5e-3 ** 2) + (np.pi * 1e-3 * 2.5e-3)
        a0 = float(ifit_params[0][0])  # ion saturation slope
        a1 = float(ifit_params[0][1])  # voltage offset / V_f
        a2 = float(ifit_params[0][2])  # ion saturation current
        a3 = float(ifit_params[0][3])  # slope at zstatsero volts
        a4 = float(ifit_params[0][4])  # current offset
        ifit_trace = a0 * (V_TRACE - a1) + a2 * np.tanh((V_TRACE - a1) / (a3)) + a4
        zero_x = V_TRACE[75:85]
        zero_y = ifit_trace[75:85]
        zero_slope, zero_intercept, zero_rval, zero_pval, std_err = \
            stats.linregress(zero_x, zero_y)
        di_dv = a0 + a2 / (a3 * np.cosh(a1 / a3)**2)
        t_e = np.abs(a2) / (2. * di_dv)
    #    print("t_e = {:6.2e} eV ".format(t_e))
        n_e = (a2 / (1.6e-19 * tipArea)) * (6.6e-26 /
                                            (1.602e-19 * t_e))**0.5
    #    print("n_e = {:6.2e} m^2 ".format(n_e))
    #    print("isat = {:6.2e} Amps".format(a2))
    #    print("slope at V_0 = {:6.2e} Amps/Volt".format(di_dv))
        radial_position = indx * 2
    #    print("radial position = " , radial_position, "mm.")
        temp_dens_data = [[radial_position], [t_e], [n_e], [a2], [di_dv]]
        temp_dens_array = np.array(temp_dens_data)

        return temp_dens_array

    data = temp_dens_get(0)
    for i in range(1, 26):
        tmp_data = data
        data = temp_dens_get(i)
        temp_dens_data = np.hstack((tmp_data, data))
        data = temp_dens_data
    temp_dens_av = temp_dens_data.mean(axis=1)

    return (filename, temp_dens_data, temp_dens_av, wavelengths, spect_array, spect_av)


filename, temp_dens_data, temp_dens_av, wavelengths, spect_array, spect_mean = retrieve_ALEXIS_data(
    spect_dir)
temp_av = np.sqrt(np.mean(np.square(temp_dens_data[1, 0:15])))
dens_av = np.sqrt(np.mean(np.square(temp_dens_data[2, 0:15])))
temp_min = np.max(temp_dens_data[1, 0:15])
temp_max = np.min(temp_dens_data[1, 0:15])

# Find the ratio of the 763 / 852 nm lines.

header_array = np.array([['radius'], ['t_e'], ['n_e'], ['isat'], ['di_dv0']])


def get_spect_ratios(wavelengths, spect):
    import peakutils

    def find_nearest(array, value):
        idx = (np.abs(array - value)).argmin()
        return idx
    values = np.array([763.5, 852.1, 919.0])

    tmp_ratio = np.zeros(3)
    tmp_value = np.zeros(4).astype(int)

    # temp_spect = spect[877:1614]    ## Uncomment these lines for  negligible speed increase
    # temp_wavelengths = wavelengths[877:1614]  ## Adjust slice to focus on
    # region of interest in spectrum.
    peaks = peakutils.indexes(spect, thres=0.05 / max(spect), min_dist=5)
    temp_spect = spect[peaks]
    temp_spect = np.reshape(temp_spect, (np.size(peaks), 1))
    temp_wavelengths = wavelengths[peaks]
    for j in range(np.size(values)):
        tmp_value[j] = find_nearest(temp_wavelengths, values[j])
    for k in range(np.size(tmp_ratio)):
        tmp_ratio[k] = temp_spect[tmp_value[k]] / temp_spect[tmp_value[k + 1]]
    ratio_763_852 = tmp_ratio[0] * 1.05305396186 / 1.04566612136
    ratio_852_919 = tmp_ratio[1] * 1.04566612136 / 1.04267580497

    return (ratio_763_852, ratio_852_919)


ratio_763_852, ratio_852_919 = get_spect_ratios(wavelengths, spect_mean)
temp_ratio = ratio_763_852
ratio = ratio_763_852
dens_ratio = ratio_852_919
limit = 16
temp_av = np.sqrt(np.mean(np.square(temp_dens_data[1, 0:limit])))
dens_av = np.sqrt(np.mean(np.square(temp_dens_data[2, 0:limit])))
temp_min = np.max(temp_dens_data[1, 0:limit])
temp_max = np.min(temp_dens_data[1, 0:limit])
# print(ratio)

# Find the metastable fraction that best fits the expermental line ratio, average temp.
#    This function loops over three hundred possible meatstable fractions until it finds the
#    one with a line ratio and temperature that best fit the experimental data.


def get_meta_scale(temp, ratio, lval, uval):
    from astropy.io import fits
    from scipy.optimize import curve_fit

    def find_nearest(array, value):
        idx = (np.abs(array - value)).argmin()
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
    
# Spectrometer Resolution (nm)
    lambda_d = 1.7
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

    t_indx, t_val = find_nearest(pec_temps, temp_av)
    d_indx = 19
    
    # Doppler Broaden the array for the specified t_indx and d_indx
    for i in range(n_wavel):
        for j in range(pec_wave.size):
            for k in range(size_meta):
                if flux_lam[i] > (pec_wave[j]-fwhm*5.0) \
                and flux_lam[i] < (pec_wave[j]+fwhm*5.0):
                    flux[k][i] = flux[k][i] + pec_meta[k][t_indx, \
                    d_indx, j] * (2 * np.pi) ** (-0.5)/lambda_d * \
                    np.exp(-abs((flux_lam[i] - pec_wave[j]) / \
                    lambda_d) ** 2.0)

                    
                    
    gscale = 2000
    mscale = np.linspace(0.1, 30, 300)

    # pec_av = np.zeros(300)
    lower, lval = find_nearest(pec_wave, lval)
    upper, uval = find_nearest(pec_wave, uval)

    pec_meta[0] = np.mean(pec_meta[0], axis=1)
    pec_meta[1] = np.mean(pec_meta[1], axis=1)
    pec_meta[2] = np.mean(pec_meta[2], axis=1)

    pec_scaled_array = list()
    for i in range(np.size(mscale)):
        pec_ratio_weighted = (gscale * pec_meta[0][:, lower] + mscale[i]
                              * pec_meta[1][:, lower] + pec_meta[2][:, lower]) / (gscale *
                                                                                  pec_meta[0][:, upper] +
                                                                                  mscale[i] *
                                                                                  pec_meta[1][:,
                                                                                              upper]
                                                                                  + pec_meta[2][:, upper])
        pec_scaled_array.append(pec_ratio_weighted)
    pec_scaled_array = np.array(pec_scaled_array)

    indx, val = find_nearest(pec_temps, temp)

    if (temp - val) > 0.0:
        pec_indx_low = indx
        pec_indx_high = indx + 1
        pec_temp_low = pec_temps[indx]
        pec_temp_high = pec_temps[indx + 1]
        # print("tag")
    else:
        pec_indx_low = indx - 1
        pec_indx_high = indx
        pec_temp_low = pec_temps[indx - 1]
        pec_temp_high = pec_temps[indx]
        # print("you're it")

    pec_x = np.linspace(0.5, 20, 1000)
    from scipy.interpolate import interp1d

    min_vals = list()
    r_list = list()
    for i in range(np.size(mscale)):
        pec_fit = interp1d(pec_temps, pec_scaled_array[i], kind='cubic')
        fit_indx, fit_val = find_nearest(pec_fit(pec_x), ratio)
        if (pec_x[fit_indx] > pec_temp_low) and (pec_x[fit_indx] < pec_temp_high):
            min_vals.append([i, fit_indx, np.abs(fit_val - ratio)])
    min_vals = np.array(min_vals)

    indx, arg_min = find_nearest(min_vals[:, 2], 0.)
    m_indx = int(min_vals[indx, 0])
    # print(min_vals)
    # print(indx, m_indx, mscale[m_indx])

    return (m_indx, mscale[m_indx], pec_temps, pec_scaled_array[m_indx], flux_lam, flux)


meta_indx, meta_scale, pec_temps, pec_array, flux_lam, flux = get_meta_scale(temp_av, ratio, 763.5, 852.1)
filetag = filename[52:116]
print(filename[52:116] + "|{:6.2f}".format(temp_av) + "|{:6.2e}".format(dens_av) +
      "| {:6.2f}".format(ratio) + "|{:5.1f}".format(meta_scale))

print(meta_scale)
###################################    PLOTS  #############################################
# %%

pec_x = np.linspace(0.5, 20, 1000)
from scipy.interpolate import interp1d
terp = interp1d(pec_temps, pec_array, kind='cubic')
pec_fit = terp(pec_x)

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]
indx, av_temp = find_nearest(pec_x, temp_av)
theor_ratio = pec_fit[indx]

# Plot temp, dens profiles, and the line ratio plot of best metastable ratio fit.
def plot_ALEXIS_data():
    os.chdir('C:/Users/ivan/Dropbox/Research/ar_OES/DATA/763_852/meta_scale_images')
    fig = plt.subplots(2, 2, figsize=(9, 6))
    plt.subplot(221)
    label = "Average Temp = {:5.2f} (eV)".format(temp_av)
    plt.plot(temp_dens_data[0, :], temp_dens_data[1, :], label=label)
    plt.ylim(2,8)
    plt.xlim(0,40)
    plt.legend(loc=3)
    plt.subplot(222)
    label = "Average Dens = {:5.2e} m^-3".format(dens_av)
    plt.plot(temp_dens_data[0, :], temp_dens_data[2, :], label=label)
    plt.ylim(5e15,2e16)
    plt.xlim(0,40)
    plt.legend(loc=3)
    plt.subplot(212)
    m_label='mscale = ' + str(meta_scale)
    plt.text(0.1,7, m_label, size=10)
    plt.plot(pec_temps, pec_array, 'rx', label='weighted PEC ratio')
    plt.axvline(temp_av, color="orange", label="exp temp")
    plt.axhline(ratio, color='red', label="exp ratio")
    plt.axhline(theor_ratio, color='k', ls='--', label="theor ratio")
    
    plt.legend(loc=4)
    plt.suptitle('Temperature, Density and Metastable Ratio: ALEXIS', weight = 'bold')
#    plt.savefig(filetag + ".png", dpi=300)
    plt.savefig("theor_mfrac.png", dpi=300)
    plt.show(fig)


plot_ALEXIS_data()

# %%


                    
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
    plt.plot(flux_lam[flux_indx-20:flux_indx+20],\
              flux_sum[flux_indx-20:flux_indx+20] / \
              flux_sum[flux_windx], 'm.', ms=2)
plt.title("ALEXIS: Synthetic emission compared w/ experiment for strong lines.") 
plt.legend(loc=1)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Normalized Emission (arbitrary units)")
plt.xlim(750, 860)
plt.savefig("syn_spect.png", dpi=300)

plt.show(fig1)


# %%


def plot_vport_E1():
    from scipy import stats
    os.chdir("c:/Users/ivan/Dropbox/Research/ar_OES/DATA/763_852")
    BC1E1_dat = np.loadtxt('BC1E_1.dat', delimiter='|')
    slope, intercept, r_value, p_value, std_err = stats.linregress(BC1E1_dat[:, 2], BC1E1_dat[:, 3])
    lin_fit = BC1E1_dat[:, 2] * slope + intercept
    r_squared = r_value ** 2

    meta_av_E = np.mean(BC1E1_dat[:, 3])
    label = ("Meta_Av = {:6.2f} ".format(meta_av_E))

    fig1 = plt.figure(figsize=(9, 6))
    plt.title("BC1E_1")
    plt.plot(BC1E1_dat[:, 2], BC1E1_dat[:, 3], 'r*', label=label)
    plt.plot(BC1E1_dat[:, 2], lin_fit, 'r', label=("R^2= {:6.2f}".format(r_squared)))
    
    plt.legend(loc=2)
    plt.savefig('vport_E1.png', dpi=300)
    plt.show(fig1)


plot_vport_E1()


# %%


def plot_vport_E2():
    from scipy import stats
    os.chdir("C:/Users/ivan/Dropbox/Research/ar_OES/DATA/763_852")
    BC1E2_dat = np.loadtxt('BC1E_2.dat', delimiter='|')
    slope, intercept, r_value, p_value, std_err = stats.linregress(BC1E2_dat[:, 2], BC1E2_dat[:, 3])
    lin_fit = BC1E2_dat[:, 2] * slope + intercept
    r_squared = r_value ** 2
    meta_av_E = np.mean(BC1E2_dat[:, 3])
    label = ("Meta_Av = {:6.2f} ".format(meta_av_E))
    fig2 = plt.figure(figsize=(9, 6))
    plt.title("BC1E_2")
    plt.plot(BC1E2_dat[:, 2], BC1E2_dat[:, 3], 'r*', label=label)
    plt.plot(BC1E2_dat[:, 2], lin_fit, 'r', label=("R^2= {:6.2f}".format(r_squared)))
    plt.legend(loc=2)
    plt.savefig('vport_E2.png', dpi=300)
    plt.show(fig2)


plot_vport_E2()

# %%


def plot_vport_C():
    from scipy import stats
    os.chdir("C:/Users/ivan/Dropbox/Research/ar_OES/DATA/763_852")
    BC1C_dat = np.loadtxt('BC1C.dat', delimiter='|')
    slope, intercept, r_value, p_value, std_err = stats.linregress(BC1C_dat[:, 2], BC1C_dat[:, 3])
    lin_fit = BC1C_dat[:, 2] * slope + intercept
    r_squared = r_value ** 2
    meta_av_E = np.mean(BC1C_dat[:, 3])

    fig3 = plt.figure(figsize=(9, 6))
    plt.title("BC1C")
    plt.plot(BC1C_dat[:, 2], BC1C_dat[:, 3], 'r*', label="Viewport C")
    plt.plot(BC1C_dat[:, 2], lin_fit, 'r', label=("R^2= {:6.2f}".format(r_squared)))
    plt.legend(loc=2)
    plt.savefig('vport_C.png', dpi=300)
    plt.show(fig3)


plot_vport_C()

# %%


def plot_vport_A():
    from scipy import stats
    os.chdir("C:/Users/ivan/Dropbox/Research/ar_OES/DATA/763_852")
    BC1A_dat = np.loadtxt('BC1A.dat', delimiter='|')
    slope, intercept, r_value, p_value, std_err = stats.linregress(BC1A_dat[:, 2], BC1A_dat[:, 3])
    lin_fit = BC1A_dat[:, 2] * slope + intercept
    r_squared = r_value ** 2
    meta_av_E = np.mean(BC1A_dat[:, 3])
    label = ("Meta_Av = {:6.2f} ".format(meta_av_E))

    fig3 = plt.figure(figsize=(9, 6))
    plt.title("BC1A")
    plt.plot(BC1A_dat[:, 2], BC1A_dat[:, 3], 'r*', label=label)
    plt.plot(BC1A_dat[:, 2], lin_fit, 'r', label=("R^2= {:6.2f}".format(r_squared)))
    plt.legend(loc=2)
    plt.savefig('vport_A.png', dpi=300)
    plt.show(fig3)


plot_vport_A()

# %%


def plot_allports():
    from scipy import stats
    os.chdir("C:/Users/ivan/Dropbox/Research/ar_OES/DATA/763_852")
    all_dat = np.loadtxt('ALLPORTS.dat', delimiter='|')
    slope, intercept, r_value, p_value, std_err = stats.linregress(all_dat[:, 2], all_dat[:, 3])
    lin_fit = all_dat[:, 2] * slope + intercept
    print("slope, intercept = ", slope, intercept)
    r_squared = r_value ** 2
    meta_av_E = np.mean(all_dat[:, 2])

    fig3 = plt.figure(figsize=(9, 6))
    plt.title("All Viewports")
    plt.plot(all_dat[:, 2], all_dat[:, 3], 'r*', label="All Viewports")
    plt.plot(all_dat[:, 2], lin_fit, 'r', label=("R^2= {:6.2f}".format(r_squared)))
    plt.xlabel("Experimental line ratio")
    plt.ylabel("Theoretical metastable fraction (1s5/1s3)")
    plt.title("ALEXIS: Calculated 1s5/1s3 metastable fractions using 852.1 nm line", weight="bold")
    plt.legend(loc=2)
    plt.savefig('vport_ALL_852.png', dpi=300)
    plt.show(fig3)


plot_allports()

##########################################################################
# %% Use  experimental line ratio and metastable fraction from the linear fit from allports plot to
#    calculate temperature.


def ar_temp_finder(line_ratio):
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    os.chdir('C:/Users/ivan/Dropbox/Research/ar_OES/DATA')

    def find_nearest(array, value):
        idx = (np.abs(array - value)).argmin()
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
        pec_meta[i] = pec_pec[:, :, n_pec * i:n_pec * (i + 1)]

    gscale = 2000
    # slope, intercept =  2.86026783707 -7.82438767561
    mscale = 2.86026783707 * line_ratio - 7.82438767561
    # mscale = 0.249341860455 * line_ratio + 3.64323291233

    lower, lval = find_nearest(pec_wave, 763.5)
    upper, uval = find_nearest(pec_wave, 852.5)

    pec_ratio_weighted = (gscale * pec_meta[0][:, :, lower] + mscale \
            * pec_meta[1][:, :, lower] + pec_meta[2][:, :, lower])  \
            / (gscale * pec_meta[0][:, :, upper] + mscale * pec_meta[1][:, \
             :, upper] + pec_meta[2][:, :, upper])

    pec_array = np.mean(pec_ratio_weighted, axis=1)
    from scipy.interpolate import interp1d
    pec_fit = interp1d(pec_temps, pec_array, kind='cubic')

    indx, val = find_nearest(pec_array, line_ratio)

    pec_x = np.linspace(0.5, 20, 1000)
    indx2, val2 = find_nearest(pec_fit(pec_x), line_ratio)

    indx, val = find_nearest(pec_array, line_ratio)
    pec_x = np.linspace(0.5, 20, 1000)
    indx2, val2 = find_nearest(pec_fit(pec_x), line_ratio)

    # Uncomment to generate plots ==============================================
    # fig2 = plt.figure(figsize=(9,4), facecolor='white')
    # plt.title('Temperature Ratio: mscale = {:6.2f}'.format(mscale))
    # plt.xlabel('Temperature (eV)')
    # plt.plot(pec_temps, pec_array, 'r*')
    # plt.plot(pec_x, pec_fit(pec_x), 'r')
    # plt.axhline(line_ratio)
    # plt.axvline(pec_x[indx2], label = 'Temp = {:6.2f}'.format(pec_x[indx2]))
    # plt.legend(loc = 4)
    # plt.xlim(0, 20)
    # plt.ylim(0,10)
    # plt.savefig("temp_finder_" + str(line_ratio) + "_763_826.png", dpi=300)
    # plt.show(fig2)
    # ==========================================================================
    # print("ratio, temp, m_frac  = {:6.2f}".format(line_ratio), " ; {:6.2f}".format(pec_x[indx2]), \
    #        " ; {:6.2f}".format(mscale) )

    return(pec_temps[indx], pec_x[indx2], mscale, pec_x, pec_array)

# %% Calls the previous function in a loop over all measured line ratios. Stores output in temp_array.


def get_all_temps_763_852():
    import numpy as np
    import os
    os.chdir('C:/Users/ivan/Dropbox/Research/ar_OES/DATA/763_852')
    all_dat = np.loadtxt('ALLPORTS.dat', delimiter='|')
    ratios = all_dat[:, 2]
    ratios.size
    temp_pec = list()
    temp_int = list()
    mfit = list()
    for i in range(np.size(ratios)):
        pec_temp, int_temp, mscale, pec_x, pec_array = ar_temp_finder(ratios[i])
        temp_pec.append(pec_temp)
        temp_int.append(int_temp)
        mfit.append(mscale)
    temp_pec = np.array(temp_pec)
    temp_int = np.array(temp_int)
    mfit = np.array(mfit)
    temp_array = np.array([temp_pec, temp_int, mfit])
    return (temp_array, temp_int, mfit)


temp_array, temp_int, mfit = get_all_temps_763_852()
print(temp_array.T)
##########################################################################
##########################################################################
#    How well do these results compare with the results from the 763/826 line ratio?
# %% Plot metastable fractions 763/826 ratio vs 763/852 ratio.


def meta_comp_826_852():
    os.chdir('C:/Users/ivan/Dropbox/Research/ar_OES/DATA/763_852')
    from scipy import stats
    m_comp = np.loadtxt('m_compare.dat', delimiter='|')
    slope, intercept, r_value, p_value, std_err = stats.linregress(m_comp[:, 0], m_comp[:, 1])
    lin_fit = m_comp[:, 0] * slope + intercept
    r_squared = r_value ** 2
    fig1 = plt.figure(figsize=(9, 6))
    plt.title("m_compare")
    plt.xlabel("m_frc (763/826)")
    plt.ylabel("m_frc (763/852)")
    plt.plot(m_comp[:, 0], m_comp[:, 1], 'r*')
    plt.plot(m_comp[:, 0], lin_fit, 'r', label=("R^2= {:6.2f}".format(r_squared)))
    plt.legend(loc=2)
    plt.savefig('temp_compare.png', dpi=300)
    plt.show(fig1)


meta_comp_826_852()

# %% Plot temps from 763/826 ratio vs 763/852 ratio.


def temp_comp_826_852():
    from scipy import stats
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    os.chdir('C:/Users/ivan/Dropbox/Research/ar_OES/DATA/763_852')
    temp_comp = np.loadtxt('temp_compare.dat', delimiter='|')
    slope, intercept, r_value, p_value, std_err = stats.linregress(temp_comp[:, 0], temp_comp[:, 1])
    lin_fit = temp_comp[:, 0] * slope + intercept
    r_squared = r_value ** 2
    fig1 = plt.figure(figsize=(9, 6))
    plt.title("temp_compare")
    plt.xlabel("temp_int (763/826)")
    plt.ylabel("temp_int (763/852)")
    plt.plot(temp_comp[:, 0], temp_comp[:, 1], 'r*')
    plt.plot(temp_comp[:, 0], lin_fit, 'r', label=("R^2= {:6.2f}".format(r_squared)))
    plt.legend(loc=2)
    plt.savefig('temp_compare.png', dpi=300)
    plt.show(fig1)


temp_comp_826_852()

# %%


def plot_ratio_temp():
    import os
    import numpy as np
    from scipy import stats
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    os.chdir('C:/Users/ivan/Dropbox/Research/ar_OES/DATA/763_852')
    all_dat = np.loadtxt('ALLPORTS.dat', delimiter='|')

    def square_fit(xx, aa, bb, cc):
        return (aa * xx**2 + bb * xx + cc)
    popt, pcov = curve_fit(square_fit, all_dat[:, 2], temp_int)
    # print(fit_data)
    fig3 = plt.figure(figsize=(9, 6))
    aa = popt[0]
    bb = popt[1]
    cc = popt[2]
    xx = np.linspace(3.0, 8, 1000)
    fit = aa * xx ** 2 + bb * xx + cc
    perr = np.sqrt(np.diag(pcov))
    print("standard deviation, quadratic fit (a, b, c) :",  perr)
    plt.plot(all_dat[:, 2], temp_int,  'g*')
    plt.plot(xx, fit,  'g')
    plt.title("ALEXIS: Temp vs Line Ratio")
    plt.xlabel("Experimental line ratio 762/852")
    plt.ylabel("Predicted temp (eV) from model")
    plt.xlim(4.5, 8.0)
    plt.show(fig3)


plot_ratio_temp()


# %%
def exp_theor_temp_comp():
    import os
    import numpy as np
    from scipy import stats
    import matplotlib.pyplot as plt
    os.chdir("C:/Users/ivan/Dropbox/Research/ar_OES/DATA/763_852")
    exp_temp_data = np.loadtxt("ALLPORTS.dat", delimiter="|")
    model_dat = np.loadtxt("temp_compare.dat", delimiter='|')

    slope, intercept, r_value, p_value, std_err = stats.linregress(
        model_dat[:, 0], exp_temp_data[:, 0])
    lin_fit = model_dat[:, 0] * slope + intercept
    r_squared = r_value ** 2
    fig1 = plt.figure(figsize=(9, 6))
    # plt.plot(model_dat[:,0], lin_fit, 'r', label = ("R^2= {:6.2f}".format(r_squared)))
    # plt.legend(loc = 1)
    plt.xlim(0,7)
    plt.ylim(0,7)
    plt.plot(model_dat[:, 0], exp_temp_data[:, 0], 'r*')
    plt.savefig('exp_vs_model_temps.png', dpi=300)
    plt.show(fig1)
    return(exp_temp_data, model_dat)

exp_dat, model_dat = exp_theor_temp_comp()


# %%
def exp_temp_ratio():
    import os
    import numpy as np
    from scipy import stats
    import matplotlib.pyplot as plt
    os.chdir("C:/Users/ivan/Dropbox/Research/ar_OES/DATA/763_852")
    exp_temp_data = np.loadtxt("ALLPORTS.dat", delimiter="|")
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        exp_temp_data[:, 2], exp_temp_data[:, 0])
    lin_fit = exp_temp_data[:, 2] * slope + intercept
    r_squared = r_value ** 2

    fig2 = plt.figure(figsize=(9, 6))
    plt.plot(exp_temp_data[:, 2], exp_temp_data[:, 0], 'r*')
    # plt.plot(exp_temp_data[:,2], lin_fit, 'r', label = ("R^2= {:6.2f}".format(r_squared)))
    # plt.legend(loc = 1)
    plt.show(fig2)


exp_temp_ratio()

# %%


def mfit_dens():
    import os
    import numpy as np
    from scipy import stats
    import matplotlib.pyplot as plt
    os.chdir("C:/Users/ivan/Dropbox/Research/ar_OES/DATA/763_852")
    exp_temp_data = np.loadtxt("BC1A.dat", delimiter="|")
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        exp_temp_data[:, 1], exp_temp_data[:, 3])
    r_squared = r_value ** 2
    lin_fit = exp_temp_data[:, 1] * slope + intercept
    fig2 = plt.figure(figsize=(9, 6))
    plt.plot(exp_temp_data[:, 1], exp_temp_data[:, 3], 'r*')
    plt.plot(exp_temp_data[:, 1], lin_fit, 'r', label=("R^2= {:6.2f}".format(r_squared)))
    plt.legend(loc=1)
    plt.savefig('mfit_dens.png', dpi=300)
    print(slope)
    plt.show(fig2)


mfit_dens()

# %%


def mfit_dens():
    import os
    import numpy as np
    from scipy import stats
    import matplotlib.pyplot as plt
    os.chdir("C:/Users/ivan/Dropbox/Research/ar_OES/DATA/763_852")
    exp_temp_data = np.loadtxt("BC1A.dat", delimiter="|")
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        exp_temp_data[:, 1], exp_temp_data[:, 2])
    r_squared = r_value ** 2
    lin_fit = exp_temp_data[:, 1] * slope + intercept
    fig2 = plt.figure(figsize=(9, 6))
    plt.plot(exp_temp_data[:, 1], exp_temp_data[:, 2], 'r*')
    plt.plot(exp_temp_data[:, 1], lin_fit, 'r', label=("R^2= {:6.2f}".format(r_squared)))
    plt.legend(loc=1)
    plt.savefig('mfit_dens.png', dpi=300)
    plt.show(fig2)


mfit_dens()


