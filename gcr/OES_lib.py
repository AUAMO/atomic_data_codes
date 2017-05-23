import numpy as np

# %% open and extract spectra and temp_dens data from an ALEXIS generated fits file (with GUI).
#    Argument is the directory where spectrum is located.
#    return (temp_dens_data, temp_dens_av, wavelengths, spect_array, spect_av)
#    temp_dens_data: Double probe data for the range r = 0 to r = 50 mm in 2 mm increments
#                    data in in np.array(data, radial position),
#                    data = [['radius'],['t_e'],['n_e'],['isat'],['di_dv0']]
#    temp_dens_ave: np.array of average radial position, average temp, etc...
#    wavelengths: wavelength array with n elements
#    spect_array: np.array(10,n); line of sight spectral data taken at times corresponding
#                 last ten probe measurements
#    spect_av: average spectral data
def retrieve_ALEXIS_data(spect_dir):
    from astropy.io import fits
    from scipy.optimize import curve_fit
    from scipy import stats
    import os
    os.chdir(spect_dir)
    import datetime; now = datetime.datetime.now()
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

        spect_av = spect.mean(axis=0)
        return (wavelengths, spect_array, spect_av)
    wavelengths, spect_array, spect_av = spect_get()

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
        temp_dens_array = temp_dens_array
        return temp_dens_array


    data = temp_dens_get(0)
    for i in range(1,26):
        tmp_data = data
        data = temp_dens_get(i)
        temp_dens_data = np.hstack((tmp_data, data))
        data = temp_dens_data
    temp_dens_av = temp_dens_data.mean(axis=1)


    return (filename, temp_dens_data, temp_dens_av, wavelengths, spect_array, spect_av)

# %% From argon spectra extract the line ratios (763/826 ; 826/852 ; 853/919).
#    Arguments are two 1D array for wavelengths and spectra of matching size.
#    return (ratio_763_826, ratio_826_852, ratio_852_919)
def get_spect_ratios(wavelengths, spect):
    import peakutils
    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return idx
    values = np.array([763.5, 852.1, 919.0])

    tmp_ratio = np.zeros(3)
    tmp_value = np.zeros(4).astype(int)

    # temp_spect = spect[877:1614]    ## Uncomment these lines for  negligible speed increase
    # temp_wavelengths = wavelengths[877:1614]  ## Adjust slice to focus on region of interest in spectrum.
    peaks = peakutils.indexes(spect, thres = 0.05/max(spect), min_dist = 5)
    temp_spect = spect[peaks]
    temp_spect = np.reshape(temp_spect,(np.size(peaks),1))
    temp_wavelengths = wavelengths[peaks]
    for j in range(np.size(values)):
        tmp_value[j] = find_nearest(temp_wavelengths, values[j])
    for k in range(np.size(tmp_ratio)):
        tmp_ratio[k] = temp_spect[tmp_value[k]] / temp_spect[tmp_value[k+1]]
    ratio_763_852 = tmp_ratio[0] * 1.05305396186 / 1.04566612136
    ratio_852_919 = tmp_ratio[1] * 1.04566612136 / 1.04267580497

    return (ratio_763_852, ratio_852_919)


def get_meta_frac(pec_file, pec_ratio, wl_low, wl_high):
    from astropy.io import fits
    # Import the .fits file which contains all the pec dat and
    # extract the information.
    hdulist = fits.open(pec_file)
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

    pec_meta = dict()
    for i in range(size_meta):
        pec_meta[i] = pec_pec[:, :, n_pec*i:n_pec*(i+1)]


    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return idx, array[idx]

    lower, lval = find_nearest(pec_wave, wl_low)
    upper, uval = find_nearest(pec_wave, wl_high)
    gscale = 2000
    mscale = np.linspace(0.1, 30, 300)
    pec_av = np.zeros(300)

    for i in range(np.size(mscale)):
        pec_ratio_weighted = (gscale * pec_meta[0][:, :, lower] + mscale[i] \
            * pec_meta[1][:, :, lower] + pec_meta[2][:, :, lower]) / (gscale * \
            pec_meta[0][:, :, upper] + mscale[i] * pec_meta[1][:, :, upper] \
            + pec_meta[2][:, :, upper])
        pec_av[i] = np.sqrt(np.mean(np.square(pec_ratio_weighted[0:10])))

    r_indx, r_val = find_nearest(pec_av, pec_ratio)
    m_ratio = mscale[r_indx]
    return (lval, uval, r_val, m_ratio)



# %% function to read all qcd blocks in to array and store in a dict.
#    Arguments are filenmae of qcd file and number of metastables.
#    return (qcd_dat)
#    Returns python dictionary with keys 'qcd_ij', where i and j are metastable indices.
def qcd_reader(filename, nmeta):
    with open(filename,'r') as f:
        qcd_dat = {}
        nblocks = nmeta * (nmeta - 1)
        for block in range(nblocks):
              line_dat = f.readline()
              while line_dat[0] != '=':
                  line_dat = f.readline()

              line_dat = f.readline()
              ii = int(line_dat[6:8])
              jj = int(line_dat[53:55])
              name = 'qcd_' + str(ii) + '_' + str(jj)
              f.readline()
              line_dat = f.readline()
              ntemps = int(line_dat[8:10])
              ndens = int(line_dat[13:15])

              for line in range(3):
                  f.readline()

              qcd_array = np.zeros((ntemps, ndens))

              dens_quot, dens_rem =  divmod(ndens,8)
              dens_lines = dens_quot
              if dens_rem > 0:
                  dens_lines = dens_lines +1

              temp_quot, temp_rem =  divmod(ntemps,8)
              temp_lines = temp_quot
              if temp_rem > 0:
                  temp_lines = temp_lines + 1

              dens_array = np.array(0)
              for i in range(dens_lines):
                  dat = str.strip(f.readline()).replace('D','E').split()
                  dat = np.array(dat)
                  dens_array = np.hstack((dens_array,dat))
              temp_array = np.array(0)
              for i in range(temp_lines):
                  dat = str.strip(f.readline()).replace('D','E').split()
                  dat = np.array(dat)
                  temp_array = np.hstack((temp_array,dat))
              temp_array = temp_array[1:36]
              f.readline()

              qcd_array = dens_array
              for i in range(ntemps):
                  ldat = qcd_array
                  cdat = np.array(0)
                  for j in range(dens_lines):
                    dat = str.strip(f.readline()).replace('D','E').split()
                    dat = np.array(dat)
                    cdat = np.hstack((cdat,dat))
                  qcd_array = np.vstack((ldat,cdat))
              qcd_array[1:36,0] = temp_array
              qcd_array = qcd_array.astype(np.float)
              qcd_array[1:36,0] = qcd_array[1:36,0] / 11604.505
              qcd_dat[name] = qcd_array
        print('\n', 'ntemps = ', ntemps, "; ndens = ", ndens)
        print(' shape of each qcd array = ', np.shape(qcd_array), '\n')
        return (qcd_dat)

# %%  function to surface plot the qcd's
def qcd_plotter(filename, nmeta):
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter
    qcd_dat = qcd_reader(filename, nmeta)
    qcd_keys = list(qcd_dat.keys())
    print(qcd_keys)
    for i in range(np.size(qcd_keys)):
        fig = plt.figure(figsize=(8, 6), facecolor='white')
        ax = fig.add_subplot(111, projection='3d')

        XX = qcd_dat[qcd_keys[i]][1:36,0]
        YY = qcd_dat[qcd_keys[i]][0,1:11]
        XX, YY = np.meshgrid(XX,YY)
        ZZ = qcd_dat[qcd_keys[i]][1:36,1:11]
        ax.set_title(qcd_keys[i])
        ax.set_xlabel('Temp (eV)', fontsize=8)
        ax.set_ylabel('Density (cm^-3)', fontsize=8)
        plt.setp(ax.get_zticklines(), visible=False)
        plt.setp(ax.get_yticklines(), visible=False)
        plt.setp(ax.get_xticklines(), visible=False)
        ax.plot_surface(XX, YY, np.transpose(ZZ), cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%4.0e'))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%4.1f'))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%4.2e'))

        plt.savefig('pec_ratio.png', dpi=300)
        plt.show(fig)
# %%

def qcd_eq_solver(filename, nmeta):
    qcd_dat = qcd_reader(filename, nmeta)
    qcd_keys = list(qcd_dat.keys())
    qcd_keys

    temp_array = qcd_dat['qcd_1_2'][1:,0]
    dens_array = qcd_dat['qcd_1_2'][0,1:]
    ntemps = np.size(temp_array)
    ndens = np.size(dens_array)
    pop_grd_matrix = np.zeros((ntemps, ndens))
    pop_1s5_matrix = np.zeros((ntemps, ndens))
    pop_1s3_matrix = np.zeros((ntemps, ndens))

    for i in range(ntemps):
        for j in range(ndens):
            # print(i, j)
            q_12 = qcd_dat['qcd_1_2'][i+1][j+1] * dens_array[j]  # ; print(q_12)
            q_13 = qcd_dat['qcd_1_3'][i+1][j+1] * dens_array[j]  # ; print(q_13)
            q_21 = qcd_dat['qcd_2_1'][i+1][j+1] * dens_array[j]  # ; print(q_21)
            q_23 = qcd_dat['qcd_2_3'][i+1][j+1] * dens_array[j]  # ; print(q_23)
            q_31 = qcd_dat['qcd_3_1'][i+1][j+1] * dens_array[j]  # ; print(q_31)
            q_32 = qcd_dat['qcd_3_2'][i+1][j+1] * dens_array[j]  # ; print(q_32)
            rate_matrix = np.array([[1,1,1], [q_12, -(q_21 + q_23), q_32], \
                                   [q_13, q_23, -(q_31 + q_32)]])
            pop_matrix = np.zeros([3,1])
            b_matrix = [[1],[0],[0]]

            rate_inv = np.linalg.inv(rate_matrix)
            pop_matrix = np.dot(rate_inv,b_matrix)
            # print(np.shape(pop_grd_matrix))

            pop_grd_matrix[i][j] = pop_matrix[0]
            pop_1s5_matrix[i][j] = pop_matrix[1]
            pop_1s3_matrix[i][j] = pop_matrix[2]
    return(dens_array, temp_array, pop_grd_matrix, pop_1s5_matrix, pop_1s3_matrix)

# %%
def qcd_td_solver(filename, nmeta, t_indx, d_indx):
    qcd_dat = qcd_reader(filename, nmeta)
    temp_array = qcd_dat['qcd_1_2'][1:,0]
    dens_array = qcd_dat['qcd_1_2'][0,1:]
    ntemps = np.size(temp_array)
    ndens = np.size(dens_array)
    qcd_dat = qcd_reader(filename, nmeta)
    qcd_keys = list(qcd_dat.keys())

    q_12 = qcd_dat['qcd_1_2'][t_indx][d_indx] * dens_array[d_indx]  # ; print(q_12)
    q_13 = qcd_dat['qcd_1_3'][t_indx][d_indx] * dens_array[d_indx]  # ; print(q_13)
    q_21 = qcd_dat['qcd_2_1'][t_indx][d_indx] * dens_array[d_indx]  # ; print(q_21)
    q_23 = qcd_dat['qcd_2_3'][t_indx][d_indx] * dens_array[d_indx]  # ; print(q_23)
    q_31 = qcd_dat['qcd_3_1'][t_indx][d_indx] * dens_array[d_indx]  # ; print(q_31)
    q_32 = qcd_dat['qcd_3_2'][t_indx][d_indx] * dens_array[d_indx]  # ; print(q_32)

    rate_matrix = np.array([[-(q_12 + q_13), q_21, q_31], \
                            [q_12, -(q_21 + q_23), q_32], \
                            [q_13, q_23, -(q_31 + q_32)]])
    qcd_max = np.max(rate_matrix)

    b_matrix = np.array([[1],[0],[0]])
    t_max = 1.0
    t_min = 0.0
    delta_t = 1./(10 * qcd_max)
    nsteps = (t_max - t_min) / delta_t
    t_matrix = np.linspace(t_min, t_max, nsteps)
    np.max([qcd_dat[i] for i in qcd_dat])

    td_pop_matrix = np.zeros((3, int(nsteps)))
    td_pop_matrix[0,0] = 1.

    for i in range(1, np.size(t_matrix)):
        tmp_matrix = delta_t * np.dot(rate_matrix, td_pop_matrix[:,i-1])
        td_pop_matrix[:,i] = td_pop_matrix[:,i-1] + tmp_matrix
        if (i > 100):
            chuckles = np.var(td_pop_matrix[0,i-100:i])
            if(chuckles < 1e-10):
                break
    t_steady = i * delta_t
    print('time step: ', delta_t)
    print('time to steady state: ', t_steady)
    print('number of steps: ', i)
    td_pop_matrix = td_pop_matrix[:,0:i]
    t_matrix= t_matrix[0:i]
    return (t_matrix, td_pop_matrix)

#%%
def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return (idx, array(idx))

# %%
def get_meta_scale(temp, ratio, lval, uval):

    from astropy.io import fits

    def find_nearest(array,value):
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
        # flux[i] = np.zeros(n_wavel)

    gscale = 2000
    mscale = np.linspace(0.1, 30, 300)

    # pec_av = np.zeros(300)
    lower, lval = find_nearest(pec_wave, lval)
    upper, uval = find_nearest(pec_wave, uval)

    pec_meta[0] = np.mean(pec_meta[0], axis = 1)
    pec_meta[1] = np.mean(pec_meta[1], axis = 1)
    pec_meta[2] = np.mean(pec_meta[2], axis = 1)

    pec_scaled_array = list()
    for i in range(np.size(mscale)):
        pec_ratio_weighted = (gscale * pec_meta[0][:, lower] + mscale[i] \
             * pec_meta[1][:, lower] + pec_meta[2][:, lower]) / (gscale * \
             pec_meta[0][:, upper] + mscale[i] * pec_meta[1][:, upper] \
             + pec_meta[2][:, upper])
        pec_scaled_array.append(pec_ratio_weighted)
    pec_scaled_array = np.array(pec_scaled_array)

    indx, val = find_nearest(pec_temps, temp)
    if (temp - val) > 0.0:
        pec_indx_low = indx
        pec_indx_high = indx+1
        pec_temp_low = pec_temps[indx]
        pec_temp_high = pec_temps[indx+1]
        # print("tag")
    else:
        pec_indx_low = indx - 1
        pec_indx_high = indx
        pec_temp_low = pec_temps[indx - 1]
        pec_temp_high = pec_temps[indx]
        # print("you're it")

    min_vals = list()
    r_list = list()
    for i in range(np.size(mscale)):
        indx, val = find_nearest(pec_scaled_array[i], ratio)
        if (indx == pec_indx_high) or (indx == pec_indx_low):
            x = [pec_temp_low, pec_temp_high]
            y = [pec_scaled_array[i, pec_indx_low], pec_scaled_array[i, pec_indx_high]]
            slope = (y[1] - y[0]) / (x[1] - x[0])
            midpoint = (pec_temp_high + pec_temp_low) / 2
            intercept = y[0] - slope * x[0]
            mid_ratio = slope * midpoint + intercept
            dif_low = np.abs(ratio - y[0])
            dif_mid = np.abs(ratio - mid_ratio)
            dif_hi = np.abs(ratio  - y[1])
            if (dif_low > dif_mid) and (dif_low > dif_hi):
                dif_ratio = slope * x[0] + intercept
                rtag = 'dif_low'
            elif (dif_mid > dif_low) and (dif_mid > dif_hi):
                dif_ratio = slope * midpoint + intercept
                rtag = 'dif_mid'
            else:
                dif_ratio = slope * x[1] + intercept
                rtag = 'dif_hi'
            min_vals.append([i,ratio - dif_ratio, y[0], y[1]])
            r_list.append([i, rtag])
    min_vals = np.array(np.abs(min_vals))


    min_indx, min_val = find_nearest(min_vals[:,1], 0.0)
    meta_indx = int(min_vals[min_indx,0])
    meta_scale = mscale[meta_indx]
    pec_array = pec_scaled_array[meta_indx]

    return (meta_indx, meta_scale, pec_temps, pec_array)
