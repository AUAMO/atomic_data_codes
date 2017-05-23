#####################################################################################################
# %% Use  experimental line ratio and metastable fraction from the linear fit from allports plot to
#    calculate temperature.
def ar_temp_finder(temp_ratio, density):
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    os.chdir('C:/Users/ivan/OneDrive/Research/ar_OES/DATA')

    def find_nearest(array, value):
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

    gscale = 2000
    # slope, intercept =  2.86026783707 -7.82438767561
    mscale =  2.86026783707  * temp_ratio - 7.82438767561
    # mscale = 0.249341860455 * line_ratio + 3.64323291233

    lower, lval = find_nearest(pec_wave, 763.5)
    upper, uval = find_nearest(pec_wave, 852.5)

    pec_ratio_weighted = (gscale * pec_meta[0][:, :, lower] + mscale \
        * pec_meta[1][:, :, lower] + pec_meta[2][:, :, lower]) / (gscale * \
        pec_meta[0][:, :, upper] + mscale * pec_meta[1][:, :, upper] \
        + pec_meta[2][:, :, upper])

    pec_array = np.mean(pec_ratio_weighted, axis=1)
    from scipy.interpolate import interp1d
    pec_fit = interp1d(pec_temps, pec_array, kind = 'cubic')

    indx, val = find_nearest(pec_array, temp_ratio)
    pec_x = np.linspace(0.5,20,1000)
    indx2, val2 = find_nearest(pec_fit(pec_x), temp_ratio)
    pec_x = np.linspace(0.5,20,1000)
    indx2, val2 = find_nearest(pec_fit(pec_x), temp_ratio)

    temp = pec_x[indx2]

    temp_indx, temp_val = find_nearest(pec_temps, temp)



    pec_dens_weighted = (gscale * pec_meta[0][temp_indx, :, lower] + mscale \
        * pec_meta[1][temp_indx, :, lower] + pec_meta[2][temp_indx, :, lower]) / (gscale * \
        pec_meta[0][temp_indx, :, upper] + mscale * pec_meta[1][temp_indx, :, upper] \
        + pec_meta[2][temp_indx, :, upper])



    # Uncomment to generate plots ==============================================
    fig2 = plt.figure(figsize=(9,4), facecolor='white')
    plt.title('Weighted Line Ratio: mscale = {:6.2f}'.format(mscale))
    plt.xlabel('Temperature (eV)')
    plt.plot(pec_temps, pec_array, 'r*')
    plt.plot(pec_x, pec_fit(pec_x), 'r')
    plt.ylabel('Weighted Line Ratio')
    plt.axhline(temp_ratio, color="red")
    plt.axvline(pec_x[indx2], label = 'Temp = {:6.2f}'.format(pec_x[indx2]), color="orange")
    plt.legend(loc = 4)
    plt.xlim(0, 20)
    plt.ylim(0,10)
    plt.savefig("temp_finder_" + str(temp_ratio) + "_763_826.png", dpi=300)
    plt.show(fig2)

#    fig2 = plt.figure(figsize=(9,4), facecolor='white')
#    plt.title('Density: mscale = {:6.2f}'.format(mscale))
#    plt.xlabel('Density (m^-3)')
#    plt.plot(pec_dens, pec_dens_weighted, 'r*')
#    # plt.plot(pec_x, pec_fit(pec_x), 'r')
#    plt.ylabel('Weighted Line Ratio: Density')
#    # plt.axhline(dens_ratio)
#    # plt.axvline(pec_x[indx2], label = 'Temp = {:6.2f}'.format(pec_x[indx2]))
#    # plt.legend(loc = 4)
#    # plt.xlim(0, 20)
#    # plt.ylim(0,10)
#    plt.savefig("temp_finder_" + str(dens_ratio) + "_763_826.png", dpi=300)
#    plt.show(fig2)
    # ==========================================================================
    print('temp_indx = ', temp_indx)
    label = ("line ratio, temp, m_frac  =" " {:6.2f}".format(temp_ratio) + " ; {:6.2f} (eV)".format(pec_x[indx2]) + \
           " ; {:6.2f}".format(mscale) )
    print(label)
    return(temp, mscale, temp_ratio, label)

temp, mscale, temp_ratio, label = ar_temp_finder(5.96, 2.0)
