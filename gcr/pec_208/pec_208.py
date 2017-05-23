from astropy.io import fits
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np

hdulist = fits.open('pec_dat_650_950_lowdens.fits')
                   
# VARIABLES

# Create numpy arrays from the pages in the .fits file
pec_temps = hdulist[0].data       # The array of temps in eV
n_meta = hdulist[1].data          # Array listing the metastable #
pec_dens = hdulist[2].data        # Density array
pec_wave = hdulist[3].data        # Wavelengths corresponiding to each PEC
pec_pec = hdulist[4].data.T       # 3-D array containing all PEC's

# Various useful things.                 
size_meta = n_meta.size           # Number of metastables
n_pec = pec_wave.size             # number of PEC's for each metastable.
wl_min = int(pec_wave[0])         # smallest wavelength
wl_max = pec_wave[n_pec - 1]      # largest wavelength

linenum = float(input('Desired PEC wavelength = ? ')) # Which PEC ?


# Array of wavelengths for PEC's.                 
wavelist = np.zeros((np.size(pec_wave),2))
for i in range(0,np.size(pec_wave)):
    wavelist[i,0] = i
    wavelist[i,1] = pec_wave[i]
    
# Function to find the nearest value in an array. returns value and index.
def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

# Create dictionaries to hold the invidual meta PEC's and Broadened flux
pec_meta = dict()
for i in range(size_meta):
    pec_meta[i] = pec_pec[:, :, n_pec*i:n_pec*(i+1)]

# Find the nearest PEC wavelength to the desired input.
indx, wavelength = find_nearest(wavelist[:,1], linenum)

print("Nearest PEC wavelength = ", wavelength, " nm")
print("PEC index # = ", indx)

# %% 3D plot of the PEC's for chosen wavelength. \
#    Displayed in increasing metastable order.

X = np.array(pec_temps)
Y = np.array(pec_dens)
X, Y = np.meshgrid(X,Y)


for i in range(size_meta):
    fig1 = plt.figure()
    Z = pec_meta[i][:,:,indx]    
    
    ax = fig1.gca(projection='3d')
    # Plot the surface.
    surf = ax.plot_surface(X, Y, np.transpose(Z), cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    plt.xlabel('Temperature (eV)')
    plt.ylabel('Density (cm-3)')
    plt.show()
