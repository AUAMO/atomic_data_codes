# %%

import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir('c:/Users/ivan/OneDrive/Research/ar_OES')
# all the magic is contained the the qcd_solver library
from qcd_solver import *


# %%

dens_array, temp_array, pop_grd_matrix, pop_1s5_matrix, \
pop_1s3_matrix = qcd_eq_solver('qcd208.pass', 3)

pop1 = pop_grd_matrix[0:25] / pop_1s3_matrix[0:25]
pop2 = pop_grd_matrix[0:25] / pop_1s5_matrix[0:25]
pop3 = pop_1s5_matrix[0:25] / pop_1s3_matrix[0:25]
# %%
t_matrix, td_pop_matrix = qcd_td_solver('qcd208.pass',\
                        nmeta = 3, t_indx = 8,  d_indx = 4)


# %%

fig2 = plt.figure(figsize=(8, 6), facecolor='white')
plt.semilogy(t_matrix, td_pop_matrix[0,:], 'c', label="ground")
plt.semilogy(t_matrix, td_pop_matrix[1,:], 'g', label="1s5")
plt.semilogy(t_matrix, td_pop_matrix[2,:], 'r', label="1s3")

plt.title("ALEXIS: Time dependent metastable populations", \
          weight = 'bold', fontsize=14)
plt.xlabel("Time (s)", fontsize=14)
plt.ylabel(r"Rate Coefficient $(cm^3 * s)^{-3}$", fontsize=14)

plt.legend(loc=4)

plt.savefig("td_pop.png", dpi=300)

plt.show(fig2)


fig2 = plt.figure(figsize=(8, 6), facecolor='white')
plt.plot(t_matrix, td_pop_matrix[1,:] / td_pop_matrix[2,:], 'b')

plt.title(r"ALEXIS: 1s5/1s3 time dependent metastable fraction", \
          weight = 'bold', fontsize=14)
plt.xlabel("Time (s)", fontsize=14)
plt.ylabel("1s5/1s3 metastable fraction", fontsize=14)

plt.legend(loc=4)

plt.savefig("td_frac.png", dpi=300)

plt.show(fig2)

therm = np.sqrt(1.36e-23 * 300 / 6.67e-26)
print('therm=', therm)
print(400*0.0007)
# %%

fig2 = plt.figure(figsize=(8, 6), facecolor='white')
ratio_mat = np.divide(td_pop_matrix[1,1:], \
                      td_pop_matrix[2,1:])

grd_ratio_mat =np.divide(td_pop_matrix[0,1:], \
                      td_pop_matrix[2,1:])
plt.plot(t_matrix[1:], grd_ratio_mat)
plt.title("Time dependent metastable fraction (ground/1s3)", \
          weight='bold', fontsize=10)
plt.xlabel("Time (s)", fontsize=14)
plt.ylabel("Ratio (ground/1s3)", fontsize=14)

plt.savefig("ground_ratio.png", dpi=300)

plt.show(fig2)

indx = find_nearest(ratio_mat, 13)


# %%
def qcd_temp_plot(d_indx):
    import matplotlib.pyplot as plt
    fig2 = plt.figure(figsize=(8, 6), facecolor='white')
    plt.plot(temp_array,  pop_1s5_matrix[:,d_indx], 'g')
    plt.plot(temp_array,  pop_1s3_matrix[:,d_indx], 'r')
    plt.show(fig2)
qcd_temp_plot(3)

# %%
def qcd_dens_plot(t_indx):
    import matplotlib.pyplot as plt
    fig1 = plt.figure(figsize=(8, 6), facecolor='white')
    plt.plot(dens_array, pop_1s3_matrix[t_indx,:], 'g')
    plt.plot(dens_array, pop_1s5_matrix[t_indx,:], 'r')
    plt.show(fig1)
qcd_dens_plot(8)

# %%

pop_grd_relground = np.divide(pop_grd_matrix, pop_grd_matrix)
pop_1s5_relground = np.divide(pop_1s5_matrix, pop_grd_matrix)
pop_1s3_relground = np.divide(pop_1s3_matrix, pop_grd_matrix)

np.savetxt('pop_1s5_relground.txt',pop_1s5_relground,\
           fmt = '%6.3e')
np.savetxt('pop_1s3_relground.txt',pop_1s3_relground,\
           fmt = '%6.3e')
pop_grd_matrix


# %%


print(dens_label)
fig1 = plt.figure(figsize=(8,6), facecolor="white")
plt.title("Neutral argon: Metastable fractions at equlibrium", \
          weight = 'bold', fontsize=14)
plt.xlabel("Electron Temperature (eV)", fontsize=14)
plt.ylabel(r"Rate Coefficient $(cm^3 * s)^{-3}$", fontsize=14)

plt.semilogy(temp_array, pop_grd_matrix[:,0], 'c+', label="ground")
plt.semilogy(temp_array, pop_1s5_matrix[:,0], 'gx', label="1s5")
plt.semilogy(temp_array, pop_1s3_matrix[:,1], 'rx', label="1s3")

plt.semilogy(temp_array, pop_grd_matrix[:,1:19], 'c+')
plt.semilogy(temp_array, pop_1s5_matrix[:,1:19], 'gx')
plt.semilogy(temp_array, pop_1s3_matrix[:,1:19], 'rx')

plt.legend(loc=4)
plt.xlim(0.45,20)
plt.ylim(1e-5,5)
plt.savefig("equib_frac.png", dpi=300)
plt.show(fig1)
