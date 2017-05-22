#! /usr/bin/python3
""" Author: Ivan
This script is designed to extract branching ratios from an adf04 file

"""
import matplotlib.pyplot as plt
from read_adf04 import read_adf04
import numpy as np
from branch import csec

# initialize some variables
ZZERO = 1e-10
LAST_LEVEL = 57
NTRANS = 1596
NERGS = 12
NTEMPS = 964

DXL = 1
DXU = 15
IMETA = 1

DXL = int(input('Lower level = ? '))
DXU = int(input('Upper level = ? '))

# %%
# Get the index values needed for the branching ratio code.
INDX_INIT = DXL - 1.0
INDX_COL = INDX_INIT * (LAST_LEVEL - 1) - \
    (INDX_INIT * (INDX_INIT - 1.0)) / 2. + 1.
INDX_TRANS = INDX_COL + (DXU - DXL) - 1.
INDX_INIT = INDX_INIT + 1.0
INDX_TRANS_FINAL = INDX_INIT * \
    (LAST_LEVEL - 1) - (INDX_INIT * (INDX_INIT - 1.0)) / 2.

# Read in the adf04 file.
DATA = read_adf04(LAST_LEVEL, NERGS, 'adf04')

# Initialize a numpy array to store A-values.
AVAL_MAT = np.zeros((NTRANS, 3))
AVAL_MAT[:, 0] = DATA.lower
AVAL_MAT[:, 1] = DATA.upper
AVAL_MAT[:, 2] = DATA.aval

# Initialize array to store the sum off all transitional A-values
# for each level in adf04.
ATOT = np.zeros(np.size(DATA.wa) + 1)
for i in range(1, np.size(DATA.wa) + 1):
    ATOT[i] = sum(DATA.aval[DATA.upper == i])

# This array stores branching ratios for every transition in DATA_mat.
# BRAT(i-->j) = AVAL(i-->j) / sum(AVAL(i-->k))
# ratio_mat = AVAL_MAT[(AVAL_MAT[:, 0] >= float(DXU)) & \
#             (AVAL_MAT[:, 1] >= float(DXU + 1)), :]
ratio_mat = AVAL_MAT[(AVAL_MAT[:, 0] >= float(DXU)), :]

for i in range(0, np.size(ratio_mat[:, 1]) - 1):
    ratio_mat[i, 2] = ratio_mat[i, 2] / ATOT[int(ratio_mat[i, 1])]
# %%

# Now the tough part. For the level above direct excitation (DXU + 1)
# through LAST_LEVEL we want to create branching FRACTIONS. As an example,
# Let's assume that we are exciting from level 2 into level 10 (DXL=2;DXU=10)
# and there are only 13 total levels.
#
# BEGIN LOOP
# temp_mat1 = [ 10  11  Br_10,11]
# temp_mat2 = [ 11  12  Br_11,12
# 	            11  13  Br_11,13]
#
# calculate new branching ratios, ' indicates a newly calculated branch

# temp_mat2 = [ 10  12  Br_10,12'
#               10  13  Br_10,13']
#
# store branching ratios into DATA_mat
#
# DATA_mat =  [ 10  11   Br_10,11
# 	            10  12   Br_10,12'
# 		        10  13   Br_10,13']
#
# the next loop requires temp_mat1 has an upper level of 12
# store everything from DATA_mat with an upper level of 12
# in temp_mat3, assuming that these new ratios are gt zzero
#
# temp_mat3 =  [10  12   Br_10,12']
# LOOP

# temp_mat1 will now check temp_mat3 AND ratio_mat for upper levels = 12
#
# temp_mat1 =  [10  12   Br_10,12     <------ this one's in ratio_mat
#               10  12   Br_10,12']   <------ this one's in temp_mat3
#                                                  from previous loop
# temp_mat2 =  [12  13   Br_12,13]
#
# caluculate new branching ratios
#
# temp_mat2 = [ 10  13  Br_10,13''
# 	            10  13  Br_10,13'']
#
# cat ratios from last loop with new  branching ratios into DATA_mat
#
# DATA_mat = [ 10  11   Br_10,11
# 	           10  12   Br_10,12'
# 		       10  13   Br_10,13'
#              10  13   Br_10_13''   <------ These two are DIFFERENT branches
#              10  13   Br_10,13'']  <------ calculated in loop 2
#
# temp_mat3 = [ 10  13   Br_10,13'
#               10  13   Br_10,13''
#               10  13   Br_10,13'']
#
# LOOP   (this last loop only pulls out the last braching ratio
#                                            --> Br_10,last_level
# DATA_mat =  [ 10  11   Br_10,11
#               10  12   Br_10,12
#               10  13   Br_10,13
# 	            10  12   Br_10,12'
# 		        10  13   Br_10,13''
#               10  13   Br_10,13'']
#
# There is only one branch from 10 to 11
# There are two branches from 10 to 12
# There are four branches from 10 to 13
#
# 11  12  13   12  13  13  13
# 10  10  10   11  12  11  12
#              10  10  10  11
#                          10
#
# DATA_mat stores all of them.
#==============================================================================

temp_mat1 = np.zeros((1, 3))
temp_mat2 = np.zeros((1, 3))
temp_mat3 = np.zeros((1, 3))
DATA_mat = np.zeros((1, 3))

for i in range(1, LAST_LEVEL - DXU + 1):
    temp_mat1 = 0.
    temp_mat2 = 0.
    tmp = ratio_mat[(ratio_mat[:, 0] == DXU) & (ratio_mat[:, 1] == DXU + i), :]
    temp_mat1 = np.concatenate((tmp, temp_mat3), axis=0)
    temp_mat2 = ratio_mat[(ratio_mat[:, 0] == DXU + i), :]

#    temp_mat1 = temp_mat1[np.argsort(temp_mat1[:, 1])]
#    temp_mat2 = temp_mat2[np.argsort(temp_mat2[:, 1])]

    for j in range(0, np.size(temp_mat1[:, 0])):
        for k in range(0, np.size(temp_mat2[:, 0])):
            temp_mat2[k, 2] = temp_mat2[k, 2] * temp_mat1[j, 2]
        temp_mat2[:, 0] = temp_mat1[0, 0]

    DATA_mat = np.concatenate(
        (temp_mat1[0, :].reshape(1, 3), temp_mat2, DATA_mat), axis=0)

    temp_mat3 = 0.
#
    temp_mat3 = DATA_mat[((DATA_mat[:, 1] == DXU + i + 1) &
                          (DATA_mat[:, 2] >= ZZERO)), :]
#    temp_mat3 = DATA_mat[(DATA_mat[:,1] == DXU+i+1),:]
    temp_mat3 = temp_mat3[np.argsort(temp_mat3[:, 1])]

    DATA_mat = DATA_mat[(DATA_mat[:, 2] >= ZZERO), :]
    DATA_mat = DATA_mat[np.argsort(DATA_mat[:, 1])]

#%%

# Sort the DATA array by uper level
DATA_mat = DATA_mat[np.argsort(DATA_mat[:, 0])]
# print(DATA_mat)


# DATA_mat stores all the various non-zero branching fractions. Sum all
# branching fractions with a common upper level and store in branch_mat.
branch_mat = np.zeros([LAST_LEVEL - DXU + 1, 2])
branch_mat[0, 1] = 1.0
branch_mat[0, 0] = DXU
for i in range(1, LAST_LEVEL - DXU + 1):
    branch_mat[i, 0] = i + DXU
    tmp_mat = DATA_mat[(DATA_mat[:, 1] == DXU + i), :]
    branch_sum = sum(tmp_mat[:, 2])
    branch_mat[i, 1] = branch_sum

# %%

level_dat = np.loadtxt('paschen.dat', dtype=str)
for i in range(0, np.size(level_dat[:, 0])):
    for j in range(0, np.size(level_dat[0, :])):
        mystring = level_dat[i, j]
        tick = len(mystring) - 1
        level_dat[i, j] = mystring[2:tick]

IMETA = int(input('Metastable index = ?'))

filename4 = "bcon_" + str(DXL) + "_" + str(DXU) + "_" + level_dat[DXL - 1, 1] \
    + "_" + level_dat[DXU - 1, 1]
branch_condensed = branch_mat[(branch_mat[:, 1] > 1e-10), :]
branch_condensed = branch_condensed[branch_condensed[:, 1].argsort()[::-1]]
branch_condensed = branch_condensed[1:, :]
bcon_indx = branch_condensed[:, 0].astype(int)
bcon_string = np.empty(np.size(bcon_indx), dtype=str)
with open(filename4, 'wb') as f:
    f.write(bytes("#" + filename4 + "\n", "UTF-8"))
    for i in range(0, np.size(bcon_indx)):
        con_string = str(bcon_indx[i]) + '  ' + "{:.4e}".format(branch_condensed[i, 1]) \
            + '  ' + str(DATA.config_data[bcon_indx[i] - 1]) + '  ' \
            + str(DATA.xja[bcon_indx[i] - 1]) + '  ' \
            + "{:.4f}".format(DATA.wa[bcon_indx[i] - 1] * 1.2398e-4)
        f.write(bytes(con_string + "\n", "UTF-8"))

# Store the branching fractions in a file compatible with cross-section
# code 'branch.f'
header1 = "# NTRMN  = " + str(int(INDX_TRANS))
header2 = "# NTRANS = " + str(int(INDX_TRANS_FINAL))
filename = "bratio_" + str(DXL) + "_" + str(DXU) + "_" + level_dat[DXL - 1, 1] + "_" + \
    level_dat[DXU - 1, 1]

with open(filename, 'wb') as f:
    f.write(bytes(header1 + "\n" + header2 + "\n", "UTF-8"))
    np.savetxt(f, branch_mat, fmt='%i +%.6e')
f.close()
with open('bratin', 'wb') as f:
    f.write(bytes(header1 + "\n" + header2 + "\n", "UTF-8"))
    np.savetxt(f, branch_mat, fmt='%i +%.6e')
f.close()

CSEC_DATA = np.asarray(
    csec(int(INDX_TRANS), int(INDX_TRANS_FINAL), IMETA, NTEMPS))
# %%

CSEC_DATA = CSEC_DATA.transpose()

tmp_up = len(level_dat[DXU - 1, 2])
JUP = level_dat[DXU - 1, 2]
JUP = JUP[tmp_up - 4:tmp_up - 3]
JUP = int(JUP)
tmp_down = len(level_dat[DXL - 1, 2])
JDOWN = level_dat[DXL - 1, 2]
JDOWN = JDOWN[tmp_down - 4:tmp_down - 3]
JDOWN = int(JDOWN)
DELTAJ = abs(JUP - JDOWN)

title = level_dat[DXL - 1, 1] + r'$\rightarrow$' + \
    level_dat[DXU - 1, 1] + "  ($\Delta$J = " + str(DELTAJ) + ")"
filename2 = level_dat[DXL - 1, 1] + "_" + \
    level_dat[DXU - 1, 1] + "_J" + str(DELTAJ)


min_jump = 30
threshold = 5
min_erg = 4
for i in range(1, NTEMPS - 1):
    test1 = np.abs(CSEC_DATA[i, 1] - CSEC_DATA[i - 1, 1])
    test2 = np.abs(CSEC_DATA[i, 1] - CSEC_DATA[i + 1, 1])
    test3 = test1 - test2
    if test1 > min_jump and test3 < threshold and CSEC_DATA[i, 0] > min_erg:
        CSEC_DATA[i, 1] = CSEC_DATA[i - 1, 1]
for i in range(1, NTEMPS - 1):
    test1 = np.abs(CSEC_DATA[i, 2] - CSEC_DATA[i - 1, 2])
    test2 = np.abs(CSEC_DATA[i, 2] - CSEC_DATA[i + 1, 2])
    test3 = test1 - test2
    if test1 > min_jump and test3 < threshold and CSEC_DATA[i, 0] > min_erg:
        CSEC_DATA[i, 2] = CSEC_DATA[i - 1, 2]

# %%
# Plot the cross section data
fig2 = plt.figure(figsize=(9, 6), facecolor='white')
plt.plot(CSEC_DATA[:, 0], CSEC_DATA[:, 1], 'r')
plt.plot(CSEC_DATA[:, 0], CSEC_DATA[:, 2], 'b--')
plt.plot(CSEC_DATA[:, 0], CSEC_DATA[:, 3], 'g-.')

plt.title(title)
plt.xlabel('Temperature (eV)')
plt.ylabel('cross-section (MB)')
plt.savefig(filename2 + ".png", format='png', dpi=300)
plt.show(fig2)

with open(filename2 + ".dat", 'wb') as f:
    f.write(bytes(header1 + "\n" + header2 + "\n", "UTF-8"))
    np.savetxt(f, CSEC_DATA)

# %%
print('\n', 'Branching ratio output saved to ', filename4)
print('\n', 'Cross section output saved to ', filename2 + ".dat")
print('\n', 'Plot output saved to ', filename2 + ".png")
