
import numpy as np
filename = 'qcd208.pass'
nmeta = 3

# %% function to read all qcd blocks in to array and store in a dict
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
    return idx  
