import numpy as np

def reject_outliers(data, column, m=2):
    return data[(abs(data[:,column] - np.mean(data[:,column]))\
     < m * np.std(data[:,column]))]

def read_adf04(nlevels, num_ergs, filename):
        # If inf_on > 1 the file uses an infinite energy point
        ## User inputs
        ntrans = int((nlevels * (nlevels-1))/2)  # Number of total transitions
                                                 # (nlevels * nlevels-1)/2

        ## We create these datatypes to help us strip the avals and config info from the adf04 file
        ## The creators of this file thought it wise to store this info in one long string??
        dt_config = np.dtype('U24')
        dt_aval = np.dtype('U7')

        ## We'll store all the information in a python dict {}
        ## The elements of the dict ('keys') can be accessed by adf04_dat['key']
        adf04_dat ={}

        f = open(filename,'r')

        ## The first line has ion information, which we store in adf04_dat['header']
        adf04_dat['header'] = str.strip(f.readline())

        ## We need to create arrays to store the data from the adf04 file
        adf04_dat['config_data'] = np.zeros(nlevels,dtype=dt_config)
        adf04_dat['isa'] = np.zeros(nlevels,dtype='int')
        adf04_dat['ila'] = np.zeros(nlevels,dtype='int')
        adf04_dat['xja'] = np.zeros(nlevels,dtype='float')
        adf04_dat['wa'] = np.zeros(nlevels,dtype='float')

        ## This for loop reads the level data into the arrays in dict adf04_dat
        for i in range(0,nlevels):
            tmp_dat =f.readline()
            adf04_dat['config_data'][i] = tmp_dat[6:29]
            adf04_dat['isa'][i] = tmp_dat[41]
            adf04_dat['ila'][i] = tmp_dat[43]
            adf04_dat['xja'][i] = tmp_dat[46]
            adf04_dat['wa'][i] = tmp_dat[60:72]

        ## Skip a line
        f.readline()

        ## The next line in the adf04 file stores some valuable info
        tmp_dat = f.readline()
        adf04_dat['whats_this'] = float(tmp_dat[2:5]) # I don't know what this flag is for
        adf04_dat['filetype']  = int(tmp_dat[8:11])  # I believe this flag indicates col strength / cross section etc...
            # The adf04 designers, in an effort to save one byte of storage, thought it was wise
            # to store the temps in a single string, stripped of the character 'e'
            # The next few lines pull out the temperature string, split it into individual temps
            # and reinsert the 'e'. This was trickier than I thought it would be in python
        temps = str.split(tmp_dat[12:])
        tmp1=[x[0:3] for x in temps]
#        print("tmp1= ", tmp1)
        tmp2 = ["e" + x[4:7] for x in temps]
#        print("tmp2= ", tmp2)
        tmp_data=[x+y for x,y in zip(tmp1,tmp2)]
        adf04_dat['temps']=np.array([float(x) for x in tmp_data])
        del(temps)

        ## We need to create arrays to store the c-section/col-strength data
        adf04_dat['upper'] = np.zeros(ntrans, dtype='int')
        adf04_dat['lower'] = np.zeros(ntrans, dtype='int')
        adf04_dat['aval'] = np.zeros(ntrans, dtype='float')
        adf04_dat['gamma'] = np.zeros([ntrans, num_ergs], dtype = 'float')
        adf04_dat['inf_erg'] = np.zeros(ntrans, dtype = 'float')


        for i in range(0,ntrans):
            tmp_dat = f.readline()
            adf04_dat['upper'][i] = int(tmp_dat[0:4])
            adf04_dat['lower'][i] = int(tmp_dat[4:8])
            # pull the aval and reinsert the 'e'
            availa = tmp_dat[9:16]
#            print("availa=", availa, availa[0:4], availa[4:7])
            adf04_dat['aval'][i] = float(availa[0:4] + "e" + availa[4:7])
            # do the same for all the temps
            tmp_csec = str.split(tmp_dat[17:])
            tmp1 = [x[0:3] for x in tmp_csec]
            tmp2 = ["e" + x[4:7] for x in tmp_csec]
            tmp_csec = [x+y for x,y in zip(tmp1,tmp2)]
            tmp_csec = np.array([float(x) for x in tmp_csec])
            if np.size(tmp_csec) == num_ergs+1:
                    adf04_dat['inf_erg'][i] = tmp_csec[np.size(tmp_csec)-8]
                    adf04_dat['gamma'][i,:np.size(tmp_csec)-1] = tmp_csec[:np.size(tmp_csec)-1]
            elif np.size(tmp_csec) != num_ergs:
                    print('ERROR!!!!  The numbers of temps is not the same as the number of data points \
                                    in the cross section array!!!!!')
            else:
                    adf04_dat['gamma'][i:] = tmp_csec

        class Struct(object):
                def __init__(self, **entries):
                        self.__dict__.update(entries)

        adf04dat = Struct(**adf04_dat)

        return (adf04dat)

def branch(istart, ifinish, imeta, ntemps):

    from branch import csec
    CSEC_DATA = csec(istart, ifinish, imeta, ntemps)
    CSEC_DATA = CSEC_DATA.transpose()

    return (CSEC_DATA)
