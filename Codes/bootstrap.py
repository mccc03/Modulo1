## Import libraries

import string
import numpy as np


## I/O directory and files

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/'
input_parameters=open(directory+'input/parameters.txt',"r")
input_bta=open(directory+'input/bta.txt',"r")
output_data=open(directory+'output/data1.txt',"a")


## Read input parameters

input_par_list = input_parameters.readlines()
input_bta_list = input_bta.readlines()

bta = float(input_bta_list[0])- 0.0005
resamplings = int(input_par_list[1])
Nlatt = int(input_par_list[3])
thermal_len = int(input_par_list[6])

input_parameters.close()
input_bta.close()


## Initialize RNG

rng = np.random.default_rng()


## Define functions

def susceptibility(array):
    return bta*Nlatt*Nlatt*np.var(array)

def heat_capacity(array):
    return bta*bta*Nlatt*Nlatt*np.var(array)

def mean2(array):
    array2=[]
    for i in range(len(array)):
        array2.append(array[i]*array[i])
    return np.mean(array2)

def mean4(array):
    array4=[]
    for i in range(len(array)):
        array4.append(array[i]*array[i]*array[i]*array[i])
    return np.mean(array4)

def bootstrap_dev(observable,array,resamplings):
    dev_list = []
    len_array = len(array)
    bin_size = 32
    while(bin_size <= 1+len_array/10):
        observable_bs = []
        for _ in range(resamplings):
            array_res = []
            for _ in range(int(len_array/bin_size)):
                i = rng.integers(len_array)
                array_res.extend(array[min(i,len_array-1-bin_size):min(i+bin_size,len_array-1)])
            observable_bs.append(observable(array_res))
        dev_list.append(np.std(observable_bs))
        bin_size*=2
    return max(dev_list)


## User interaction
print('Running bootsrap algorithm for bta = %f\n'%(bta))

## Load magnetization and energy arrays

mag_list = np.loadtxt(directory+'magnetization.txt', unpack=True,skiprows=thermal_len)
ene_list = np.loadtxt(directory+'energy.txt', unpack=True,skiprows=thermal_len)


## Compute observables and associated errors

#ene = np.mean(ene_list)
#ene_dev = bootstrap_dev(np.mean,ene_list,resamplings)
mag = np.mean(mag_list)
mag_dev = bootstrap_dev(np.mean,mag_list,resamplings)
sus = susceptibility(mag_list)
sus_dev = bootstrap_dev(susceptibility,mag_list,resamplings)
heatc = heat_capacity(ene_list)
heatc_dev = bootstrap_dev(heat_capacity,ene_list,resamplings)

#m2 = mean2(mag_list)
#m2_dev = bootstrap_dev(mean2,mag_list,resamplings)
#m4 = mean4(mag_list)
#m4_dev = bootstrap_dev(mean4,mag_list,resamplings)
#binder = 1-m4/(3*m2*m2)
#binder_dev = np.sqrt(((m4_dev*m4_dev)/((3*m2*m2)*(3*m2*m2)))+((2*m4)*(2*m4)*(m2_dev*m2_dev)/((3*m2*m2*m2)*(3*m2*m2*m2))))


## Write results onto output file

output_data.write(str(Nlatt)+'\t'+str(bta)+'\t'+str(mag)+'\t'+str(mag_dev)+'\t'+str(sus)+'\t'+str(sus_dev)+'\t'+str(heatc)+'\t'+str(heatc_dev)+'\t'+'\n')

output_data.close()
