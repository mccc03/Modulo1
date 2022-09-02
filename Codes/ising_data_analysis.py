## Import libraries

import string
import numpy as np
import matplotlib.pyplot as plt


## Define functions

def mean_func(array):
    return float(sum(array))/float(len(array))

def variance(array, mean):
    return np.sqrt(np.sum((array-mean)**2)/len(array))


## I/O directory and files

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/'
input_parameters=open(directory+'input/parameters.txt',"r")
input_bta=open(directory+'input/bta.txt',"r")
output_data=open(directory+'output/resampling_output_data.txt',"a") # access mode for appending lines


## Read input parameters

input_par_list = input_parameters.readlines()
input_bta_list = input_bta.readlines()

bta = float(input_bta_list[0])- 0.001
Nlatt = int(input_par_list[4])


## Create numpy lists from input data

en_index, en = np.loadtxt(directory+'ene_res.txt', unpack=True)
mag_index, mag = np.loadtxt(directory+'mag_res.txt', unpack=True)
sus_index, sus = np.loadtxt(directory+'sus_res.txt', unpack=True)
heatc_index, heatc = np.loadtxt(directory+'heatc_res.txt', unpack=True)


## Compute mean and variance

mean_en = mean_func(en)
var_en = variance(en,mean_en)

mean_mag = mean_func(mag)
var_mag = variance(mag,mean_mag)

mean_sus = mean_func(sus)
var_sus = variance(sus,mean_sus)

mean_heatc = mean_func(heatc)
var_heatc = variance(heatc,mean_heatc)


## Print results

print('Simulation results for beta = %f , Nlatt = %d \n' % (bta, Nlatt))
print('energy = %f +- %f' % (mean_en, var_en))
print('magnetization = %f +- %f' % (mean_mag, var_mag))
print('susceptibility = %f +- %f' % (mean_sus, var_sus))
print('heat capacity =  %f +- %f \n' % (mean_heatc, var_heatc))


## Store results onto output files

output_data.write(str(Nlatt)+'\t'+str(bta)+'\t'+str(mean_en)+'\t'+str(var_en)+'\t'+str(mean_mag)+'\t'+str(var_mag)+'\t'+str(mean_sus)+'\t'+str(var_sus)+'\t'+str(mean_heatc)+'\t'+str(var_heatc)+'\n')
