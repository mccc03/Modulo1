import numpy as np
import matplotlib.pyplot as plt
import random

def mean_func(array):
    return float(sum(array))/float(len(array))

N = 1000000
resamples = 500

directory='/home/exterior//Documents/Physics/MetodiNumerici/Modulo1/_data/'
ene_array = np.loadtxt(directory+'energy.txt')
out_file = open(directory+'energy_r.txt', "w+")
ene_array_res = [0]*N

for k in range(resamples):
    ene_array_res
    for i in range(len(ene_array)):
        j = random.randint(0, N-1)
        ene_array_res[i] = ene_array[j]
    ene_mean = mean_func(ene_array_res)
    out_file.write(str(k)+'\t'+str(ene_mean)+'\n')

out_file.close()
