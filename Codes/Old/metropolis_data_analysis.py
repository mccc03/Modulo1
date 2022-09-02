import numpy as np
import matplotlib.pyplot as plt

# Mean function (just for tidyness)
def mean_func(array):
    return float(sum(array))/float(len(array))

x, array = np.loadtxt('/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/MetDistr.txt', unpack=True)

# Calculating mean value and variance
l, mean_values = np.loadtxt('/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/mean_values.txt', unpack=True)
mean = mean_func(mean_values)
variance = np.sqrt(float(np.sum((mean_values-mean)**2))/float(len(mean_values)))


print('mean = %f, variance = %f' % (mean, variance))

# Plotting the data
fig1, axs1 = plt.subplots(2,1)
axs1[0].set_xlabel('Number of steps')
axs1[0].set_ylabel('x')
axs1[0].grid(color='gray')
axs1[0].plot(x,array)
axs1[1].set_xlabel('Bins')
axs1[1].set_ylabel('Binned array')
axs1[1].grid(color='gray')
axs1[1].hist(array,bins=50)
fig2, axs2 = plt.subplots(2,1)
axs2[0].set_xlabel('Resampling size')
axs2[0].set_ylabel('Mean')
axs2[0].grid(color='gray')
axs2[0].plot(l,mean_values, 'co')
axs2[1].set_xlabel('Bins')
axs2[1].set_ylabel('Binned array')
axs2[1].grid(color='gray')
axs2[1].hist(mean_values,bins=50)

plt.show()
