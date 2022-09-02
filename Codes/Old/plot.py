import numpy as np
import matplotlib.pyplot as plt

N = 100000

x, array = np.loadtxt('/home/mccc/Documents/Physics/MetodiNumerici/Modulo1/_data/MetDistr.txt', unpack=True)
mean = float(array.sum())/float(N)
variance = (np.sqrt(float(np.sum((array-mean)**2))/float(N-1))/float(N))
print('mean = %f, variance = %f' % (mean, variance))

plt.figure(1)
plt.grid(color='gray')
plt.xlabel('Number of steps'); plt.ylabel('x')
plt.plot(x,array)
plt.show()
