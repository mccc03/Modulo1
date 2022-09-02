import numpy as np
import matplotlib.pyplot as plt

# We begin by setting the parameters for the pseudo-random number generator
m = 2**64
a = 6364136223846793005
c = 1442695040888963407

# Employing the seed
seed = 2

# Selecting the starting point
x = 2.0

# Setting the parameters of the probability distribution
x0 = 5.
sigma = 1.

# Defining the step
D = .1

# Setting the number of iterations
N = 1024*1024*128

# Number of resamplings
k = 512

# Defining the pseudo-random number generator
def rand(seed1):
    global seed
    seed = (a*seed1 + c)%m
    return float(seed)/float(m)

# Defining the normal probability distribution
def p(x):
    return float(np.e**(-(x-x0)**2/(2*sigma**2)))

# Mean function (just for tidyness)
def mean_func(array):
    return float(sum(array))/float(len(array))

# Setting in-out data
Directory = '/home/mccc/Documents/Physics/MetodiNumerici/Modulo1/_data/'
#out_file = open(Directory+'MetDistr.txt', "w+")
#out_bootstrap = open(Directory+'mean_values.txt', "w+")


# Applying the Hasting's algorithm and saving the data in MetDistr.txt
"""
for i in range(N):
    x1 = x+(2*(rand(seed))-1)*D
    u = (rand(seed))
    if u < p(x1)/p(x):
        x = x1
    out_file.write(str(i)+' '+str(x)+'\n')
out_file.close()
"""
x, array = np.loadtxt('/home/mccc/Documents/Physics/MetodiNumerici/Modulo1/_data/MetDistr1.txt', unpack=True)

print("Starting bootstrap algorithm...\n")

# Building the bootstrap algorithm
"""
for l in range(1,k+1,1):
    array_l = []
    for i in range(int(N/l)):
        j = int(N*rand(seed))
        for s in range(l):
            if(j+l>N):
                array_l.append(array[N-s-1])
            else:
                array_l.append(array[i+s])
    mean_l = mean_func(array_l)
    out_bootstrap.write(str(l+1)+' '+str(mean_l)+'\n')
out_bootstrap.close()
"""
# Calculating mean value and variance
#l, mean_values = np.loadtxt('/home/mccc/Documents/Physics/MetodiNumerici/Modulo1/_data/mean_values.txt', unpack=True)
#mean = mean_func(mean_values)
#variance = np.sqrt(float(np.sum((mean_values-mean)**2))/float(len(mean_values)))


#print('mean = %f, variance = %f' % (mean, variance))

# Plotting the data
plt.figure(1)
plt.grid(color='gray')
plt.xlabel('Number of steps'); plt.ylabel('x')
plt.plot(x,array)
#plt.figure(2)
#plt.xlabel('Resampling size'); plt.ylabel('Mean')
#plt.plot(l,mean_values)
plt.show()
