import numpy as np
import matplotlib.pyplot as plt

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/'

#Nlatt, bta, ene, dev_ene, mag, dev_mag, sus, dev_sus, heatc, dev_heatc = np.loadtxt(directory+'output/data.txt', unpack=True)
mag1 = np.loadtxt(directory+'magnetization.txt', unpack=True)
tau = np.linspace(0, 9999, 10000)

plt.figure(1)
plt.title('Magnetic susceptibility, Nlatt=%d' % (40))
plt.ylabel(r'$\chi$')
plt.xlabel(r'$\beta$')
plt.grid(color = 'gray')
#plt.errorbar(bta,sus, yerr=dev_sus, color='green',fmt='.')
plt.plot(tau,mag1,color='green',marker=',',linewidth=1)

plt.show()
