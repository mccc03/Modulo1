## Import libraries

import numpy as np
import matplotlib.pyplot as plt


## Set input data directory

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/'


## Load input data

Nlatt, bta, ene, dev_ene, mag, dev_mag, sus, dev_sus, heatc, dev_heatc = np.loadtxt(directory+'output/data.txt', unpack=True)

## Split data arrays into subarrays, one for each value of Nlatt, creating lists of lists

# For each value of Nlatt, we have run 150 simulations for varying bta
bta_split = np.split(bta, len(bta)/150)
sus_split = np.split(sus, len(sus)/150)
dev_sus_split = np.split(dev_sus, len(dev_sus)/150)


## Plot magnetization for different Nlatt values

plt.figure(1)

plt.title('Magnetization, Nlatt= %d, %d, %d, %d, %d' % (Nlatt[0], Nlatt[150],Nlatt[300],Nlatt[450],Nlatt[600]))
plt.ylabel(r'$M$')
plt.xlabel(r'$\beta$')
plt.grid(color = 'gray')
plt.errorbar(bta_split[0],mag_split[0], yerr=dev_mag_split[0], color='gray',fmt='.')
plt.errorbar(bta_split[1],mag_split[1], yerr=dev_mag_split[1], color='orange',fmt='o')
plt.errorbar(bta_split[2],mag_split[2], yerr=dev_mag_split[2], color='red',fmt='<')
plt.errorbar(bta_split[3],mag_split[3], yerr=dev_mag_split[3], color='cyan',fmt='>')
plt.errorbar(bta_split[4],mag_split[4], yerr=dev_mag_split[4], color='green',fmt='v')


## Plot energy for different Nlatt values

plt.figure(2)

plt.title('Energy density, Nlatt= %d, %d, %d, %d, %d' % (Nlatt[0], Nlatt[150],Nlatt[300],Nlatt[450],Nlatt[600]))
plt.ylabel(r'$\epsilon$')
plt.xlabel(r'$\beta$')
plt.grid(color = 'gray')
plt.errorbar(bta_split[0],ene_split[0], yerr=dev_ene_split[0], color='gray',fmt='.')
plt.errorbar(bta_split[1],ene_split[1], yerr=dev_ene_split[1], color='orange',fmt='o')
plt.errorbar(bta_split[2],ene_split[2], yerr=dev_ene_split[2], color='red',fmt='<')
plt.errorbar(bta_split[3],ene_split[3], yerr=dev_ene_split[3], color='cyan',fmt='>')
plt.errorbar(bta_split[4],ene_split[4], yerr=dev_ene_split[4], color='green',fmt='v')


## Plot magnetic susceptibility for different Nlatt values

plt.figure(3)

plt.title('Magnetic susceptibility, Nlatt= %d, %d, %d, %d, %d' % (Nlatt[0], Nlatt[150],Nlatt[300],Nlatt[450],Nlatt[600]))
plt.ylabel(r'$\chi$')
plt.xlabel(r'$\beta$')
plt.grid(color = 'gray')
plt.errorbar(bta_split[0],sus_split[0], yerr=dev_sus_split[0], color='gray',fmt='.')
plt.errorbar(bta_split[1],sus_split[1], yerr=dev_sus_split[1], color='orange',fmt='o')
plt.errorbar(bta_split[2],sus_split[2], yerr=dev_sus_split[2], color='red',fmt='<')
plt.errorbar(bta_split[3],sus_split[3], yerr=dev_sus_split[3], color='cyan',fmt='>')
plt.errorbar(bta_split[4],sus_split[4], yerr=dev_sus_split[4], color='green',fmt='v')


## Plot heat capacity for different Nlatt values

plt.figure(4)

plt.title('Heat capacity, Nlatt= %d, %d, %d, %d, %d' % (Nlatt[0], Nlatt[150],Nlatt[300],Nlatt[450],Nlatt[600]))
plt.ylabel(r'$C$')
plt.xlabel(r'$\beta$')
plt.grid(color = 'gray')
plt.errorbar(bta_split[0],heatc_split[0], yerr=dev_heatc_split[0], color='gray',fmt='.')
plt.errorbar(bta_split[1],heatc_split[1], yerr=dev_heatc_split[1], color='orange',fmt='o')
plt.errorbar(bta_split[2],heatc_split[2], yerr=dev_heatc_split[2], color='red',fmt='<')
plt.errorbar(bta_split[3],heatc_split[3], yerr=dev_heatc_split[3], color='cyan',fmt='>')
plt.errorbar(bta_split[4],heatc_split[4], yerr=dev_heatc_split[4], color='green',fmt='v')


plt.show()
