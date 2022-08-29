import numpy as np
import matplotlib.pyplot as plt

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/'

Nlatt, bta, ene, var_ene, mag, var_mag, sus, var_sus, heatc, var_heatc = np.loadtxt(directory+'resampling_output_data.txt', skiprows=300, unpack=True)

plt.figure(1)
plt.title('Magnetic susceptibility, Nlatt=10')
plt.ylabel(r'$\chi$')
plt.xlabel(r'$\beta$')
plt.grid(color = 'gray')
plt.errorbar(bta,sus,yerr=var_sus, color='green',fmt='.')

plt.show()
