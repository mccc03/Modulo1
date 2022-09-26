## Import libraries

import numpy as np
import matplotlib.pyplot as plt


## Set input data directory

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/'


## Load input data

Nlatt, bta, ene, dev_ene, mag, dev_mag, sus, dev_sus, heatc, dev_heatc = np.loadtxt(directory+'output/data.txt', unpack=True)


## Split data arrays into subarrays, one for each value of Nlatt, creating lists of lists, so that the order of the file data.txt is unimportant when plotting, allowing for additional simulations to be run

bta_split=[]
ene_split=[]
dev_ene_split=[]
mag_split=[]
dev_mag_split=[]
sus_split=[]
dev_sus_split=[]
heatc_split=[]
dev_heatc_split=[]
Nlatt_split=[]

for _ in range(6):
    bta_split.append([])
    ene_split.append([])
    dev_ene_split.append([])
    mag_split.append([])
    dev_mag_split.append([])
    sus_split.append([])
    dev_sus_split.append([])
    heatc_split.append([])
    dev_heatc_split.append([])

for i in range(len(Nlatt)):
    for j in range(6):
        if(Nlatt[i]==(10*j+10)):
            bta_split[j].append(bta[i])
            ene_split[j].append(ene[i])
            dev_ene_split[j].append(dev_ene[i])
            mag_split[j].append(mag[i])
            dev_mag_split[j].append(dev_mag[i])
            sus_split[j].append(sus[i])
            dev_sus_split[j].append(dev_sus[i])
            heatc_split[j].append(heatc[i])
            dev_heatc_split[j].append(dev_heatc[i])

for k in range(6):
    Nlatt_split.append(k*10 + 10)


## Plot magnetization for different Nlatt values

plt.figure(1)

plt.title('Magnetization')
plt.ylabel(r'$M$')
plt.xlabel(r'$\beta$')
plt.grid(color = 'gray')
plt.errorbar(bta_split[0],mag_split[0], yerr=dev_mag_split[0], color='gray',fmt=',',label="Nlatt=%d"%(Nlatt_split[0]))
plt.errorbar(bta_split[1],mag_split[1], yerr=dev_mag_split[1], color='orange',fmt='.',label="Nlatt=%d"%(Nlatt_split[1]))
plt.errorbar(bta_split[2],mag_split[2], yerr=dev_mag_split[2], color='red',fmt='<',label="Nlatt=%d"%(Nlatt_split[2]))
plt.errorbar(bta_split[3],mag_split[3], yerr=dev_mag_split[3], color='cyan',fmt='>',label="Nlatt=%d"%(Nlatt_split[3]))
plt.errorbar(bta_split[4],mag_split[4], yerr=dev_mag_split[4], color='green',fmt='v',label="Nlatt=%d"%(Nlatt_split[4]))
plt.errorbar(bta_split[5],mag_split[5], yerr=dev_mag_split[5], color='blue',fmt='^',label="Nlatt=%d"%(Nlatt_split[5]))
plt.legend(loc="upper left")


## Plot energy for different Nlatt values

plt.figure(2)

plt.title('Energy density')
plt.ylabel(r'$\epsilon$')
plt.xlabel(r'$\beta$')
plt.grid(color = 'gray')
plt.errorbar(bta_split[0],ene_split[0], yerr=dev_ene_split[0], color='gray',fmt=',',label="Nlatt=%d"%(Nlatt_split[0]))
plt.errorbar(bta_split[1],ene_split[1], yerr=dev_ene_split[1], color='orange',fmt='.',label="Nlatt=%d"%(Nlatt_split[1]))
plt.errorbar(bta_split[2],ene_split[2], yerr=dev_ene_split[2], color='red',fmt='<',label="Nlatt=%d"%(Nlatt_split[2]))
plt.errorbar(bta_split[3],ene_split[3], yerr=dev_ene_split[3], color='cyan',fmt='>',label="Nlatt=%d"%(Nlatt_split[3]))
plt.errorbar(bta_split[4],ene_split[4], yerr=dev_ene_split[4], color='green',fmt='v',label="Nlatt=%d"%(Nlatt_split[4]))
plt.errorbar(bta_split[5],ene_split[5], yerr=dev_ene_split[5], color='blue',fmt='^',label="Nlatt=%d"%(Nlatt_split[5]))
plt.legend(loc="upper right")


## Plot magnetic susceptibility for different Nlatt values

plt.figure(3)

plt.title('Magnetic susceptibility')
plt.ylabel(r'$\chi$')
plt.xlabel(r'$\beta$')
plt.grid(color = 'gray')
plt.errorbar(bta_split[0],sus_split[0], yerr=dev_sus_split[0], color='gray',fmt=',',label="Nlatt=%d"%(Nlatt_split[0]))
plt.errorbar(bta_split[1],sus_split[1], yerr=dev_sus_split[1], color='orange',fmt='.',label="Nlatt=%d"%(Nlatt_split[1]))
plt.errorbar(bta_split[2],sus_split[2], yerr=dev_sus_split[2], color='red',fmt='<',label="Nlatt=%d"%(Nlatt_split[2]))
plt.errorbar(bta_split[3],sus_split[3], yerr=dev_sus_split[3], color='cyan',fmt='>',label="Nlatt=%d"%(Nlatt_split[3]))
plt.errorbar(bta_split[4],sus_split[4], yerr=dev_sus_split[4], color='green',fmt='v',label="Nlatt=%d"%(Nlatt_split[4]))
plt.errorbar(bta_split[5],sus_split[5], yerr=dev_sus_split[5], color='blue',fmt='^',label="Nlatt=%d"%(Nlatt_split[5]))
plt.legend(loc="upper left")


## Plot heat capacity for different Nlatt values

plt.figure(4)

plt.title('Heat capacity')
plt.ylabel(r'$C$')
plt.xlabel(r'$\beta$')
plt.grid(color = 'gray')
plt.errorbar(bta_split[0],heatc_split[0], yerr=dev_heatc_split[0], color='gray',fmt=',',label="Nlatt=%d"%(Nlatt_split[0]))
plt.errorbar(bta_split[1],heatc_split[1], yerr=dev_heatc_split[1], color='orange',fmt='.',label="Nlatt=%d"%(Nlatt_split[1]))
plt.errorbar(bta_split[2],heatc_split[2], yerr=dev_heatc_split[2], color='red',fmt='<',label="Nlatt=%d"%(Nlatt_split[2]))
plt.errorbar(bta_split[3],heatc_split[3], yerr=dev_heatc_split[3], color='cyan',fmt='>',label="Nlatt=%d"%(Nlatt_split[3]))
plt.errorbar(bta_split[4],heatc_split[4], yerr=dev_heatc_split[4], color='green',fmt='v',label="Nlatt=%d"%(Nlatt_split[4]))
plt.errorbar(bta_split[5],heatc_split[5], yerr=dev_heatc_split[5], color='blue',fmt='^',label="Nlatt=%d"%(Nlatt_split[5]))
plt.legend(loc="upper left")


plt.show()
