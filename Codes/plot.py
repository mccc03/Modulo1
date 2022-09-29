## Import libraries

import numpy as np
import matplotlib.pyplot as plt


## Set input data directory

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/'


## Load input data

Nlatt, bta, ene, dev_ene, mag, dev_mag, sus, dev_sus, heatc, dev_heatc, binder, dev_binder = np.loadtxt(directory+'output/data.txt', unpack=True)


## Split data arrays into subarrays, one for each value of Nlatt, creating lists of lists, so that the order of the file data.txt is unimportant when plotting, allowing for additional simulations to be run

## Create lists of lists, first index

Nlatt_split=[]
bta_split=[]
ene_split=[]
dev_ene_split=[]
mag_split=[]
dev_mag_split=[]
sus_split=[]
dev_sus_split=[]
heatc_split=[]
dev_heatc_split=[]
binder_split=[]
dev_binder_split=[]

## Create lists of lists, second index

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
    binder_split.append([])
    dev_binder_split.append([])

## Fill arrays from data

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
            binder_split[j].append(binder[i])
            dev_binder_split[j].append(dev_binder[i])

for k in range(6):
    Nlatt_split.append(k*10 + 10)


## Find list index of maximum of susceptibility
p=0
for l in range(len(sus_split[5])-1):
    if(sus_split[5][l+1]>sus_split[5][l]):
        p=l+1

## Approximate bta_c with point of maximum of susceptibility
bta_c = bta_split[5][p]

## Create reduced temperature variable
t = []
chi = []
dev_chi=[]
for m in range(15):
    t.append(np.log(bta_c*(1.0/bta_split[5][p-15+m])-1))
    chi.append(np.log(sus_split[5][p-15+m]))
    dev_chi.append(dev_sus_split[5][p-15+m]/sus_split[5][p-15+m])

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

## Plot Binder cumulant for different Nlatt values

plt.figure(5)

plt.title('Binder cumulant')
plt.ylabel(r'$1-\frac{\langle m^4\rangle}{3\langle m^2\rangle ^2}$')
plt.xlabel(r'$\beta$')
plt.grid(color = 'gray')
plt.errorbar(bta_split[0],binder_split[0], yerr=dev_binder_split[0], color='gray',fmt=',',label="Nlatt=%d"%(Nlatt_split[0]))
plt.errorbar(bta_split[1],binder_split[1], yerr=dev_binder_split[1], color='orange',fmt='.',label="Nlatt=%d"%(Nlatt_split[1]))
plt.errorbar(bta_split[2],binder_split[2], yerr=dev_binder_split[2], color='red',fmt='<',label="Nlatt=%d"%(Nlatt_split[2]))
plt.errorbar(bta_split[3],binder_split[3], yerr=dev_binder_split[3], color='cyan',fmt='>',label="Nlatt=%d"%(Nlatt_split[3]))
plt.errorbar(bta_split[4],binder_split[4], yerr=dev_binder_split[4], color='green',fmt='v',label="Nlatt=%d"%(Nlatt_split[4]))
plt.errorbar(bta_split[5],binder_split[5], yerr=dev_binder_split[5], color='blue',fmt='^',label="Nlatt=%d"%(Nlatt_split[5]))
plt.legend(loc="upper left")

## Plot chi vs t

plt.figure(6)

plt.title('Susceptibility near critical point')
plt.ylabel(r'$M$')
plt.xlabel(r'$\beta$')
plt.grid(color = 'gray')
plt.errorbar(t,chi, yerr=dev_chi, color='red')

plt.show()
