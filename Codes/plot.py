## Import libraries

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


## Define functions for curve fits

def FitPower(x,g,m,c):
    return m*x**(-g) + c

def FitParabola(x,m,b,a):
    return m - b*(x-a)*(x-a)

def FitLine(x,m,c):
    return m*x + c

def FitBta(L, b_max, x_max, a):
    return b_max - x_max*(L**(-a))

def Fit(fit_func,x,y,dy,initial_values):
    arr_x = np.array(x)
    arr_y = np.array(y)
    arr_dy = np.array(dy)

    popt, pcov = curve_fit(fit_func, arr_x, arr_y, sigma=arr_dy, p0=initial_values, absolute_sigma=False)
    chisq = (((arr_y - fit_func(arr_x, *popt)) / arr_dy)**2).sum()

    return popt, pcov, chisq

## Define function to convert data into logarithmic scale

def Logalize(x_array, y_array, dev_y_array):
    x_log_list = []
    y_log_list = []
    dev_y_log_list = []

    for i in range(len(x_array)):
        x_log_list.append(np.log(x_array[i]))
        y_log_list.append(np.log(y_array[i]))
        dev_y_log_list.append(dev_y_array[i]/y_array[i])

    x_log = np.array(x_log_list)
    y_log = np.array(y_log_list)
    dev_y_log = np.array(dev_y_log_list)

    return x_log, y_log, dev_y_log

## Define function to find maxima and associated errors and store them into lists

def FindMaxMulti(x_matrix, y_matrix, dev_y_matrix):
    obs_max = []
    dev_obs_max = []
    x_max = []
    dev_x_max =[]

    for i in range(len(y_matrix)):
        tmp = 0.0
        dev_tmp = 0.0
        x_tmp = 0.0
        for j in range(len(y_matrix[i])):
            if(y_matrix[i][j]>tmp):
                tmp = y_matrix[i][j]
                dev_tmp = dev_y_matrix[i][j]
                x_tmp = x_matrix[i][j]
        obs_max.append(tmp)
        dev_obs_max.append(dev_tmp)
        x_max.append(x_tmp)
        dev_x_max.append(2*x_matrix[0][1]-2*x_matrix[0][0])

    return x_max, dev_x_max, obs_max, dev_obs_max

def FindMax(x_array, y_array, dev_y_array):
    tmp = 0.0
    dev_tmp = 0.0
    x_tmp = 0.0
    dev_x_tmp = x_array[1] - x_array[0]
    for j in range(len(y_array)):
        if(y_array[j]>tmp):
            tmp = y_array[j]
            dev_tmp = dev_y_array[j]
            x_tmp = x_array[j]
    return x_tmp, dev_x_tmp, tmp, dev_tmp



## Define function that creates lists of measures after critical point, considering bta_c = 0.44

def CreateFitArray(x_matrix, y_matrix, dev_y_matrix,size):
    t_list = []
    y_list = []
    dy_list = []

    x_max, dev_x_max, y_max, dev_y_max = FindMaxMulti(x_matrix, y_matrix, dev_y_matrix)

    i=0

    while(i<(len(x_matrix[size]))):
        if(x_matrix[size][i] > x_max[size] + 0.019):

            t_list.append(1-x_max[size]/x_matrix[size][i])
            y_list.append(y_matrix[size][i])
            dy_list.append(dev_y_matrix[size][i])

            if(len(t_list)==30):
                t = np.array(t_list)
                y = np.array(y_list)
                dy = np.array(dy_list)

                return t, y, dy
        i = i+1

## Set input data directory

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/'


## Load input data for plots and maxima fit

Nlatt, bta, ene, dev_ene, mag, dev_mag, sus, dev_sus, heatc, dev_heatc, binder, dev_binder = np.loadtxt(directory+'output/data.txt', unpack=True)

## Load input data for power laws fit
Nlatt1, bta1, mag1, dev_mag1, sus1, dev_sus1, heatc1, dev_heatc1 = np.loadtxt(directory+'output/data1.txt', unpack=True)

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

## Maxima analysis

bta_max, dev_bta_max, chi_max, dev_chi_max = FindMaxMulti(bta_split, sus_split, dev_sus_split)


## Start fit to find bta_c
initial_values = (0.44, 0.33, 1.0)
popt_bta, pcov_bta, chisq_bta = Fit(FitBta, Nlatt_split, bta_max, dev_bta_max, initial_values)

## Print fit values for bta_c

print('Fit values for beta_c studying maxima points in susceptibility\n')

print(f'bta_c = {popt_bta[0]:.4f} +/- {np.sqrt(pcov_bta[0, 0]):.4f}')
print(f'max_fss = {popt_bta[1]:.4f} +/- {np.sqrt(pcov_bta[1, 1]):.4f}')
print(f'nu = {popt_bta[2]:.4f} +/- {np.sqrt(pcov_bta[2, 2]):.4f}')
print(f'chisq = {chisq_bta:.4f}\n')

## Start fit to find gamma/nu

initial_values = (-1.75, 0.1, 1.0)
popt_gammanu, pcov_gammanu, chisq_gammanu = Fit(FitPower, Nlatt_split, chi_max, dev_chi_max, initial_values)

## Print fit values of gamma/nu

print('Fit values for gamma/nu studying maxima points in susceptibility\n')

print(f'gamma/nu = {-popt_gammanu[0]:.4f} +/- {np.sqrt(pcov_gammanu[0, 0]):.4f}')
print(f'chisq = {chisq_gammanu:.4f}')


## Using bta_c = 0.44, find gamma and beta exponents using fit in power laws

## Find gamma from susceptibility

t, chi, dev_chi = CreateFitArray(bta_split, sus_split, dev_sus_split, 5)

initial_values=(-1.75, 0.003, 0.0)
popt_gamma, pcov_gamma, chisq_gamma = Fit(FitPower, t, chi, dev_chi, initial_values)

## Print fit values of gamma

print('Fit values for gamma/nu studying maxima points in susceptibility\n')

print(f'gamma = {popt_gamma[0]:.4f} +/- {np.sqrt(pcov_gamma[0, 0]):.4f}')
print(f'chisq = {chisq_gamma:.4f}')


## Finite Size Scaling

## Create list of lists for fss variable

x_fss = []
for _ in range(len(Nlatt_split)):
    x_fss.append([])

## Fill fss variable, using bta_c = 0.44 and nu = 1

for i in range(len(bta_split)):
    for j in range(len(bta_split[i])):
        x_fss[i].append(Nlatt_split[i]*(bta_split[i][j] - 0.44))

## Create fss data

mag_fss = []
dev_mag_fss = []
sus_fss = []
dev_sus_fss = []
heatc_fss = []
dev_heatc_fss = []
binder_fss = []
dev_binder_fss = []

for _ in range(len(Nlatt_split)):
    mag_fss.append([])
    dev_mag_fss.append([])
    sus_fss.append([])
    dev_sus_fss.append([])
    heatc_fss.append([])
    dev_heatc_fss.append([])
    binder_fss.append([])
    dev_binder_fss.append([])

for i in range(len(mag_split)):
    for j in range(len(mag_split[i])):
        mag_fss[i].append((Nlatt_split[i]**(0.125))*mag_split[i][j])
        dev_mag_fss[i].append((Nlatt_split[i]**(0.125))*dev_mag_split[i][j])
        sus_fss[i].append((Nlatt_split[i]**(-1.75))*sus_split[i][j])
        dev_sus_fss[i].append((Nlatt_split[i]**(-1.75))*dev_sus_split[i][j])
        heatc_fss[i].append((Nlatt_split[i]**(-1.75))*heatc_split[i][j])
        dev_heatc_fss[i].append((Nlatt_split[i]**(-1.75))*dev_heatc_split[i][j])
        binder_fss[i].append(binder_split[i][j])
        dev_binder_fss[i].append(dev_binder_split[i][j])


## Plot magnetization for different Nlatt values

plt.figure(1)

plt.title('Magnetization')
plt.ylabel(r'$\langle | M |\rangle$')
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
plt.ylabel(r'$\langle\epsilon\rangle$')
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

## Plot FSS magnetization

plt.figure(6)

plt.title('Susceptibility - FSS')
plt.ylabel(r'$\chi N_{latt} ^{-\gamma / \nu}$')
plt.xlabel(r'$N_{latt} (\beta -\beta _c)$')
plt.grid(color = 'gray')
plt.errorbar(x_fss[0],sus_fss[0], yerr=dev_sus_fss[0], color='gray',fmt=',',label="Nlatt=%d"%(Nlatt_split[0]))
plt.errorbar(x_fss[1],sus_fss[1], yerr=dev_sus_fss[1], color='orange',fmt='.',label="Nlatt=%d"%(Nlatt_split[1]))
plt.errorbar(x_fss[2],sus_fss[2], yerr=dev_sus_fss[2], color='red',fmt='<',label="Nlatt=%d"%(Nlatt_split[2]))
plt.errorbar(x_fss[3],sus_fss[3], yerr=dev_sus_fss[3], color='cyan',fmt='>',label="Nlatt=%d"%(Nlatt_split[3]))
plt.errorbar(x_fss[4],sus_fss[4], yerr=dev_sus_fss[4], color='green',fmt='v',label="Nlatt=%d"%(Nlatt_split[4]))
plt.errorbar(x_fss[5],sus_fss[5], yerr=dev_sus_fss[5], color='blue',fmt='^',label="Nlatt=%d"%(Nlatt_split[5]))
plt.legend(loc="upper left")

plt.figure(7)

plt.title('Magnetization - FSS')
plt.ylabel(r'$M N_{latt} ^{\beta / \nu}$')
plt.xlabel(r'$N_{latt} (\beta -\beta _c)$')
plt.grid(color = 'gray')
plt.errorbar(x_fss[0],mag_fss[0], yerr=dev_mag_fss[0], color='gray',fmt=',',label="Nlatt=%d"%(Nlatt_split[0]))
plt.errorbar(x_fss[1],mag_fss[1], yerr=dev_mag_fss[1], color='orange',fmt='.',label="Nlatt=%d"%(Nlatt_split[1]))
plt.errorbar(x_fss[2],mag_fss[2], yerr=dev_mag_fss[2], color='red',fmt='<',label="Nlatt=%d"%(Nlatt_split[2]))
plt.errorbar(x_fss[3],mag_fss[3], yerr=dev_mag_fss[3], color='cyan',fmt='>',label="Nlatt=%d"%(Nlatt_split[3]))
plt.errorbar(x_fss[4],mag_fss[4], yerr=dev_mag_fss[4], color='green',fmt='v',label="Nlatt=%d"%(Nlatt_split[4]))
plt.errorbar(x_fss[5],mag_fss[5], yerr=dev_mag_fss[5], color='blue',fmt='^',label="Nlatt=%d"%(Nlatt_split[5]))
plt.legend(loc="upper left")

plt.figure(8)

plt.title('Binder - FSS')
plt.ylabel(r'$U$')
plt.xlabel(r'$N_{latt} (\beta -\beta _c)$')
plt.grid(color = 'gray')
plt.errorbar(x_fss[0],binder_fss[0], yerr=dev_binder_fss[0], color='gray',fmt=',',label="Nlatt=%d"%(Nlatt_split[0]))
plt.errorbar(x_fss[1],binder_fss[1], yerr=dev_binder_fss[1], color='orange',fmt='.',label="Nlatt=%d"%(Nlatt_split[1]))
plt.errorbar(x_fss[2],binder_fss[2], yerr=dev_binder_fss[2], color='red',fmt='<',label="Nlatt=%d"%(Nlatt_split[2]))
plt.errorbar(x_fss[3],binder_fss[3], yerr=dev_binder_fss[3], color='cyan',fmt='>',label="Nlatt=%d"%(Nlatt_split[3]))
plt.errorbar(x_fss[4],binder_fss[4], yerr=dev_binder_fss[4], color='green',fmt='v',label="Nlatt=%d"%(Nlatt_split[4]))
plt.errorbar(x_fss[5],binder_fss[5], yerr=dev_binder_fss[5], color='blue',fmt='^',label="Nlatt=%d"%(Nlatt_split[5]))
plt.legend(loc="upper left")


plt.figure(9)
plt.grid(color = 'gray')
plt.errorbar(t,chi,yerr=dev_chi)
plt.plot(t,FitPower(t,*popt_gamma))

plt.show()
