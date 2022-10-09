## Import libraries

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


## Define functions for curve fits

def FitPower(x,g,m,c):
    return m*x**(-g) + c

def FitParabola(x,m,b,a):
    return m - b*(x-a)*(x-a)

def FitLinear(x,m,c):
    return m*x + c

def FitBta(L, b_max, x_max, a):
    return b_max - x_max*(L**(-a))

## This funtions returns the values for the best fit parameters as a list, for the covariance matrix as a matrix and for the normalized \chi^2 of the fit.

def Fit(fit_func,x,y,dy,initial_values):
    arr_x = np.array(x)
    arr_y = np.array(y)
    arr_dy = np.array(dy)

    popt, pcov = curve_fit(fit_func, arr_x, arr_y, sigma=arr_dy, p0=initial_values, absolute_sigma=False)

    chisq = (((arr_y - fit_func(arr_x, *popt)) / arr_dy)**2).sum()

    ndof = len(arr_y)-len(popt)
    chisq = chisq/ndof

    return popt, pcov, chisq

## Define function to convert list/array of linear data into logarithmic scale as numpy arrays

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

## Define function that creates lists of measures after critical point, considering bta_c = 0.44
## Variable size corresponds to Nlatt by Nlatt = size*10 + 10

def CreateFitArrayPower(x_matrix, y_matrix, dev_y_matrix, size, beta_c, l):
    t_list = []
    y_list = []
    dy_list = []

    i=0

    while(i<(len(x_matrix[size]))):
        if(x_matrix[size][i] > beta_c + 0.016): # Distance from beta_c is chosen arbitrarily, knowing that the data can not be too close (finite size effects) nor too far (phase transition effects not visible) from the critical point

            t_list.append(1-beta_c/x_matrix[size][i]) # Reduced temperature
            y_list.append(y_matrix[size][i])
            dy_list.append(dev_y_matrix[size][i])

            if(len(t_list)==l):
                t,y,dy = Logalize(t_list,y_list,dy_list)
                return t, y, dy
        i = i+1




## Set input data directory

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/'

## Set fit results file

output_fit = open(directory+'output/fit_results.txt',"w")

## Load input data for plots and fits

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

dist = 0 # Array length for each lattice size
for i in range(len(Nlatt)-1):
    if(Nlatt[i]<Nlatt[i+1]):
        Nlatt_split.append(Nlatt[i])
        if(dist==0):
            dist = i+1
Nlatt_split.append(Nlatt[len(Nlatt)-1])


## Create lists of lists, second index

for _ in range(len(Nlatt_split)):
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
    bta_split[i//dist].append(bta[i])
    ene_split[i//dist].append(ene[i])
    dev_ene_split[i//dist].append(dev_ene[i])
    mag_split[i//dist].append(mag[i])
    dev_mag_split[i//dist].append(dev_mag[i])
    sus_split[i//dist].append(sus[i])
    dev_sus_split[i//dist].append(dev_sus[i])
    heatc_split[i//dist].append(heatc[i])
    dev_heatc_split[i//dist].append(dev_heatc[i])
    binder_split[i//dist].append(binder[i])
    dev_binder_split[i//dist].append(dev_binder[i])




## Maxima analysis



## Start fit to find bta_c, maximum of phi_chi (not needed) and nu

bta_max, dev_bta_max, chi_max, dev_chi_max = FindMaxMulti(bta_split, sus_split, dev_sus_split)

initial_values = (0.44, 0.33, 1.0)
popt_bta, pcov_bta, chisq_bta = Fit(FitBta, Nlatt_split, bta_max, dev_bta_max, initial_values)

## Write fit values for bta_c, nu onto output file

output_fit.write('Fit values for beta_c, nu studying maxima points in susceptibility\n')

output_fit.write('bta_c\tdev_bta_c\tnu\tdev_nu\txbar\tdev_xbar\tchisq_norm\n')

output_fit.write(str(popt_bta[0])+'\t'+str(np.sqrt(pcov_bta[0, 0]))+'\t'+str(popt_bta[2])+'\t'+str(np.sqrt(pcov_bta[2, 2]))+'\t'+str(popt_bta[1])+'\t'+str(np.sqrt(pcov_bta[1, 1]))+'\t'+str(chisq_bta)+'\n')



## Start fit to find gamma/nu

# Nlatt_log, chi_max_log, dev_chi_max_log = Logalize(Nlatt_split, chi_max, dev_chi_max)

initial_values = (-1.75, 0.1, 0.0)
popt_gammanu, pcov_gammanu, chisq_gammanu = Fit(FitPower, Nlatt_split, chi_max, dev_chi_max, initial_values)

## Write fit values of gamma/nu onto output file

output_fit.write('Fit values for gamma/nu studying maxima points in susceptibility\n')

output_fit.write('gamma/nu\tdev_gamma/nu\tchisq_norm\n')

output_fit.write(str(-popt_gammanu[0])+'\t'+str(np.sqrt(pcov_gammanu[0, 0]))+'\t'+str(chisq_gammanu)+'\n')



## Using bta_c = 0.44, find gamma and beta exponents using fit in power laws

## Find gamma from susceptibility

t, chi, dev_chi = CreateFitArrayPower(bta_split, sus_split, dev_sus_split, 5, 0.4406868, 35)

initial_values=(-1.75, 0.003)
popt_gamma, pcov_gamma, chisq_gamma = Fit(FitLinear, t, chi, dev_chi, initial_values)

## Write fit values for gamma onto output file

output_fit.write('Fit values for gamma from power law fit\n')

output_fit.write('gamma\tdev_gamma\tchisq_norm\n')

output_fit.write(str(-popt_gamma[0])+'\t'+str(np.sqrt(pcov_gamma[0, 0]))+'\t'+str(chisq_gamma)+'\n')


## Find beta from magnetization

tb, m, dev_m = CreateFitArrayPower(bta_split, mag_split, dev_mag_split, 5, 0.4406868, 15)

initial_values=(0.125, 0.003)
popt_b, pcov_b, chisq_b = Fit(FitLinear, tb, m, dev_m, initial_values)

## Write fit values for beta onto output file

output_fit.write('Fit values for beta from power law fit\n')

output_fit.write('beta\tdev_beta\tchisq_norm\n')

output_fit.write(str(popt_b[0])+'\t'+str(np.sqrt(pcov_b[0, 0]))+'\t'+str(chisq_gamma)+'\n')


output_fit.close()




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




## Plots

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
plt.plot(t,FitLinear(t,*popt_gamma))

plt.figure(10)
plt.grid(color = 'gray')
plt.errorbar(tb,m,yerr=dev_m)
plt.plot(tb,FitLinear(tb,*popt_b))

plt.show()
