# Imports
import numpy as np
import pandas as pd

# We begin by setting the parameters for the pseudo-random number generator
m = 2**64
a = 6364136223846793005
c = 1442695040888963407

# Defining the pseudo-random number generator
def rand(seed1):
    global seed
    seed = (a*seed1 + c)%m
    return float(seed)/float(m)

# Employing the seed
seed = 2

# Defining the dimensions
Nlatt = 10
Vol = Nlatt*Nlatt

# Setting the input-output options
Directory = '/home/mccc/Documents/Physics/MetodiNumerici/Modulo1/_data/'
Lattice = open(Directory+'lattice.txt', 'r')
Parameters = open(Directory+'parameters.txt','r')
Data = open(Directory+'data.txt','w+')

iflag = int(Parameters.readline())
measures = int(Parameters.readline())
i_decorrel = int(Parameters.readline())
extfield = float(Parameters.readline())
beta = float(Parameters.readline())

# Move functions
npp=[]
nmm=[]

for i in range(Nlatt):
    npp.append(i+1)
    nmm.append(i-1)
npp[Nlatt-1]=0
nmm[0]=Nlatt-1

# Defining the Metropolis algorithm
def Update_Metropolis():
    global seed
    for ivol in range(Vol):
        i = int(rand(seed)*Nlatt)
        j = int(rand(seed)*Nlatt)
        jm=nmm[j];jp=npp[j];im=nmm[i];ip=npp[i]
        force = field[i][jm] + field[i][jp] + field[im][j] + field [ip][j]
        force = beta*(force+extfield)
        phi = field[i][j]
        p_rat = np.e**(-2*phi*force)
        x = rand(seed)
        if x<p_rat:
            field[i][j] = -phi

# Magnetization measurement function
def Magnetization():
    global xmagn
    xmagn = 0
    for i in range(Nlatt):
        for j in range(Nlatt):
            xmagn = xmagn+field[i][j]
    return xmagn/float(Vol)

# Energy measurement function
def Energy():
    global xene
    xene = 0
    for i in range(Nlatt):
        for j in range(Nlatt):
            i = int(rand(seed)*Nlatt)
            j = int(rand(seed)*Nlatt)
            jm=nmm[j];jp=npp[j];im=nmm[i];ip=npp[i]
            force = field[i][jm] + field[i][jp] + field[im][j] + field [ip][j]
            xene = xene - 0.5*force*field[i][j] - extfield*field[i][j]
    return xene/float(Vol)

# Main program
LatticeIni(iflag,seed)
Lattice.close()

for l in range(measures):
    for k in range(i_decorrel):
        Update_Metropolis()
    xene = Energy()
    xmagn = Magnetization()
    Data.write(str(l)+' '+str(xmagn)+' '+str(xene)+'\n')
Lattice = open(Directory+'lattice.txt', 'w+')
for i in range(Nlatt):
    for j in range(Nlatt):
        Lattice.write(str(field[i][j])+'\n')
Data.close()
Lattice.close()
Parameters.close()
