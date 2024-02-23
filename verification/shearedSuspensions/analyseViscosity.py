import math
import numpy as np
import os

phi = np.linspace(0.05,0.25,5)
tau = 0.7
Re = 0.1
nu = 0.001
rho = 500.0
lc = 0.0005
d = lc

phiMax = phi[-1]
edimMin = 10.0*d
Vp = 4.0/3.0*math.pi*pow(d/2.0,3.0)
VboxMin = edimMin**3
NP = int(phiMax*VboxMin/Vp + 0.5)

Ndisc = 11.0
dx = d/(Ndisc-1.0)

stress_xy = []


for indexPhi in range(len(phi)):

    V_p = (4.0/3.0)*math.pi*pow(d/2.0,3.0)*NP
    V_domain = V_p/phi[indexPhi]
    edim = pow(V_domain,1.0/3.0) # End size of domain for create input script
    Nnodes = int(edim/dx+1.5) # Make sure that the length of edim results in an integer when divided by the lattice cell size
    edim = float(Nnodes-1)*dx
    edim = edim/2.0
    shearRate = nu/(d*d)*Re

    caseCreate = 'phi0p'+str(math.fmod(phi[indexPhi],1.0))[2:4]
    filename = caseCreate+'/'+caseCreate+'.stresses'
    file_exists = os.path.isfile(filename)
    if(file_exists):
        with open (filename, 'r') as datafile:
            lineID = 0
            for line in datafile:
                if(lineID > 1):
                    line = line.strip()
                    columns = line.split()
                    stress_xy.append(float(columns[2]))
                lineID += 1

    visc = stress_xy[-1]/shearRate
    mu = nu*rho   # dynamic viscosity
    visc_dimless = abs(visc/mu) + 1


    resultsFileName = './phi-viscosity.dat'
    file_exists = os.path.isfile(resultsFileName)

    with open (resultsFileName, 'a') as datafile:
        header = "# 1: Phi | 2: Dimensionless viscosity eta/eta_f]\n"

        if not file_exists:
            datafile.write(header)

        datafile.write( str(phi[indexPhi]) + " " + str(abs(visc_dimless)) + "\n" )

    datafile.close()
