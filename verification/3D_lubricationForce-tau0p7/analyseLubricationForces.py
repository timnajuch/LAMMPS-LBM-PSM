import math
import numpy as np
import os
import matplotlib.pyplot as plt

## Input parameters

Re = 0.1
nu = 0.0001
rho = 1000.0
lc = 0.0005

Uc = Re*nu/lc

mu = nu*rho
r1 = lc/2.0
r2 = lc/2.0

N_discretisation = 11.0
dx = lc/N_discretisation

## Read simulation forces

force_simulation_1 = []
force_simulation_2 = []
pos_simulation_1 = []
pos_simulation_2 = []

filename = 's.p'
file_exists = os.path.isfile(filename)
if(file_exists):
    with open (filename, 'r') as datafile:
        lineID = 1
        blockID = 1
        for line in datafile:
            if(lineID == 10*blockID+1*(blockID-1) or lineID == 10*blockID+1*(blockID-1)+1):
                line = line.strip()
                columns = line.split()
                if(int(columns[0]) == 1):
                    force_simulation_1.append(float(columns[10]))
                    pos_simulation_1.append([float(columns[4]), float(columns[5]), float(columns[6])])
                elif(int(columns[0]) == 2):
                    force_simulation_2.append(float(columns[10]))
                    pos_simulation_2.append([float(columns[4]), float(columns[5]), float(columns[6])])
                if(lineID == 10*blockID+1*(blockID-1)+1):
                    blockID += 1
            lineID += 1

N = len(force_simulation_1)
ccdx = np.zeros(N)
ccdy = np.zeros(N)
ccd = np.zeros(N)
gap = np.zeros(N)
gap_dimensionless = np.zeros(N)
gap_latticeCells = np.zeros(N)

fn1 = np.asarray(force_simulation_1, dtype=float)
fn2 = np.asarray(force_simulation_2, dtype=float)
pos_simulation_1 = np.asarray(pos_simulation_1, dtype=float)
pos_simulation_2 = np.asarray(pos_simulation_2, dtype=float)


# Calculate dimensionless parameters / output
for i in range(N):
    ccdx[i] = pos_simulation_1[i,0] - pos_simulation_2[i,0]
    ccdy[i] = pos_simulation_1[i,1] - pos_simulation_2[i,1]
    ccd[i] = math.sqrt(ccdx[i]*ccdx[i] + ccdy[i]*ccdy[i])


gap[:] = ccd[:] - r1 - r2
gap_dimensionless[:] = gap[:]/r1 #/(r1+r2)*2.0
gap_latticeCells[:] = gap[:]/dx
fn1[:] /= (6.0*math.pi*r1*mu*Uc)
fn2[:] /= (6.0*math.pi*r2*mu*Uc)


## Compute the analytical solution

X_A_11 = np.zeros(N)
flubanalyticaln = np.zeros(N)
errorn = np.zeros(N)

beta0 = r2/r1
beta1 = 1.0+beta0

b0p2 = beta0*beta0
b1p2 = beta1*beta1
b1p3 = beta1*b1p2

for i in range(N-1):
    if(gap_dimensionless[i] > 0.0):
        # Jeffrey (1982)
        flubanalyticaln[i] = 0.5/gap_dimensionless[i] - 9.0/20.0*math.log(gap_dimensionless[i]) + 1.34558 - 3.0/56.0*gap_dimensionless[i]*math.log(gap_dimensionless[i]) + 0.19*gap_dimensionless[i]
        errorn[i] = abs((abs(flubanalyticaln[i]) - abs(fn1[i])))/abs(flubanalyticaln[i])*100.0



### PLOTS (uncomment and comment the plot you need) ###

pltIndex = int(0.0*N)

plt.subplot(211)
plt.loglog(gap_latticeCells[pltIndex:], abs(fn1[pltIndex:]), 'bs', gap_latticeCells[pltIndex:], flubanalyticaln[pltIndex:], 'k-')
plt.axvline(x = 1, color = 'r', linestyle = '--')
plt.xlim(0.1,10)
plt.xticks([0.1, 1, 10])
plt.xlabel('Number of lattice cells in gap')
plt.ylabel('Force/Force_{Stokes}')
plt.subplot(212)
plt.plot(gap_latticeCells[pltIndex:], errorn[pltIndex:], 'bs--')
plt.axhline(y = 5, color = 'r', linestyle = '--')
plt.axvline(x = 1, color = 'r', linestyle = '--')
plt.xlim(0.0,10)
plt.ylim(0.0,30)
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
plt.xlabel('Number of lattice cells in gap')
plt.ylabel('Error [%]')
plt.show()

