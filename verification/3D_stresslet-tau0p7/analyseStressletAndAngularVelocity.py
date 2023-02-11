import math
import numpy as np
import os

tau = 0.7
Re = 0.1
nu = 0.0001
rho = 1000.0
lc = 0.0005

Uc = Re*nu/lc
h = 0.005
sr = Uc*2.0/h

Exy = 0.5*sr

stresslet = 20.0/3.0*math.pi*nu*rho*(lc/2.0)**3.0*Exy

omegaAnalytical = -Uc/h

stresslet_simulation = []

filename = 'stressOutputFIX'
file_exists = os.path.isfile(filename)
if(file_exists):
    with open (filename, 'r') as datafile:
        lineID = 0
        for line in datafile:
            if(lineID > 1):
                line = line.strip()
                columns = line.split()
                stresslet_simulation.append(float(columns[5])*h**3.0)
            lineID += 1

print("Difference: {}, {}, {}".format(stresslet_simulation[-1]/stresslet, stresslet_simulation[-1], stresslet))

omega_simulation = []

filename = 's.p'
file_exists = os.path.isfile(filename)
if(file_exists):
    with open (filename, 'r') as datafile:
        lineID = 1
        blockID = 1
        for line in datafile:
            if(lineID == 10*blockID):
                line = line.strip()
                columns = line.split()
                omega_simulation.append(float(columns[18]))
                blockID += 1
            lineID += 1

print("Difference Omega: {}, {}, {}".format(omega_simulation[-1]/omegaAnalytical, omega_simulation[-1], omegaAnalytical))

resultsFileName = './ErrorSummary.dat'
file_exists = os.path.isfile(resultsFileName)

with open (resultsFileName, 'a') as datafile:
    header = "# 1: Tau | 2: StressletError[%] | 3: RotationalVelocityError[%]\n"

    if not file_exists:
        datafile.write(header)

    datafile.write( str(tau) + " " + str(abs(abs(stresslet_simulation[-1])-stresslet)/stresslet*100.0) + " " + str(abs(abs(omega_simulation[-1])-abs(omegaAnalytical))/abs(omegaAnalytical)*100.0) + "\n" )

datafile.close()
