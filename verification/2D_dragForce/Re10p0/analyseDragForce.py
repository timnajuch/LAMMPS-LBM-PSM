import math
import numpy as np
import os

Re = 10.0
nu = 0.0001
rho = 1000.0
lc = 0.0005
Uc = Re*nu/lc

dx = lc/10.0
A = lc*dx

force_simulation = []

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
                force_simulation.append(float(columns[10]))
                blockID += 1
            lineID += 1


force_simulation = np.array(force_simulation, float)
meanFd = np.mean(force_simulation[int(len(force_simulation)*0.8):len(force_simulation)])
Cd = abs(meanFd) / (0.5*rho*pow(Uc, 2)*A)  # drag coefficient

resultsFileName = '../dragCoefficients.dat'
file_exists = os.path.isfile(resultsFileName)
with open (resultsFileName, 'a') as datafile:
    header = "# 1: Re | 2: DragCoefficient\n"

    if not file_exists:
        datafile.write(header)

    datafile.write( str(Re) + " " + str(Cd) + "\n" )

datafile.close()
