import math
import numpy as np
import os

tau = 0.7
Re = 0.1
nu = 0.0001
rho = 1000.0
lc = 0.0005

Uc = Re*nu/lc

stokesForceAnalytical = 6.0*math.pi*nu*rho*lc/2.0*Uc


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

print("Difference in forces: {}, {}, {}".format(force_simulation[-1]/stokesForceAnalytical, force_simulation[-1], stokesForceAnalytical))

resultsFileName = './ForceError.dat'
file_exists = os.path.isfile(resultsFileName)

with open (resultsFileName, 'a') as datafile:
    header = "# 1: Tau | 2: StokesForceError[%]\n"

    if not file_exists:
        datafile.write(header)

    datafile.write( str(tau) + " " + str(abs(abs(force_simulation[-1])-stokesForceAnalytical)/stokesForceAnalytical*100.0) + "\n" )

datafile.close()
