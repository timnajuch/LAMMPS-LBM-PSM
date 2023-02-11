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

torqueAnalytical = 8.0*math.pi*nu*rho*(lc/2.0)**3.0*sr*0.5

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


torque_simulation = []

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
                torque_simulation.append(float(columns[15]))
                blockID += 1
            lineID += 1

print("Difference Torque: {}, {}, {}".format(torque_simulation[-1]/torqueAnalytical, torque_simulation[-1], torqueAnalytical))

resultsFileName = './ErrorSummary.dat'
file_exists = os.path.isfile(resultsFileName)

with open (resultsFileName, 'a') as datafile:
    header = "# 1: Tau | 2: StressletError[%] | 3: TorqueError[%]\n"

    if not file_exists:
        datafile.write(header)

    datafile.write( str(tau) + " " + str(abs(abs(stresslet_simulation[-1])-stresslet)/stresslet*100.0) + " " + str(abs(abs(torque_simulation[-1])-abs(torqueAnalytical))/abs(torqueAnalytical)*100.0) + "\n" )

datafile.close()
