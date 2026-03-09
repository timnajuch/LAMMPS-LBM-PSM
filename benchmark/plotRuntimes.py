import matplotlib.pyplot as plt
import re
import numpy as np

numberOfProcessors = [1, 2, 4, 8, 12, 18]

logFilename = 'log.lammps'

runTimeInSeconds = []
relativeSpeedUp = [1] # Related to first list item

idealSpeedUp = [1] # Linear scaling in runtime reduction

for i, case in enumerate(numberOfProcessors):
    pathLog = 'Torque3D_{}/{}'.format(case, logFilename)

    with open(pathLog, "r") as logFile:
        for line in logFile:
            if(re.search("Total wall time", line)):
                lineList = line.strip().split()
                timeSplit = lineList[-1].split(":")
                runTimeInSeconds.append(int(timeSplit[0])*3600 + int(timeSplit[1])*60 + int(timeSplit[2]))
                print("Run time for case {}: {}".format(case, runTimeInSeconds[-1]))

                if (i > 0):
                    relativeSpeedUp.append(float(runTimeInSeconds[0])/float(runTimeInSeconds[i]))
                    idealSpeedUp.append(case/numberOfProcessors[0])

with open("benchmark-results.dat", "w") as outFile:
    outFile.write("# 1: NumberCores | 2: Run-time [s]\n")
    for noCores, runTime in zip(numberOfProcessors, runTimeInSeconds):
        outFile.write("{} {}\n".format(noCores, runTime))


plt.plot(numberOfProcessors, relativeSpeedUp, 'rx',\
         numberOfProcessors, idealSpeedUp, linestyle='dashed')

plt.legend(['Simulation results', 'Ideal speed-up'])
plt.xlabel('Number of cores')
plt.ylabel('Relative speed-up')

plt.xticks(np.arange(min(numberOfProcessors), max(numberOfProcessors)+1, 1))
plt.yticks(np.arange(1, int(idealSpeedUp[-1])+1, 1))
plt.grid()

plt.savefig("benchmarkResult.png")
plt.show()
