import os
import shutil

print ("Run benchmark")

numberOfProcessors = [18, 12, 8, 4, 2, 1]
pathToLAMMPS = '~/LAMMPS/lammps/src/lmp_mpi'

for processors in numberOfProcessors:
    print ('Current number of cores', processors)

    simulationFolder = 'Torque3D_'+str(processors)
    print('Run case: ', simulationFolder)
    
    if os.path.isdir(simulationFolder):
        shutil.rmtree(simulationFolder)
    os.mkdir(simulationFolder)
    shutil.copy('in.lbmpsm_'+str(processors), simulationFolder)

    os.chdir(simulationFolder)
    os.system("mpirun -np {} {} < in.lbmpsm_{}".format(processors, pathToLAMMPS, processors))
    os.chdir("..")

print('Finished')
