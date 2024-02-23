#!/usr/bin/python

import numpy as np
import math
import subprocess 
import os, sys
import shutil

print ("Create particle configurations")

phi = np.linspace(0.05,0.25,5)

NP = 500.0
random = 24000
d = 0.0005
rho = 500.0

phiMax = phi[-1]
edimMin = 10.0*d
Vp = 4.0/3.0*math.pi*pow(d/2.0,3.0)
VboxMin = edimMin**3
NP = int(phiMax*VboxMin/Vp + 0.5)
print("Number of particles in system: {}".format(NP))

deff = d

skin = 1.5*d

kn = 1e5
kt = 2.0/7.0*kn
gamman = 10
gammat = 0.0
xmur = 0.0

m = 4.0/3.0*pow(d/2.0,3.0)*rho
meff = m
#dt_crit = math.pi/math.sqrt(kn/meff - pow(gamman/2.0,2.0)) 
dt_crit = math.pi/math.sqrt(kn/meff) 
dt = dt_crit/125
dt_create = dt_crit/200

N_timesteps_creation = 10000000
N_timesteps_relaxation = 50000000

N_particle_dump = 10000

bdim = d*50.0 # Beginning size of domain to create suspensions with above specified solid fraction
bdim = bdim/2

pathToLAMMPS = '~/lammps-stable_29Sep2021_update3/src/lmp_mpi'
processors = 8
nodes = 1
time = '24:00:00'
partition = 'short'


for indexPhi in range(len(phi)):
    print ('Current solid fraction:', phi[indexPhi])

    V_p = (4.0/3.0)*math.pi*pow(d/2.0,3.0)*NP
    V_domain = V_p/phi[indexPhi]
    edim = pow(V_domain,1.0/3.0) # End size of domain for create input script
    edim = edim/2-d

    caseCreate = 'phi0p'+str(math.fmod(phi[indexPhi],1.0))[2:4]
    simulationFolder = caseCreate
    print('Creating initial assembly: ', caseCreate)

    if( os.path.exists(simulationFolder) == False ):
        os.mkdir(simulationFolder)
    
    shutil.copy('create.lammps', simulationFolder)
    
    jobSubmissionFile = simulationFolder+'/job_submission_'+str(caseCreate)+'.sh'
 
    parseStrings = [
    '#!/bin/bash',
    '#SBATCH -o output.out',
    '#SBATCH -e output.err',
    '#SBATCH --job-name='+str(caseCreate),
    '#SBATCH --nodes='+str(nodes),
    '#SBATCH --time='+time,
    '#SBATCH --ntasks-per-node='+str(processors),
    '#SBATCH --partition='+partition,                        
    'srun --mpi=pmix '+pathToLAMMPS+\
    " -var NP "+str(int(NP))+" -var random "+str(int(random))+" -var d "+str(d)+" -var rho "+str(rho)+ \
    " -var kn "+str(kn)+" -var kt "+str(kt)+" -var gamman "+str(gamman)+" -var gammat "+str(gammat)+" -var xmur "+str(xmur)+" -var dt_create "+str(dt_create)+" -var skin "+str(skin)+ \
    " -var N_timesteps_creation "+str(int(N_timesteps_creation))+" -var N_timesteps_relaxation "+str(int(N_timesteps_relaxation))+ \
    " -var bdim "+str(bdim)+" -var edim "+str(edim)+" -var caseCreate "+caseCreate+ \
    " < create.lammps > "+caseCreate+".log"
    ]

    try:
        f = open(jobSubmissionFile,'x')
    except:
        f = open(jobSubmissionFile,'w')
    for iLine in range(len(parseStrings)):
        f.write(parseStrings[iLine]+'\n')
    f.close()

####   Executing LAMMPS create.lammps
    os.chdir(simulationFolder)
    os.system("sbatch job_submission_"+str(caseCreate)+".sh")
    os.chdir("..")

print('Finished')
