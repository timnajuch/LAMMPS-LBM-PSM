/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#ifndef LBMPSM_PARTDATALATTICE_H
#define LBMPSM_PARTDATALATTICE_H

#include <iostream>
#include <vector>

#include "lmptype.h"

using namespace std;

class ParticleDataOnLattice{
  public:
    ParticleDataOnLattice();
    ~ParticleDataOnLattice();

    vector<LAMMPS_NS::tagint> particleID;
    vector<double> solidFraction;
    vector<double> particleVelocity;
    vector<double> hydrodynamicForce;
};

#endif
