/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#ifndef PARTDATALATTICE_H
#define PARTDATALATTICE_H

#include <vector>
#include <iostream>

using namespace std;

class ParticleDataOnLattice{
  public:
    ParticleDataOnLattice();
    ~ParticleDataOnLattice();

    vector<int> particleID;
    vector<double> solidFraction;
    vector<double> particleVelocity;
    vector<double> hydrodynamicForce;

};


#endif
