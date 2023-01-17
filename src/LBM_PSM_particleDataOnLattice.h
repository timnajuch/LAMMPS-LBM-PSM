/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
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
