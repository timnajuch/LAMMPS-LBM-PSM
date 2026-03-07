/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch
------------------------------------------------------*/

#ifndef LBMPSM_PARTDATALATTICE_H
#define LBMPSM_PARTDATALATTICE_H

#include <array>

#include "lmptype.h"

using namespace std;

class ParticleDataOnLattice{
  public:
    ParticleDataOnLattice();
    ~ParticleDataOnLattice();

    array<LAMMPS_NS::tagint, 2> particleID{};
    array<double, 2> solidFraction{};
    array<double, 6> particleVelocity{};
    array<double, 6> hydrodynamicForce{};

};

#endif
