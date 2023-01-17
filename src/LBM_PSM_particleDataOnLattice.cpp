/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#include "LBM_PSM_particleDataOnLattice.h"


ParticleDataOnLattice::ParticleDataOnLattice()
{
  int numberParticles = 2;
  particleID.resize(numberParticles);
  solidFraction.resize(numberParticles);
  particleVelocity.resize(numberParticles*3);
  hydrodynamicForce.resize(numberParticles*3);

  particleID[0] = 0;
  particleID[1] = 0;
  solidFraction[0] = 0.0;
  solidFraction[1] = 0.0;
  particleVelocity[0] = 0.0;
  particleVelocity[1] = 0.0;
  particleVelocity[2] = 0.0;
  particleVelocity[3] = 0.0;
  particleVelocity[4] = 0.0;
  particleVelocity[5] = 0.0;
  hydrodynamicForce[0] = 0.0;
  hydrodynamicForce[1] = 0.0;
  hydrodynamicForce[2] = 0.0;
  hydrodynamicForce[3] = 0.0;
  hydrodynamicForce[4] = 0.0;
  hydrodynamicForce[5] = 0.0;
}

ParticleDataOnLattice::~ParticleDataOnLattice() {}
