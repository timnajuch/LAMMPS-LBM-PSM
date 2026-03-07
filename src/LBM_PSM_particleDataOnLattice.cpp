/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch
------------------------------------------------------*/

#include "LBM_PSM_particleDataOnLattice.h"


ParticleDataOnLattice::ParticleDataOnLattice()
{
  particleID.fill(0);
  solidFraction.fill(0.0);
  particleVelocity.fill(0.0);
  hydrodynamicForce.fill(0.0);
}

ParticleDataOnLattice::~ParticleDataOnLattice() {}
