/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#include "particleDataOnLattice.h"


ParticleDataOnLattice::ParticleDataOnLattice()
{
  int numberParticles = 2;
  particleID.resize(numberParticles);
  solidFraction.resize(numberParticles);
  particleVelocity.resize(numberParticles*3);
  hydrodynamicForce.resize(numberParticles*3);
};

ParticleDataOnLattice::~ParticleDataOnLattice() {};
