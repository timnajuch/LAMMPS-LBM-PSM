/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#ifndef EXC_PART_DATA_H
#define EXC_PART_DATA_H

#include <algorithm>
#include <iostream>
#include <vector>

#include "atom.h"
#include "lmptype.h"

#include "lattice2D.h"
#include "unit_conversion.h"


using namespace std;

class ExchangeParticleData {
  private:
    // Particle diameter (in LB system)
    double dp;  
    // Particle position (in LB system)
    std::vector<double> xp;    
    // Particle velocity (in LB system)
    std::vector<double> us;


  public:
    //ExchangeParticleData(double dp_, std::vector<double> xp_, std::vector<double> us_);
    ExchangeParticleData();
    ~ExchangeParticleData();

  // Add methods for placing particles on lattice and obtaining force
  // So far, only one particle can be set via constructor

    //void setParticlesOnLattice(Lattice2D &lattice2D_);
    void setParticlesOnLattice(Lattice2D *lattice2D_);
    void setParticlesOnLattice(Lattice2D *lattice2D_, Unit_Conversion *unitConversion, int numberParticles, LAMMPS_NS::tagint *tag, double **xPart, double **uPart, double *rp, vector<double> boxLength, vector<double> origin);
    double calcSolidFraction(int i, int j, double xP_LB, double yP_LB, double rP_LB);
    void calculateHydrodynamicInteractions(Lattice2D *lattice2D_, Unit_Conversion *unitConversion, int numberParticles, LAMMPS_NS::tagint *tag, double **xPart, double *rp, int hydroForceID, int hydroTorqueID, int stressletID, LAMMPS_NS::Atom *atom, vector<double> origin);

};

#endif
