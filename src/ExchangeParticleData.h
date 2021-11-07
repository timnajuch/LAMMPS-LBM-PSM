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
//#include "fix_property_atom.h"

#include "lattice2D.h"
#include "unit_conversion.h"


using namespace std;

class ExchangeParticleData {
  public:
    ExchangeParticleData();
    ~ExchangeParticleData();

    void setParticlesOnLattice(Lattice2D *lattice2D_, Unit_Conversion *unitConversion, int numberParticles, LAMMPS_NS::tagint *tag, double **xPart, double **uPart, double **omega, double *rp, vector<double> boxLength, vector<double> origin);
    double calcSolidFraction(int i, int j, intk, double xP_LB, double yP_LB, double zP_LB, double rP_LB);
    void calculateHydrodynamicInteractions(Lattice2D *lattice2D_, Unit_Conversion *unitConversion, double *xPart, double rp, vector<double> &fHydro, vector<double> &tHydro, vector<double> &stresslet);

};

#endif
