/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#ifndef LBM_PSM_EXC_PART_DATA_H
#define LBM_PSM_EXC_PART_DATA_H

#include <algorithm>
#include <iostream>
#include <vector>

#include "atom.h"
#include "lmptype.h"

#include "LBM_PSM_lattice.h"
#include "LBM_PSM_unitConversion.h"


using namespace std;

class ExchangeParticleData {
  private: 
    int dimension;
    vector<double> origin;

  public:
    ExchangeParticleData(int dimension_, vector<double> origin_);
    ~ExchangeParticleData();

    void setParticlesOnLattice(LBMPSMLattice *lattice_, UnitConversion *unitConversion, int numberParticles, LAMMPS_NS::tagint *tag, double **xPart, double **uPart, double **omega, double *rp);
    double calcSolidFraction(int i, int j, int k, double xP_LB, double yP_LB, double zP_LB, double rP_LB);
    void calculateHydrodynamicInteractions(LBMPSMLattice *lattice_, UnitConversion *unitConversion, LAMMPS_NS::tagint tag, double *xPart, double rp, vector<double> &fHydro, vector<double> &tHydro, vector<double> &stresslet);

};

#endif
