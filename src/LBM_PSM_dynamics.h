/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#ifndef LBMPSMDYNAMICS_2D
#define LBMPSMDYNAMICS_2D

#include <iostream>
#include <math.h>
#include <vector>

#include "LBM_PSM_lattice.h"

using namespace std;

class LBMPSMDynamics : public LBMPSMLattice {
  
  public:
    LBMPSMDynamics(int nx_, int ny_, int nz_, int q_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_, double dx_);
    ~LBMPSMDynamics();

    double feq(int iq_, int ind_phys_1D_, int ind_phys_2D_, vector<double> &rho_, vector<double> &u_);
    double feq(int iq_, double rho, vector<double> u);

    void streamBulk(int i_, int j_, int k_, int iq_);
    void streamBC(int i_, int j_, int k_, int iq_);

};

#endif
