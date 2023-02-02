/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
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
    LBMPSMDynamics(int nx_, int ny_, int nz_, int q_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_);
    ~LBMPSMDynamics();

    double feq(int iq_, int ind_phys_1D_, int ind_phys_2D_, vector<double> &rho_, vector<double> &u_);
    double feq(int iq_, double rho, vector<double> u);

    void streamBulk(int i_, int j_, int k_, int iq_);
    void streamBC(int i_, int j_, int k_, int iq_);

    // External force according to Guo et al. (2002). When called needs to be multiplied by (1-0.5/tau)
    double F_iq(int iq_, vector<double> u, vector<double> F);
    double F_iq(int iq_, int ind_phys_2D_, vector<double>& u, vector<double> F);

};

#endif

