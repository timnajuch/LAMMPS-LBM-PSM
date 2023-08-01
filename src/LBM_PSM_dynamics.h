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
    LBMPSMDynamics(int nx_, int ny_, int nz_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_);
    ~LBMPSMDynamics();

    double feq(int iq_, int ind_phys_1D_, int ind_phys_2D_); // Accesses the density and velocities stored in the Lattice class
    double feq(int iq_, double rho_, vector<double> u_);

    // External force according to Guo et al. (2002). When called needs to be multiplied by (1-0.5/tau)
    double F_iq(int iq_, vector<double> u_, vector<double> F_);
    double F_iq(int iq_, int ind_phys_2D_, vector<double> F_); // Accesses the velocities stored in the Lattice class

    virtual void initialise_dynamics(double rho_, double ux_, double uy_, double uz_){};
    virtual void compute_macro_values(int i_, int j_, int k_, int currentStep_){};
    virtual void collisionAndStream(int i_, int j_, int k_, int iq_, int iShift_, int jShift_, int kShift_, int currentStep_, int nextStep_){};
    virtual void macroCollideStream(){};

};

#endif