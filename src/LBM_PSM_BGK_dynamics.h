/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details. 

Tim Najuch, 2022
------------------------------------------------------*/

#ifndef LBMPSM_BGK_DYNAMICS
#define LBMPSM_BGK_DYNAMICS

#include "LBM_PSM_lattice.h"
#include "LBM_PSM_dynamics.h"
#include "domain.h"

class LBMPSMBGKDynamics : public LBMPSMDynamics {
  
  private:
    double tau;
    vector<double> F_lbm;
    int dimension;

  public:
    LBMPSMBGKDynamics(double tau_, int nx_, int ny_, int nz_, int q_, vector<double> F_lbm_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_, double dx_);
    ~LBMPSMBGKDynamics();

    void compute_macro_values();

    void initialise_dynamics(double rho_, double ux_, double uy_, double uz_);

    void macroCollideStream();
    void compute_macro_values(int i_, int j_, int k_);
    void collision(int i_, int j_, int k_, int iq_);
    
};

#endif
