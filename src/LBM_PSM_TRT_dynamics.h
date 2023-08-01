/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details. 

Tim Najuch, 2022
------------------------------------------------------*/

#ifndef LBMPSM_TRT_DYNAMICS
#define LBMPSM_TRT_DYNAMICS

#include "LBM_PSM_lattice.h"
#include "LBM_PSM_dynamics.h"
#include "domain.h"

class LBMPSMTRTDynamics : public LBMPSMDynamics {
  
  private:
    double tau_p, tau_m, omega_p, omega_m;      // Relaxation parameters tau_p, tau_m, and their inverse omega_p and omega_m
    vector<double> F_lbm;
    double F_lbm_mag_pow2;
    vector<double> F_hyd_1;
    vector<double> F_hyd_2;

  public:
    LBMPSMTRTDynamics(double tau_, double tau_m_, int nx_, int ny_, int nz_, vector<double> F_lbm_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_);
    ~LBMPSMTRTDynamics();

    void initialise_dynamics(double rho_, double ux_, double uy_, double uz_);
    void compute_macro_values(int i_, int j_, int k_, int currentStep_);
    void collisionAndStream(int i_, int j_, int k_, int iq_, int iShift_, int jShift_, int kShift_, int currentStep_, int nextStep_);
    void macroCollideStream();
};

#endif
