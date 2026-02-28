/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#include "LBM_PSM_dynamics.h"

LBMPSMDynamics::LBMPSMDynamics(int nx_, int ny_, int nz_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_) :
  LBMPSMLattice(nx_, ny_, nz_, decomposition_, procCoordinates_, origin_, boxLength_, dimension_)
{
  for(int i = 0; i < nx; ++i){
    for(int j = 0; j < ny; ++j){
      for(int k = 0; k < nz; ++k){
        for(int iq = 0; iq < q; ++iq){
          int ind_phys_1D = index_1D(i, j, k);
          int ind_phys_2D = index_2D(i, j, k, 0);

          double feq_val = feq(iq, ind_phys_1D, ind_phys_2D);
          int ind_iq = index_fi(i, j, k, iq);
          f0[ind_iq]     = feq_val;
          f_curr[ind_iq] = feq_val;
          f_next[ind_iq] = feq_val;
        }
      }
    }
  }
}


LBMPSMDynamics::~LBMPSMDynamics(){};