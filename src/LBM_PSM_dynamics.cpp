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

          set_f0(i, j, k, iq, feq(iq, ind_phys_1D, ind_phys_2D) );
          set_f(i, j, k, iq, 0, feq(iq, ind_phys_1D, ind_phys_2D) );
          set_f(i, j, k, iq, 1, feq(iq, ind_phys_1D, ind_phys_2D) );
        }
      }
    }
  }
}


LBMPSMDynamics::~LBMPSMDynamics(){};


double LBMPSMDynamics::feq(int iq_, int ind_phys_1D_, int ind_phys_2D_){
  return rho[ind_phys_1D_] * w[iq_] * 
        (1.0 + ( e[3*iq_] * u[ind_phys_2D_] + e[3*iq_+1] * u[ind_phys_2D_+1] + e[3*iq_+2] * u[ind_phys_2D_+2]) / (pow(cs, 2.0))
        + pow( e[3*iq_] * u[ind_phys_2D_] + e[3*iq_+1] * u[ind_phys_2D_+1] + e[3*iq_+2] * u[ind_phys_2D_+2] , 2.0) / (2.0*pow(cs, 4.0))
        - 1.0/2.0 * ( u[ind_phys_2D_] * u[ind_phys_2D_] + u[ind_phys_2D_+1] * u[ind_phys_2D_+1] + u[ind_phys_2D_+2] * u[ind_phys_2D_+2]) / (pow(cs,2.0) ) );
}


double LBMPSMDynamics::feq(int iq_, int ind_phys_1D_, int ind_phys_2D_, vector<double> &rho_, vector<double> &u_){
  return rho_[ind_phys_1D_] * w[iq_] * 
        (1.0 + ( e[3*iq_] * u_[ind_phys_2D_] + e[3*iq_+1] * u_[ind_phys_2D_+1] + e[3*iq_+2] * u_[ind_phys_2D_+2]) / (pow(cs, 2.0))
        + pow( e[3*iq_] * u_[ind_phys_2D_] + e[3*iq_+1] * u_[ind_phys_2D_+1] + e[3*iq_+2] * u_[ind_phys_2D_+2] , 2.0) / (2.0*pow(cs, 4.0))
        - 1.0/2.0 * ( u_[ind_phys_2D_] * u_[ind_phys_2D_] + u_[ind_phys_2D_+1] * u_[ind_phys_2D_+1] + u_[ind_phys_2D_+2] * u_[ind_phys_2D_+2]) / (pow(cs,2.0) ) );
}


double LBMPSMDynamics::feq(int iq_, double rho, vector<double> u){
  return rho * w[iq_] * 
        (1.0 + ( e[3*iq_] * u[0] + e[3*iq_+1] * u[1] + e[3*iq_+2] * u[2]) / (pow(cs, 2.0))
        + pow( e[3*iq_] * u[0] + e[3*iq_+1] * u[1] + e[iq_*3+2] * u[2], 2.0) / (2.0*pow(cs, 4.0))
        - 1.0/2.0 * ( u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) / (pow(cs,2.0) ) );
}


double LBMPSMDynamics::F_iq(int iq_, vector<double> u, vector<double> F){
  return w[iq_] * 
        ( ( (e[3*iq_] - u[0]) * F[0] + (e[3*iq_+1] - u[1]) * F[1] + (e[3*iq_+2] - u[2]) * F[2]) / pow(cs, 2.0)
          + ( (e[3*iq_] * u[0] + e[3*iq_+1] * u[1] + e[iq_*3+2] * u[2]) / pow(cs, 4.0)
            * (e[3*iq_] * F[0] + e[3*iq_+1] * F[1] + e[3*iq_+2] * F[2]) ) );
}


double LBMPSMDynamics::F_iq(int iq_, int ind_phys_2D_, vector<double>& u, vector<double> F){
  return w[iq_] * 
        ( ( (e[3*iq_] - u[ind_phys_2D_+0]) * F[0] + (e[3*iq_+1] - u[ind_phys_2D_+1]) * F[1] + (e[3*iq_+2] - u[ind_phys_2D_+2]) * F[2]) / pow(cs, 2.0)
          + ( (e[3*iq_] * u[ind_phys_2D_+0] + e[3*iq_+1] * u[ind_phys_2D_+1] + e[iq_*3+2] * u[ind_phys_2D_+2]) / pow(cs, 4.0)
            * (e[3*iq_] * F[0] + e[3*iq_+1] * F[1] + e[3*iq_+2] * F[2]) ) );
}


double LBMPSMDynamics::F_iq(int iq_, int ind_phys_2D_, vector<double> F){
  return w[iq_] * 
        ( ( (e[3*iq_] - u[ind_phys_2D_+0]) * F[0] + (e[3*iq_+1] - u[ind_phys_2D_+1]) * F[1] + (e[3*iq_+2] - u[ind_phys_2D_+2]) * F[2]) / pow(cs, 2.0)
          + ( (e[3*iq_] * u[ind_phys_2D_+0] + e[3*iq_+1] * u[ind_phys_2D_+1] + e[iq_*3+2] * u[ind_phys_2D_+2]) / pow(cs, 4.0)
            * (e[3*iq_] * F[0] + e[3*iq_+1] * F[1] + e[3*iq_+2] * F[2]) ) );
}