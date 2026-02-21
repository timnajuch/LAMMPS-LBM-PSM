/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#include "LBM_PSM_BGK_dynamics.h"

LBMPSMBGKDynamics::LBMPSMBGKDynamics(double tau_, int nx_, int ny_, int nz_, vector<double> F_lbm_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_) :
  LBMPSMDynamics(nx_, ny_, nz_, decomposition_, procCoordinates_, origin_, boxLength_, dimension_), tau(tau_), omega(1.0/tau_), F_lbm(F_lbm_) { }


LBMPSMBGKDynamics::~LBMPSMBGKDynamics(){}


void  LBMPSMBGKDynamics::initialise_dynamics(double rho_, double ux_, double uy_, double uz_){
  for(int i = 0; i < nx; ++i){
    for(int j = 0; j < ny; ++j){
      for(int k = 0; k < nz; ++k){
        int ind_phys_1D = index_1D(i, j, k);
        int ind_phys_2D = index_2D(i, j, k, 0);
        rho[ind_phys_1D] = rho_;
        u[ind_phys_2D] = ux_;
        u[ind_phys_2D+1] = uy_;
        u[ind_phys_2D+2] = uz_;

        for(int iq = 0; iq < q; ++iq){
            set_f0(i, j, k, iq, feq(iq, ind_phys_1D, ind_phys_2D) );
            set_f(i, j, k, iq, 0, feq(iq, ind_phys_1D, ind_phys_2D) );
            set_f(i, j, k, iq, 1, feq(iq, ind_phys_1D, ind_phys_2D) );
        }
      }
    }
  }
  F_lbm_mag_pow2 = F_lbm[0]*F_lbm[0] + F_lbm[1]*F_lbm[1] + F_lbm[2]*F_lbm[2];
}



void LBMPSMBGKDynamics::compute_macro_values(int i_, int j_, int k_, int currentStep_){
  int ind_phys_1D = index_1D(i_, j_, k_);
  int ind_phys_2D = index_2D(i_, j_, k_, 0);

  double rho_tmp = 0.0;
  double jx = 0.0;
  double jy = 0.0;
  double jz = 0.0;

  double fi = 0.0;
  
  int ind_iq = 0;
  for(int iq = 0; iq < q; ++iq){
    ind_iq = index_fi(i_, j_, k_, iq, currentStep_);
    fi = f[ind_iq];
    rho_tmp += fi;
    jx += fi * e[3*iq];
    jy += fi * e[3*iq+1];
    jz += fi * e[3*iq+2];
  }
  rho[ind_phys_1D] = rho_tmp;

  // External force according to Guo et al. (2002)
  jx += F_lbm[0]*0.5;
  jy += F_lbm[1]*0.5;
  jz += F_lbm[2]*0.5;

  u[ind_phys_2D] = jx/rho_tmp;
  u[ind_phys_2D+1] = jy/rho_tmp;
  u[ind_phys_2D+2] = jz/rho_tmp;
}



void LBMPSMBGKDynamics::macroCollideStream(){
  int iShift = 0;
  int jShift = 0;
  int kShift = 0;

  for(int i = 0; i < nx; ++i){
    for(int j = 0; j < ny; ++j){
      for(int k = 0; k < nz; ++k){

        compute_macro_values(i, j, k, currentStep);

        for(int iq = 0; iq < q; ++iq){

          iShift = i + e[3*iq];
          if (iShift < 0 || iShift > nx-1) { continue; }
          jShift = j + e[3*iq+1];
          if (jShift < 0 || jShift > ny-1) { continue; }
          kShift = k + e[3*iq+2];
          if (kShift < 0 || kShift > nz-1) { continue; }

          collisionAndStream(i, j, k, iq, iShift, jShift, kShift, currentStep, nextStep);
        }

      }
    }
  }

  currentStep = nextStep;
  nextStep = 1 - currentStep;
}