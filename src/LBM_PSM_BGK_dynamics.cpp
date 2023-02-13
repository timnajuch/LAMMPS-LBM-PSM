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
  LBMPSMDynamics(nx_, ny_, nz_, decomposition_, procCoordinates_, origin_, boxLength_, dimension_), tau(tau_), F_lbm(F_lbm_) { }


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
  
  int ind_iq = 0;
  for(int iq = 0; iq < q; ++iq){
    ind_iq = index_fi(i_, j_, k_, iq, currentStep_);
    rho_tmp += get_f(ind_iq);
    jx += get_f(ind_iq) * e[3*iq];
    jy += get_f(ind_iq) * e[3*iq+1];
    jz += get_f(ind_iq) * e[3*iq+2];
  }
  rho[ind_phys_1D] = rho_tmp;

  // External force according to Guo et al. (2002)
  jx += F_lbm[0]/2.0;
  jy += F_lbm[1]/2.0;
  jz += F_lbm[2]/2.0;

  u[ind_phys_2D] = jx/rho_tmp;
  u[ind_phys_2D+1] = jy/rho_tmp;
  u[ind_phys_2D+2] = jz/rho_tmp;
}



void LBMPSMBGKDynamics::collisionAndStream(int i_, int j_, int k_, int iq_, int iShift_, int jShift_, int kShift_, int currentStep_, int nextStep_){
  int ind_iq = index_fi(i_, j_, k_, iq_, currentStep_); // Index for populations
  int ind_iq0 = index_fi(i_, j_, k_, iq_, 0);           // Index for equilibrium populations
  int ind_phys_1D = index_1D(i_, j_, k_);
  int ind_phys_2D = index_2D(i_, j_, k_, 0);

  set_f0(i_, j_, k_, iq_, feq(iq_, ind_phys_1D, ind_phys_2D) );

  LAMMPS_NS::tagint pID1 = getParticleDataOnLatticeNode(ind_phys_1D).particleID[0];
  LAMMPS_NS::tagint pID2 = getParticleDataOnLatticeNode(ind_phys_1D).particleID[1];
  double Btmp1 = 0.0;
  double Btmp2 = 0.0;
  vector<double> uSolid1{0.0, 0.0, 0.0};
  vector<double> uSolid2{0.0, 0.0, 0.0};

  if (pID1 > 0){
    Btmp1 = getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[0];
    uSolid1[0] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[0];
    uSolid1[1] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[1];
    uSolid1[2] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[2];
  }
  if (pID2 > 0){
    Btmp2 = getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[1];
    uSolid2[0] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[3];
    uSolid2[1] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[4];
    uSolid2[2] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[5];
  }

  double f0_solid1 = feq(iq_, get_rho(ind_phys_1D), uSolid1);
  double f0_solid2 = feq(iq_, get_rho(ind_phys_1D), uSolid2);

  double BGKcoll = ( get_f0(ind_iq0) - get_f(ind_iq) ) / tau;

  double solid_coll1 = f0_solid1 - get_f(ind_iq) + ( 1.0 - 1.0/tau) * (get_f(ind_iq) - get_f0(ind_iq0) );
  double solid_coll2 = f0_solid2 - get_f(ind_iq) + ( 1.0 - 1.0/tau) * (get_f(ind_iq) - get_f0(ind_iq0) );

  double Btot = Btmp1 + Btmp2;
  
  double B1 = Btmp1;
  double B2 = Btmp2;
  if(Btot > 1.0){
    B1 = Btmp1/Btot;
    B2 = Btmp2/Btot;
    Btot = 1.0;
  }
 
  double F_lbm_iq = 0.0;
  if (F_lbm_mag_pow2 > 0.0)
    { F_lbm_iq = (1.0-0.5/tau) * F_iq(iq_, ind_phys_2D, F_lbm); }

  set_f(iShift_, jShift_, kShift_, iq_, nextStep_, get_f(ind_iq)  + ( 1.0 - Btot ) * BGKcoll  + B1 * solid_coll1 + B2 * solid_coll2 + F_lbm_iq);

  add_Fhyd(ind_phys_1D, pID1, -B1 * solid_coll1 * e[3*iq_], 0);
  add_Fhyd(ind_phys_1D, pID1, -B1 * solid_coll1 * e[3*iq_+1], 1);
  add_Fhyd(ind_phys_1D, pID1, -B1 * solid_coll1 * e[3*iq_+2], 2);
  add_Fhyd(ind_phys_1D, pID2, -B2 * solid_coll2 * e[3*iq_], 0);
  add_Fhyd(ind_phys_1D, pID2, -B2 * solid_coll2 * e[3*iq_+1], 1);
  add_Fhyd(ind_phys_1D, pID2, -B2 * solid_coll2 * e[3*iq_+2], 2);
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