/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#include "LBM_PSM_BGK_dynamics.h"

LBMPSMBGKDynamics::LBMPSMBGKDynamics(double tau_, int nx_, int ny_, int nz_, int q_, vector<double> F_lbm_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_, double dx_) :
  LBMPSMDynamics(nx_, ny_, nz_, q_, decomposition_, procCoordinates_, origin_, boxLength_, dimension_, dx_), tau(tau_), F_lbm(F_lbm_) {};


LBMPSMBGKDynamics::~LBMPSMBGKDynamics(){};


void LBMPSMBGKDynamics::compute_macro_values(){
  for(int i = 0; i < LBMPSMLattice::nx; ++i){
    for(int j = 0; j < LBMPSMLattice::ny; ++j){
      for(int k = 0; k < LBMPSMLattice::nz; ++k){

        int ind_phys_1D = i * LBMPSMLattice::ny * LBMPSMLattice::nz + j * LBMPSMLattice::nz + k;
        int ind_phys_2D = i * LBMPSMLattice::ny * LBMPSMLattice::nz * 3 + j*LBMPSMLattice::nz*3 + k*3;

        double rho_tmp = 0.0;
        double jx = 0.0;
        double jy = 0.0;
        double jz = 0.0;
      
        for(int iq = 0; iq < LBMPSMLattice::q; ++iq){
          int ind_iq = i * LBMPSMLattice::ny * LBMPSMLattice::nz * LBMPSMLattice::q + j * LBMPSMLattice::nz*LBMPSMLattice::q + k * LBMPSMLattice::q + iq; 
          rho_tmp += LBMPSMLattice::get_f(ind_iq);
          jx += LBMPSMLattice::get_f(ind_iq) * LBMPSMLattice::e[3*iq];
          jy += LBMPSMLattice::get_f(ind_iq) * LBMPSMLattice::e[3*iq+1];
          jz += LBMPSMLattice::get_f(ind_iq) * LBMPSMLattice::e[3*iq+2];
        }
        LBMPSMLattice::rho[ind_phys_1D] = rho_tmp;

  //      jx += F_lbm[0]/2.0;
  //      jy += F_lbm[1]/2.0;
        LBMPSMLattice::u[ind_phys_2D] = jx/rho_tmp;
        LBMPSMLattice::u[ind_phys_2D+1] = jy/rho_tmp;
        LBMPSMLattice::u[ind_phys_2D+2] = jz/rho_tmp;

        //p[i][j] = rho_tmp*csPow2;
      }
    }
  }
};


void LBMPSMBGKDynamics::compute_macro_values(int i_, int j_, int k_){
  int ind_phys_1D = i_ * LBMPSMLattice::ny * LBMPSMLattice::nz + j_ * LBMPSMLattice::nz + k_;
  int ind_phys_2D = i_ * LBMPSMLattice::ny * LBMPSMLattice::nz * 3 + j_*LBMPSMLattice::nz*3 + k_*3;

  double rho_tmp = 0.0;
  double jx = 0.0;
  double jy = 0.0;
  double jz = 0.0;

  for(int iq = 0; iq < LBMPSMLattice::q; ++iq){
    int ind_iq = i_ * LBMPSMLattice::ny * LBMPSMLattice::nz * LBMPSMLattice::q + j_ * LBMPSMLattice::nz * LBMPSMLattice::q + k_*LBMPSMLattice::q + iq; 
    rho_tmp += LBMPSMLattice::get_f(ind_iq);
    jx += LBMPSMLattice::get_f(ind_iq) * LBMPSMLattice::e[3*iq];
    jy += LBMPSMLattice::get_f(ind_iq) * LBMPSMLattice::e[3*iq+1];
    jz += LBMPSMLattice::get_f(ind_iq) * LBMPSMLattice::e[3*iq+2];
  }
  LBMPSMLattice::rho[ind_phys_1D] = rho_tmp;

  LBMPSMLattice::u[ind_phys_2D] = jx/rho_tmp;
  LBMPSMLattice::u[ind_phys_2D+1] = jy/rho_tmp;
  LBMPSMLattice::u[ind_phys_2D+2] = jz/rho_tmp;
};


void LBMPSMBGKDynamics::collision(int i_, int j_, int k_, int iq_){
        
  int ind_iq = i_ * LBMPSMLattice::ny * LBMPSMLattice::nz * LBMPSMLattice::q + j_ * LBMPSMLattice::nz * LBMPSMLattice::q + k_*LBMPSMLattice::q + iq_;
  int ind_phys_1D = i_ * LBMPSMLattice::ny *LBMPSMLattice::nz + j_*LBMPSMLattice::nz + k_;
  int ind_phys_2D = i_ * LBMPSMLattice::ny *LBMPSMLattice::nz* 3 + j_*3*LBMPSMLattice::nz + k_*3;

  LBMPSMLattice::set_f0(i_, j_, k_, iq_, LBMPSMDynamics::feq(iq_, ind_phys_1D, ind_phys_2D, LBMPSMLattice::rho, LBMPSMLattice::u) );

  LAMMPS_NS::tagint pID1 = getParticleDataOnLatticeNode(ind_phys_1D).particleID[0];
  LAMMPS_NS::tagint pID2 = getParticleDataOnLatticeNode(ind_phys_1D).particleID[1];
  double Btmp1 = 0.0;
  double Btmp2 = 0.0;
  vector<double> uSolid1{0.0, 0.0, 0.0};
  vector<double> uSolid2{0.0, 0.0, 0.0};

  if (pID1 > 0){
    Btmp1 = LBMPSMLattice::getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[0];
    //double eps1 = LBMPSMLattice::getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[0];
    //Btmp1 = eps1*(tau-0.5)/(1.0-eps1+tau-0.5);
    uSolid1[0] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[0];
    uSolid1[1] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[1];
    uSolid1[2] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[2];
  }
  if (pID2 > 0){
    Btmp2 = LBMPSMLattice::getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[1];
    //double eps2 = LBMPSMLattice::getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[1];
    //Btmp2 = eps2*(tau-0.5)/(1.0-eps2+tau-0.5);
    uSolid2[0] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[3];
    uSolid2[1] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[4];
    uSolid2[2] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[5];
  }

  double f0_solid1 = LBMPSMDynamics::feq(iq_, LBMPSMLattice::get_rho(ind_phys_1D), uSolid1);
  double f0_solid2 = LBMPSMDynamics::feq(iq_, LBMPSMLattice::get_rho(ind_phys_1D), uSolid2);

  double BGKcoll = ( LBMPSMLattice::get_f0(ind_iq) - LBMPSMLattice::get_f(ind_iq) ) / tau;
/*
  // Eq8
  int iq_m = 0;
  if(iq_ > 0)
  {
    if(LBMPSMLattice::dimension == 2){
      if(iq_ < 5)
      {
        iq_m = iq_ < 3 ? iq_ + 2 : iq_ - 2;
      }
      else
      {
        iq_m = iq_ < 7 ? iq_ + 2 : iq_ - 2;
      }
    } else { //3D
      if (iq_ % 2 == 0){
        iq_m = iq_ - 1;
      } else {
        iq_m = iq_ + 1;
      }
    }
  }
  int ind_iq_m = i_ * LBMPSMLattice::ny * LBMPSMLattice::nz * LBMPSMLattice::q + j_ * LBMPSMLattice::nz * LBMPSMLattice::q + k_*LBMPSMLattice::q + iq_m;
  //LBMPSMLattice::set_f0(i_, j_, k_, iq_m, LBMPSMDynamics::feq(iq_m, ind_phys_1D, ind_phys_2D, LBMPSMLattice::rho, LBMPSMLattice::u) );

  //double solid_coll1 = f0_solid1 - LBMPSMLattice::get_f(ind_iq) + LBMPSMLattice::get_f(ind_iq_m) - LBMPSMLattice::get_f0(ind_iq_m); //Eq8
  //double solid_coll2 = f0_solid2 - LBMPSMLattice::get_f(ind_iq) + LBMPSMLattice::get_f(ind_iq_m) - LBMPSMLattice::get_f0(ind_iq_m); //Eq8
  vector<double> uF{LBMPSMLattice::get_u(ind_phys_2D), LBMPSMLattice::get_u(ind_phys_2D+1), LBMPSMLattice::get_u(ind_phys_2D+2)};
  double solid_coll1 = f0_solid1 - LBMPSMLattice::get_f(ind_iq) + LBMPSMLattice::get_f(ind_iq_m) - LBMPSMDynamics::feq(iq_m, LBMPSMLattice::get_rho(ind_phys_1D), uF); //Eq8
  double solid_coll2 = f0_solid2 - LBMPSMLattice::get_f(ind_iq) + LBMPSMLattice::get_f(ind_iq_m) - LBMPSMDynamics::feq(iq_m, LBMPSMLattice::get_rho(ind_phys_1D), uF); //Eq8
*/  
  // Eq9
  double solid_coll1 = f0_solid1 - LBMPSMLattice::get_f(ind_iq) + ( 1.0 - 1.0/tau) * (LBMPSMLattice::get_f(ind_iq) - LBMPSMLattice::get_f0(ind_iq) ); //Eq9
  double solid_coll2 = f0_solid2 - LBMPSMLattice::get_f(ind_iq) + ( 1.0 - 1.0/tau) * (LBMPSMLattice::get_f(ind_iq) - LBMPSMLattice::get_f0(ind_iq) ); //Eq9

  double Btot = Btmp1 + Btmp2;
  
  double B1 = Btmp1;
  double B2 = Btmp2;
  if(Btot > 1.0){
    B1 = Btmp1/Btot;
    B2 = Btmp2/Btot;
  }

  if(Btot > 1.0)
  { Btot = 1.0; }
    
  LBMPSMLattice::set_f(i_, j_, k_, iq_, LBMPSMLattice::get_f(ind_iq)  + ( 1.0 - Btot ) * BGKcoll  + B1 * solid_coll1 + B2 * solid_coll2 );

  if(iq_ == 0)
  {
    LBMPSMLattice::set_Fhydx(ind_phys_1D, 0, -B1 * solid_coll1 * LBMPSMLattice::e[3*iq_]);
    LBMPSMLattice::set_Fhydy(ind_phys_1D, 0, -B1 * solid_coll1 * LBMPSMLattice::e[3*iq_+1]);
    LBMPSMLattice::set_Fhydz(ind_phys_1D, 0, -B1 * solid_coll1 * LBMPSMLattice::e[3*iq_+2]);
    LBMPSMLattice::set_Fhydx(ind_phys_1D, 1, -B2 * solid_coll2 * LBMPSMLattice::e[3*iq_]);
    LBMPSMLattice::set_Fhydy(ind_phys_1D, 1, -B2 * solid_coll2 * LBMPSMLattice::e[3*iq_+1]);
    LBMPSMLattice::set_Fhydz(ind_phys_1D, 1, -B2 * solid_coll2 * LBMPSMLattice::e[3*iq_+2]);
  }else{
    LBMPSMLattice::add_Fhydx(ind_phys_1D, 0, -B1 * solid_coll1 * LBMPSMLattice::e[3*iq_]);
    LBMPSMLattice::add_Fhydy(ind_phys_1D, 0, -B1 * solid_coll1 * LBMPSMLattice::e[3*iq_+1]);
    LBMPSMLattice::add_Fhydz(ind_phys_1D, 0, -B1 * solid_coll1 * LBMPSMLattice::e[3*iq_+2]);
    LBMPSMLattice::add_Fhydx(ind_phys_1D, 1, -B2 * solid_coll2 * LBMPSMLattice::e[3*iq_]);
    LBMPSMLattice::add_Fhydy(ind_phys_1D, 1, -B2 * solid_coll2 * LBMPSMLattice::e[3*iq_+1]);
    LBMPSMLattice::add_Fhydz(ind_phys_1D, 1, -B2 * solid_coll2 * LBMPSMLattice::e[3*iq_+2]);
  }

};


void  LBMPSMBGKDynamics::initialise_dynamics(double rho_, double ux_, double uy_, double uz_){

  for(int i = 0; i < LBMPSMLattice::nx; ++i){
    for(int j = 0; j < LBMPSMLattice::ny; ++j){
      for(int k = 0; k < LBMPSMLattice::nz; ++k){
        int ind_phys_1D = i * LBMPSMLattice::ny * LBMPSMLattice::nz + j*LBMPSMLattice::nz + k;
        int ind_phys_2D = i * LBMPSMLattice::ny *LBMPSMLattice::nz * 3 + j*LBMPSMLattice::nz*3 + k*3;
        LBMPSMLattice::rho[ind_phys_1D] = rho_;
        LBMPSMLattice::u[ind_phys_2D] = ux_;
        LBMPSMLattice::u[ind_phys_2D+1] = uy_;
        LBMPSMLattice::u[ind_phys_2D+2] = uz_;

        for(int iq = 0; iq < LBMPSMLattice::q; ++iq){
            int ind_iq = i * LBMPSMLattice::ny * LBMPSMLattice::nz * LBMPSMLattice::q + j * LBMPSMLattice::nz *LBMPSMLattice::q + k*LBMPSMLattice::q + iq; 

            LBMPSMLattice::set_f0(i, j, k, iq, LBMPSMDynamics::feq(iq, ind_phys_1D, ind_phys_2D, LBMPSMLattice::rho, LBMPSMLattice::u) );
            LBMPSMLattice::set_f(i, j, k, iq, LBMPSMDynamics::feq(iq, ind_phys_1D, ind_phys_2D, LBMPSMLattice::rho, LBMPSMLattice::u) );
            LBMPSMLattice::set_fcoll(i, j, k, iq, LBMPSMDynamics::feq(iq, ind_phys_1D, ind_phys_2D, LBMPSMLattice::rho, LBMPSMLattice::u) );

        }
      }
    }
  }

};


void LBMPSMBGKDynamics::macroCollideStream(){
  for(int i = 0; i < LBMPSMLattice::nx; ++i){
    for(int j = 0; j < LBMPSMLattice::ny; ++j){
      for(int k = 0; k < LBMPSMLattice::nz; ++k){

        compute_macro_values(i, j, k);

        for(int iq = 0; iq < LBMPSMLattice::q; ++iq){
          collision(i, j, k, iq);

          if(LBMPSMLattice::dimension == 2){
            if( (i > 0) && (i < LBMPSMLattice::nx-1) && (j > 0) && (j < LBMPSMLattice::ny-1) ){
              LBMPSMDynamics::streamBulk(i, j, k, iq);
            }else{
              LBMPSMDynamics::streamBC(i,j,k,iq);
            }
          }else{
            if( (i > 0) && (i < LBMPSMLattice::nx-1) && (j > 0) && (j < LBMPSMLattice::ny-1) && (k > 0) && (k < LBMPSMLattice::nz-1) ){
              LBMPSMDynamics::streamBulk(i, j, k, iq);
            }else{
              LBMPSMDynamics::streamBC(i,j,k,iq);
            }
          }
        }

      }
    }
  }

  LBMPSMLattice::cp_fcoll_f();

};
