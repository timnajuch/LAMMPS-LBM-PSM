/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#include "LBM_PSM_dynamics.h"

LBMPSMDynamics::LBMPSMDynamics(int nx_, int ny_, int nz_, int q_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_, double dx_) :
  LBMPSMLattice(nx_, ny_, nz_, q_, decomposition_, procCoordinates_, origin_, boxLength_, dimension_, dx_)
{
  for(int i = 0; i < LBMPSMLattice::nx; ++i){
    for(int j = 0; j < LBMPSMLattice::ny; ++j){
      for(int k = 0; k < LBMPSMLattice::nz; ++k){
        for(int iq = 0; iq < LBMPSMLattice::q; ++iq){
          int ind_phys_1D = i * LBMPSMLattice::ny * LBMPSMLattice::nz + j * LBMPSMLattice::nz + k;
          int ind_phys_2D = (i * LBMPSMLattice::ny * LBMPSMLattice::nz + j * LBMPSMLattice::nz + k)*3;
          int ind_iq = i * LBMPSMLattice::ny * LBMPSMLattice::nz * LBMPSMLattice::q + j * LBMPSMLattice::nz * LBMPSMLattice::q + k*LBMPSMLattice::q + iq;

          LBMPSMLattice::set_f0(i, j, k, iq, LBMPSMDynamics::feq(iq, ind_phys_1D, ind_phys_2D, LBMPSMLattice::rho, LBMPSMLattice::u) );
          LBMPSMLattice::set_f(i, j, k, iq, LBMPSMDynamics::feq(iq, ind_phys_1D, ind_phys_2D, LBMPSMLattice::rho, LBMPSMLattice::u) );
          LBMPSMLattice::set_fcoll(i, j, k, iq, LBMPSMDynamics::feq(iq, ind_phys_1D, ind_phys_2D, LBMPSMLattice::rho, LBMPSMLattice::u) );
        }
      }
    }
  }
};


LBMPSMDynamics::~LBMPSMDynamics(){};
    

double LBMPSMDynamics::feq(int iq_, int ind_phys_1D_, int ind_phys_2D_, vector<double> &rho_, vector<double> &u_){
  return rho_[ind_phys_1D_] * LBMPSMLattice::w[iq_] * 
        (1.0 + ( LBMPSMLattice::e[3*iq_] * u_[ind_phys_2D_] + LBMPSMLattice::e[3*iq_+1] * u_[ind_phys_2D_+1] + LBMPSMLattice::e[3*iq_+2] * u_[ind_phys_2D_+2]) / (pow(LBMPSMLattice::cs, 2.0))
        + pow( LBMPSMLattice::e[3*iq_] * u_[ind_phys_2D_] + LBMPSMLattice::e[3*iq_+1] * u_[ind_phys_2D_+1] + LBMPSMLattice::e[3*iq_+2] * u_[ind_phys_2D_+2] , 2.0) / (2.0*pow(LBMPSMLattice::cs, 4.0))
        - 1.0/2.0 * ( u_[ind_phys_2D_] * u_[ind_phys_2D_] + u_[ind_phys_2D_+1] * u_[ind_phys_2D_+1] + u_[ind_phys_2D_+2] * u_[ind_phys_2D_+2]) / (pow(LBMPSMLattice::cs,2.0) ) );
}


double LBMPSMDynamics::feq(int iq_, double rho, vector<double> u){
  return rho * LBMPSMLattice::w[iq_] * 
        (1.0 + ( LBMPSMLattice::e[3*iq_] * u[0] + LBMPSMLattice::e[3*iq_+1] * u[1] + LBMPSMLattice::e[3*iq_+2] * u[2]) / (pow(LBMPSMLattice::cs, 2.0))
        + pow( LBMPSMLattice::e[3*iq_] * u[0] + LBMPSMLattice::e[3*iq_+1] * u[1] + LBMPSMLattice::e[iq_*3+2] * u[2], 2.0) / (2.0*pow(LBMPSMLattice::cs, 4.0))
        - 1.0/2.0 * ( u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) / (pow(LBMPSMLattice::cs,2.0) ) );
}


void LBMPSMDynamics::streamBulk(int i_, int j_, int k_, int iq_){
  LBMPSMLattice::set_fcoll( i_ + LBMPSMLattice::e[3*iq_] , j_ + LBMPSMLattice::e[3*iq_+1] , k_ + LBMPSMLattice::e[3*iq_+2], iq_, LBMPSMLattice::get_f(i_, j_, k_, iq_ ) );
};


void LBMPSMDynamics::streamBC(int i_, int j_, int k_, int iq_)
{
  if(   LBMPSMLattice::e[3*iq_]+i_ >= 0 && LBMPSMLattice::e[3*iq_]+i_ <= LBMPSMLattice::nx-1
    &&  LBMPSMLattice::e[3*iq_+1]+j_ >= 0 && LBMPSMLattice::e[3*iq_+1]+j_ <= LBMPSMLattice::ny-1
    &&  LBMPSMLattice::e[3*iq_+2]+k_ >= 0 && LBMPSMLattice::e[3*iq_+2]+k_ <= LBMPSMLattice::nz-1)
    {
      LBMPSMLattice::set_fcoll( i_ + LBMPSMLattice::e[3*iq_], j_ + LBMPSMLattice::e[3*iq_+1], k_ + LBMPSMLattice::e[3*iq_+2], iq_, LBMPSMLattice::get_f(i_, j_, k_, iq_ ) );
    }
}
