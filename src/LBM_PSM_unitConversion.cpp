/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#include "LBM_PSM_unitConversion.h"

UnitConversion::UnitConversion(double rhof_, double nu_, double lc_, double Re_, int N_, double tau_, int dimension_) : rhof(rhof_), nu(nu_), lc(lc_), Re(Re_), N(N_), tau(tau_), dimension(dimension_) {
   
  Uc_d = 1.0;
  Lc_d = 1.0;

  cs = 1.0/sqrt(3.0);
  csPow2 = 1.0/3.0;
  nu_lb = (tau_-0.5)*csPow2;

  dx = lc_/(double)(N_-1);
  dx_d = Lc_d/(double)(N_-1);

  Uc = Re_*nu_/lc_;
  u_lb = Re_*nu_lb/(double)(N_-1);
  dt_d = dx_d*u_lb/Uc_d;

  tc = lc_/Uc;

  // TODO check units for 2D and 3D
  if(dimension == 2){
    forceFactor = rhof_ * pow(lc_,4)/pow(tc,2) * pow(dx_d,4)/pow(dt_d,2);
    torqueFactor = rhof_ * pow(lc_,5)/pow(tc,2) * pow(dx_d,5)/pow(dt_d,2);
  }else{
    forceFactor = rhof_ * pow(lc_,4)/pow(tc,2) * pow(dx_d,4)/pow(dt_d,2);
    torqueFactor = rhof_ * pow(lc_,5)/pow(tc,2) * pow(dx_d,5)/pow(dt_d,2);
  }
};

UnitConversion::~UnitConversion(){};


double UnitConversion::get_dx(){
  return dx;
};


double UnitConversion::get_Uc(){
  return Uc;
};

double UnitConversion::get_radius_lb(double rp){
  return rp/dx;
}

double UnitConversion::get_u_lb(){
  return u_lb;
};

double UnitConversion::get_vel_lb(double vel_phys){
  return vel_phys / Uc * u_lb;
};

double UnitConversion::get_freq_lb(double freq_phys){
  return freq_phys * tc * dt_d;
};

double UnitConversion::get_pos_lb(double pos_phys){
  return pos_phys/dx;
}


double UnitConversion::get_forceFactor(){
  return forceFactor;
};

double UnitConversion::get_torqueFactor(){
  return torqueFactor;
};

std::vector<double> UnitConversion::get_force_lb_2D(std::vector<double> F_phys_2D_){
  std::vector<double> F_lb(3,0.0);
  F_lb[0] = F_phys_2D_[0]/forceFactor;
  F_lb[1] = F_phys_2D_[1]/forceFactor;
  F_lb[2] = F_phys_2D_[2]/forceFactor;

  return F_lb;
};

double UnitConversion::get_phys_time(double time_lb){
  return time_lb * dt_d * tc;
}
    
double UnitConversion::get_dt_d(){
  return dt_d;
}
