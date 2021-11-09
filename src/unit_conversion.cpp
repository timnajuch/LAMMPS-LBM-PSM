/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#include "unit_conversion.h"

Unit_Conversion::Unit_Conversion(double rhof_, double nu_, double lc_, double Re_, int N_, double tau_, int dimension_) : rhof(rhof_), nu(nu_), lc(lc_), Re(Re_), N(N_), tau(tau_), dimension(dimension_) {
   
  Uc_d = 1.0;
  Lc_d = 1.0;

  cs = 1.0/sqrt(3.0);
  csPow2 = 1.0/3.0;
  nu_lb = (tau_-0.5)*csPow2;

  dx = lc_/(double)(N_-1);
  dx_d = Lc_d/(double)N_;

  Uc = Re*nu/lc_;
  u_lb = Re_*nu_lb/(double)N_;
  //double dx_d = dx_d*u_lb/Uc_d;
  //double dt_d = lc_/Uc;
  dt_d = dx_d*u_lb/Uc_d;

  tc = lc_/Uc;

  // forceFacdoubleor is doubleo scale dp/dx (nodouble doublehe force). units of dp/dx = kg/(m^2s^2) (force F = kg/m/s^2). 
  // Hence, we need additionally "/ (pow(lc,2)*pow(dx_d,2))" to scale correctly
  // TODO check units for 2D and 3D
  //if(domain->dimension == 2){
  if(dimension == 2){
    forceFactor = rhof_ * pow(lc_,4)/pow(tc,2) * pow(dx_d,4)/pow(dt_d,2);// / (pow(lc_,3)*pow(dx_d,3));
    torqueFactor = rhof_ * pow(lc_,5)/pow(tc,2) * pow(dx_d,5)/pow(dt_d,2);// / (pow(lc_,3)*pow(dx_d,3));
  }else{
    forceFactor = rhof_ * pow(lc_,4)/pow(tc,2) * pow(dx_d,4)/pow(dt_d,2);// / (pow(lc_,3)*pow(dx_d,3));
    torqueFactor = rhof_ * pow(lc_,5)/pow(tc,2) * pow(dx_d,5)/pow(dt_d,2);// / (pow(lc_,3)*pow(dx_d,3));
  }

};

Unit_Conversion::~Unit_Conversion(){};


double Unit_Conversion::get_dx(){
  return dx;
};


double Unit_Conversion::get_Uc(){
  return Uc;
};

double Unit_Conversion::get_radius_lb(double rp){
  return rp/dx;
}

double Unit_Conversion::get_u_lb(){
  return u_lb;
};

double Unit_Conversion::get_vel_lb(double vel_phys){
  return vel_phys / Uc * u_lb;
};

double Unit_Conversion::get_freq_lb(double freq_phys){
  return freq_phys * tc * dt_d;
};

double Unit_Conversion::get_pos_lb(double pos_phys){
  //return pos_phys*N/lc;
  return pos_phys/dx;
}


double Unit_Conversion::get_forceFactor(){
  return forceFactor;
};

double Unit_Conversion::get_torqueFactor(){
  return torqueFactor;
};

std::vector<double> Unit_Conversion::get_force_lb_2D(std::vector<double> F_phys_2D_){
  std::vector<double> F_lb(3,0.0);
  F_lb[0] = F_phys_2D_[0]/forceFactor;
  F_lb[1] = F_phys_2D_[1]/forceFactor;
  F_lb[2] = F_phys_2D_[2]/forceFactor;

  return F_lb;
};

double Unit_Conversion::get_phys_time(double time_lb){
  return time_lb * dt_d * tc;
}
    
double Unit_Conversion::get_dt_d(){
  return dt_d;
}
