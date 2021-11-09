/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#ifndef UNIT_CONVERSION_H
#define UNIT_CONVERSION_H

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

class Unit_Conversion{
  private:
    // Input parameters (read in via constructor)
    // Fluid density in physical system
    double rhof;
    // Fluid viscosity in physical system
    double nu;
    // Characteristic lenght in physical system
    double lc;
    // Reynolds number
    double Re;
    // Number of lattice point over the characteristic length
    int N;
    // LBM relaxation parameter tau
    double tau;

    // Fixed dimensionless parameters and LB parameters (defined in constructor)
    //Characterisitic velocity in the dimensionless system
    double Uc_d;
    // Characteristic length in the dimensionless system
    double Lc_d;
    // Sound of speed in LB system
    double cs;           
    double csPow2;
    
    // Calculated parameters and converted characteristic values or parameters (converted to dimensionless or LB system)
    // Fluid viscosity in LB system
    double nu_lb;
    // lattice cell width in physical system
    double dx;
    // lattice cell width in dimensionless system
    double dx_d;
    // Characteristic velocity in physical system
    double Uc;
    // Characteristic LB velocity
    double u_lb;
    // Time increment in dimensionless system
    double dt_d;
    // Characteristic time in physical system
    double tc;    

    // Factor to transfer between physical and LBM system
    // forceFactor is to scale dp/dx (not the force). units of dp/dx = kg/(m^2s^2) (force F = kg/m/s^2). 
    // Hence, we need additionally "/ (pow(ly,2)*pow(dy_d,2))" to scale correctly
    double forceFactor;
    double torqueFactor;
    // Force converted to LB system
    //double F_lb[2];

    int dimension;

  public:
    Unit_Conversion(double rhof_, double nu_, double lc_, double Re_, int N_, double tau_, int dimension_);
    ~Unit_Conversion();

    double get_dx();
    double get_Uc();
    double get_radius_lb(double rp);
    double get_u_lb();
    double get_vel_lb(double vel_phys);
    double get_freq_lb(double freq_phys);
    double get_pos_lb(double pos_phys);
    double get_forceFactor();
    double get_torqueFactor();
    std::vector<double> get_force_lb_2D(std::vector<double> F_phys_2D_);
    double get_phys_time(double time_lb);
    double get_dt_d();

};

#endif
