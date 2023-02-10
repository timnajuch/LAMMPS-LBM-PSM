/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#ifndef LBM_PSM_UNIT_CONVERSION_H
#define LBM_PSM_UNIT_CONVERSION_H

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

class UnitConversion{
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

    // Factors to transfer between physical and LBM system    
    double forceFactor;
    double torqueFactor;
    double volumeForceFactor;
  
    int dimension;

  public:
    UnitConversion(double rhof_, double nu_, double lc_, double Re_, int N_, double tau_, int dimension_);
    ~UnitConversion();

    double get_dx();
    double get_Uc();
    double get_radius_lb(double rp);
    double get_u_lb();
    double get_vel_lb(double vel_phys);
    double get_freq_lb(double freq_phys);
    double get_pos_lb(double pos_phys);
    double get_forceFactor();
    double get_torqueFactor();
    double get_volumeForceFactor();
    std::vector<double> get_volume_force_lb(std::vector<double> F_phys);
    double get_phys_time(double time_lb);
    double get_dt_d();

};

#endif
