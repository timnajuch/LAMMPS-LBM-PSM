/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#ifndef LBMPSMDYNAMICS_2D
#define LBMPSMDYNAMICS_2D

#include <iostream>
#include <math.h>
#include <vector>

#include "LBM_PSM_lattice.h"

using namespace std;

class LBMPSMDynamics : public LBMPSMLattice {
  
  public:
    LBMPSMDynamics(int nx_, int ny_, int nz_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_);
    ~LBMPSMDynamics();

    inline double feq(int iq_, int ind_phys_1D_, int ind_phys_2D_); // Accesses the density and velocities stored in the Lattice class
    inline double feq(int iq_, double rho_, vector<double> u_);
    inline double feq(int iq, double rho, const std::array<double,3>& u);


    // External force according to Guo et al. (2002). When called needs to be multiplied by (1-0.5/tau)
    double F_iq(int iq_, vector<double> u_, vector<double> F_);
    double F_iq(int iq_, int ind_phys_2D_, vector<double> F_); // Accesses the velocities stored in the Lattice class

};



inline double LBMPSMDynamics::feq(int iq_, int ind_phys_1D_, int ind_phys_2D_){
  double eiq_dot_u = e[3*iq_] * u[ind_phys_2D_] + e[3*iq_+1] * u[ind_phys_2D_+1] + e[3*iq_+2] * u[ind_phys_2D_+2];
  return rho[ind_phys_1D_] * w[iq_] *
        (1.0 + eiq_dot_u * invCsPow2
        + eiq_dot_u*eiq_dot_u * 0.5 * invCsPow4
        - 0.5 * ( u[ind_phys_2D_] * u[ind_phys_2D_] + u[ind_phys_2D_+1] * u[ind_phys_2D_+1] + u[ind_phys_2D_+2] * u[ind_phys_2D_+2]) * invCsPow2 );
}


inline double LBMPSMDynamics::feq(int iq_, double rho_, vector<double> u_){
  double eiq_dot_u = e[3*iq_] * u_[0] + e[3*iq_+1] * u_[1] + e[3*iq_+2] * u_[2];
  return rho_ * w[iq_] * 
        (1.0 + eiq_dot_u * invCsPow2
        + eiq_dot_u*eiq_dot_u * 0.5 * invCsPow4
        - 0.5 * ( u_[0] * u_[0] + u_[1] * u_[1] + u_[2] * u_[2]) * invCsPow2 );



inline double LBMPSMDynamics::feq(int iq, double rho, const std::array<double,3>& u) {
  const double ex = e[3*iq];
  const double ey = e[3*iq+1];
  const double ez = e[3*iq+2];

  const double eiq_dot_u = ex*u[0] + ey*u[1] + ez*u[2];
  const double u_dot_u = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];

  return rho * w[iq] *
      (1.0
      + eiq_dot_u * invCsPow2
      + 0.5 * eiq_dot_u * eiq_dot_u * invCsPow4
      - 0.5 * u_dot_u * invCsPow2);
}


#endif

