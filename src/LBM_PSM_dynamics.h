/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch
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
    inline double F_iq(int iq_, vector<double> u_, vector<double> F_);
    inline double F_iq(int iq_, int ind_phys_2D_, vector<double> F_); // Accesses the velocities stored in the Lattice class

};



inline double LBMPSMDynamics::feq(int iq_, int ind_phys_1D_, int ind_phys_2D_){
  const double ex_i = ex[iq_];
  const double ey_i = ey[iq_];
  const double ez_i = ez[iq_];
  const double w_i = w[iq_];

  const double ux =  u[ind_phys_2D_];
  const double uy =  u[ind_phys_2D_+1];
  const double uz =  u[ind_phys_2D_+2];

  double eiq_dot_u = ex_i * ux + ey_i * uy + ez_i * uz;

  return rho[ind_phys_1D_] * w_i *
        (1.0 + eiq_dot_u * invCsPow2
        + eiq_dot_u*eiq_dot_u * 0.5 * invCsPow4
        - 0.5 * ( ux * ux + uy * uy + uz * uz) * invCsPow2 );
}


inline double LBMPSMDynamics::feq(int iq_, double rho_, vector<double> u_){
  double eiq_dot_u = ex[iq_] * u_[0] + ey[iq_] * u_[1] + ez[iq_] * u_[2];

  return rho_ * w[iq_] * 
        (1.0 + eiq_dot_u * invCsPow2
        + eiq_dot_u*eiq_dot_u * 0.5 * invCsPow4
        - 0.5 * ( u_[0] * u_[0] + u_[1] * u_[1] + u_[2] * u_[2]) * invCsPow2 );
}


inline double LBMPSMDynamics::feq(int iq, double rho, const std::array<double,3>& u) {
  const double ex_i = ex[iq];
  const double ey_i = ey[iq];
  const double ez_i = ez[iq];

  const double eiq_dot_u = ex_i*u[0] + ey_i*u[1] + ez_i*u[2];
  const double u_dot_u = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];

  return rho * w[iq] *
      (1.0
      + eiq_dot_u * invCsPow2
      + 0.5 * eiq_dot_u * eiq_dot_u * invCsPow4
      - 0.5 * u_dot_u * invCsPow2);
}



inline double LBMPSMDynamics::F_iq(int iq_, vector<double> u_, vector<double> F_){
  const double ex_i = ex[iq_];
  const double ey_i = ey[iq_];
  const double ez_i = ez[iq_];
  const double w_i = w[iq_];

  return w_i * 
        ( ( (ex_i - u_[0]) * F_[0] + (ey_i - u_[1]) * F_[1] + (ez_i - u_[2]) * F_[2]) * invCsPow2
          + ( (ex_i * u_[0] + ey_i * u_[1] + ez_i * u_[2]) * invCsPow4
            * (ex_i * F_[0] + ey_i * F_[1] + ez_i * F_[2]) ) );
}


inline double LBMPSMDynamics::F_iq(int iq_, int ind_phys_2D_, vector<double> F_){
  double ex_i = ex[iq_];
  double ey_i = ey[iq_];
  double ez_i = ez[iq_];
  double w_i = w[iq_];

  double ux =  u[ind_phys_2D_];
  double uy =  u[ind_phys_2D_+1];
  double uz =  u[ind_phys_2D_+2];

  return w_i * 
        ( ( (ex_i - ux) * F_[0] + (ey_i - uy) * F_[1] + (ez_i - uz) * F_[2]) * invCsPow2
          + ( (ex_i * ux + ey_i * uy + ez_i * uz) * invCsPow4
            * (ex_i * F_[0] + ey_i * F_[1] + ez_i * F_[2]) ) );
}


#endif

