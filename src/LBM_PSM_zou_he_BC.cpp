/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#include "LBM_PSM_zou_he_BC.h"

ZouHeBC::ZouHeBC (LBMPSMLattice *lattice_) : lattice(lattice_)
{
  oneOverThree = 1.0/3.0;
  twoOverThree = 2.0/3.0;
  oneOverSix   = 1.0/6.0;
}

ZouHeBC::~ZouHeBC () {}


void ZouHeBC::setZouHeVelBC2D_xn (int ix_, int iy0_, int iy1_, double ux_bc_, int currentStep) {
  for (int j = iy0_; j <= iy1_; ++j){
    int ind_iq = lattice->index_fi(ix_, j, 0, 0, currentStep);

    double rho_tmp = 1.0/(1.0-ux_bc_)*( lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 4) +
                          2.0*( lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 7) ) );

    lattice->set_f(ind_iq + 5, lattice->get_f(ind_iq + 7) + 0.5*(lattice->get_f(ind_iq + 4) - lattice->get_f(ind_iq + 2))
                          + rho_tmp*ux_bc_*oneOverSix );

    lattice->set_f(ind_iq + 8, lattice->get_f(ind_iq + 6) + 0.5*(lattice->get_f(ind_iq + 2) - lattice->get_f(ind_iq + 4))
                          + rho_tmp*ux_bc_*oneOverSix );

    lattice->set_f(ind_iq + 1, lattice->get_f(ind_iq + 3) + twoOverThree*rho_tmp*ux_bc_ );
  }
}


void ZouHeBC::setZouHeVelBC2D_xp (int ix_, int iy0_, int iy1_, double ux_bc_, int currentStep) {
  for (int j = iy0_; j <= iy1_; ++j){
    int ind_iq = lattice->index_fi(ix_, j, 0, 0, currentStep);

    double rho_tmp = 1.0/(1.0-ux_bc_)*( lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 4) +
                          2.0*( lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 8) ) );

    lattice->set_f(ind_iq + 7, lattice->get_f(ind_iq + 5) + 0.5*(lattice->get_f(ind_iq + 2) - lattice->get_f(ind_iq + 4))
                          - rho_tmp*ux_bc_*oneOverSix );

    lattice->set_f(ind_iq + 6, lattice->get_f(ind_iq + 8) + 0.5*(lattice->get_f(ind_iq + 4) - lattice->get_f(ind_iq + 2))
                          - rho_tmp*ux_bc_*oneOverSix );

    lattice->set_f(ind_iq + 3, lattice->get_f(ind_iq + 1) - twoOverThree*rho_tmp*ux_bc_ );
  }
}


void ZouHeBC::setZouHeDensBC2D_xp (int ix_, int iy0_, int iy1_, double rho_bc_, int currentStep) {
  for (int j = iy0_; j <= iy1_; ++j){
    int ind_iq = lattice->index_fi(ix_, j, 0, 0, currentStep);

    double ux_tmp = -1.0 + 1.0/rho_bc_*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 4) 
                      + 2.0*(lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 8)) );

    lattice->set_f(ind_iq + 6, lattice->get_f(ind_iq + 8) + 0.5*(lattice->get_f(ind_iq + 4) - lattice->get_f(ind_iq + 2))
                          - rho_bc_*ux_tmp*oneOverSix );
                          
    lattice->set_f(ind_iq + 7, lattice->get_f(ind_iq + 5) + 0.5*(lattice->get_f(ind_iq + 2) - lattice->get_f(ind_iq + 4))
                          - rho_bc_*ux_tmp*oneOverSix );
                          
    lattice->set_f(ind_iq + 3, lattice->get_f(ind_iq + 1) - twoOverThree*rho_bc_*ux_tmp );
  }
}


void ZouHeBC::setZouHeVelBC2D_yn (int iy_, int ix0_, int ix1_, double ux_bc_, int currentStep) {
  for (int i = ix0_; i <= ix1_; ++i){
    int ind_iq = lattice->index_fi(i, iy_, 0, 0, currentStep);

    double rho_tmp =  lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 1) +
                              2.0*( lattice->get_f(ind_iq + 8) + lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 7) );

    lattice->set_f(ind_iq + 6, lattice->get_f(ind_iq + 8) + 0.5*(lattice->get_f(ind_iq + 1) - lattice->get_f(ind_iq + 3))
                          - rho_tmp*ux_bc_*0.5 );

    lattice->set_f(ind_iq + 5, lattice->get_f(ind_iq + 7) + 0.5*(lattice->get_f(ind_iq + 3) - lattice->get_f(ind_iq + 1))
                          + rho_tmp*ux_bc_*0.5 );

    lattice->set_f(ind_iq + 2, lattice->get_f(ind_iq + 4) );
  }
}


void ZouHeBC::setZouHeVelBC2D_yp (int iy_, int ix0_, int ix1_, double ux_bc_, int currentStep) {
  for (int i = ix0_; i <= ix1_; ++i){
    int ind_iq = lattice->index_fi(i, iy_, 0, 0, currentStep);

    double rho_tmp =  lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 1) +
                              2.0*( lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 5) );

    lattice->set_f(ind_iq + 7, lattice->get_f(ind_iq + 5) + 0.5*(lattice->get_f(ind_iq + 1) - lattice->get_f(ind_iq + 3))
                          - rho_tmp*ux_bc_*0.5 );

    lattice->set_f(ind_iq + 8, lattice->get_f(ind_iq + 6) + 0.5*(lattice->get_f(ind_iq + 3) - lattice->get_f(ind_iq + 1))
                          + rho_tmp*ux_bc_*0.5 );

    lattice->set_f(ind_iq + 4, lattice->get_f(ind_iq + 2) );
  }
}


void ZouHeBC::setZouHeNeumannVelBC2D_yn (int iy_, int ix0_, int ix1_, int currentStep) {
  for (int i = ix0_; i <= ix1_; ++i){
    int ind_iq = lattice->index_fi(i, iy_, 0, 0, currentStep);
    int ind_u_neighbour_1D = lattice->index_1D(i, iy_+1, 0);

    double rho_tmp =  lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 1) +
                              2.0*( lattice->get_f(ind_iq + 8) + lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 7) );

    lattice->set_f(ind_iq + 6, lattice->get_f(ind_iq + 8) + 0.5*(lattice->get_f(ind_iq + 1) - lattice->get_f(ind_iq + 3))
                          - rho_tmp*lattice->get_u_at_node(ind_u_neighbour_1D, 0)*0.5 );

    lattice->set_f(ind_iq + 5, lattice->get_f(ind_iq + 7) + 0.5*(lattice->get_f(ind_iq + 3) - lattice->get_f(ind_iq + 1))
                          + rho_tmp*lattice->get_u_at_node(ind_u_neighbour_1D, 0)*0.5 );

    lattice->set_f(ind_iq + 2, lattice->get_f(ind_iq + 4) );
  }
}


void ZouHeBC::setZouHeNeumannVelBC2D_yp (int iy_, int ix0_, int ix1_, int currentStep) {
  for (int i = ix0_; i <= ix1_; ++i){
    int ind_iq = lattice->index_fi(i, iy_, 0, 0, currentStep);
    int ind_u_neighbour_1D = lattice->index_1D(i, iy_-1, 0);

    double rho_tmp =  lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 1) +
                              2.0*( lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 5) );

    lattice->set_f(ind_iq + 7, lattice->get_f(ind_iq + 5) + 0.5*(lattice->get_f(ind_iq + 1) - lattice->get_f(ind_iq + 3))
                          - rho_tmp*lattice->get_u_at_node(ind_u_neighbour_1D, 0)*0.5 );

    lattice->set_f(ind_iq + 8, lattice->get_f(ind_iq + 6) + 0.5*(lattice->get_f(ind_iq + 3) - lattice->get_f(ind_iq + 1))
                          + rho_tmp*lattice->get_u_at_node(ind_u_neighbour_1D, 0)*0.5 );

    lattice->set_f(ind_iq + 4, lattice->get_f(ind_iq + 2) );
  }
}


void ZouHeBC::setZouHeVelBC3D_xn (int ix_, int iy0_, int iy1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_, int currentStep) {
  for (int j = iy0_; j <= iy1_; ++j){
    for (int k = iz0_; k <= iz1_; ++k){
      int ind_iq = lattice->index_fi(ix_, j, k, 0, currentStep);

      double rho_tmp =  1.0/(1.0-ux_bc_)*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 4) +
                                          lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 11) +
                                          lattice->get_f(ind_iq + 17) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 12) +
                                2.0*( lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 14) + lattice->get_f(ind_iq + 8) +
                                      lattice->get_f(ind_iq + 16) + lattice->get_f(ind_iq + 10) ) );

      double Nyx = 0.5*(lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 11) + lattice->get_f(ind_iq + 17)
                        -(lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 12)))
                   - oneOverThree*rho_tmp*uy_bc_;

      double Nzx = 0.5*(lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 11)
                        -(lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 17) + lattice->get_f(ind_iq + 12)))
                   - oneOverThree*rho_tmp*uz_bc_;


      lattice->set_f(ind_iq + 1, lattice->get_f(ind_iq + 2) + oneOverThree*rho_tmp*ux_bc_ );

      lattice->set_f(ind_iq + 13, lattice->get_f(ind_iq + 14) + rho_tmp*oneOverSix*(ux_bc_ - uy_bc_) + Nyx );

      lattice->set_f(ind_iq + 7, lattice->get_f(ind_iq + 8) + rho_tmp*oneOverSix*(ux_bc_ + uy_bc_) - Nyx );

      lattice->set_f(ind_iq + 9, lattice->get_f(ind_iq + 10) + rho_tmp*oneOverSix*(ux_bc_ + uz_bc_) - Nzx );

      lattice->set_f(ind_iq + 15, lattice->get_f(ind_iq + 16) + rho_tmp*oneOverSix*(ux_bc_ - uz_bc_) + Nzx );

    }
  }
}



void ZouHeBC::setZouHeVelBC3D_xp (int ix_, int iy0_, int iy1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_, int currentStep) {
  for (int j = iy0_; j <= iy1_; ++j){
    for (int k = iz0_; k <= iz1_; ++k){
      int ind_iq = lattice->index_fi(ix_, j, k, 0, currentStep);

      double rho_tmp =  1.0/(1.0+ux_bc_)*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 4) +
                                          lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 11) +
                                          lattice->get_f(ind_iq + 17) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 12) +
                                2.0*( lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 7) + lattice->get_f(ind_iq + 13) +
                                      lattice->get_f(ind_iq + 9) + lattice->get_f(ind_iq + 15) ) );

      double Nyx = 0.5*(lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 11) + lattice->get_f(ind_iq + 17)
                        -(lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 12)))
                   - oneOverThree*rho_tmp*uy_bc_;

      double Nzx = 0.5*(lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 11)
                        -(lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 17) + lattice->get_f(ind_iq + 12)))
                   - oneOverThree*rho_tmp*uz_bc_;


      lattice->set_f(ind_iq + 2, lattice->get_f(ind_iq + 1) - oneOverThree*rho_tmp*ux_bc_ );

      lattice->set_f(ind_iq + 14, lattice->get_f(ind_iq + 13) + rho_tmp*oneOverSix*(-ux_bc_ + uy_bc_) - Nyx );

      lattice->set_f(ind_iq + 8, lattice->get_f(ind_iq + 7) + rho_tmp*oneOverSix*(-ux_bc_ - uy_bc_) + Nyx );

      lattice->set_f(ind_iq + 10, lattice->get_f(ind_iq + 9) + rho_tmp*oneOverSix*(-ux_bc_ - uz_bc_) + Nzx );

      lattice->set_f(ind_iq + 16, lattice->get_f(ind_iq + 15) + rho_tmp*oneOverSix*(-ux_bc_ + uz_bc_) - Nzx );

    }
  }
}



void ZouHeBC::setZouHeDensBC3D_xp (int ix_, int iy0_, int iy1_, int iz0_, int iz1_, double rho_bc_, double uy_bc_, double uz_bc_, int currentStep) {
  for (int j = iy0_; j <= iy1_; ++j){
    for (int k = iz0_; k <= iz1_; ++k){
      int ind_iq = lattice->index_fi(ix_, j, k, 0, currentStep);

      double ux_bc_ =  -1.0 + 1.0/rho_bc_*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 4) +
                                          lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 11) +
                                          lattice->get_f(ind_iq + 17) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 12) +
                                2.0*( lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 7) + lattice->get_f(ind_iq + 13) +
                                      lattice->get_f(ind_iq + 9) + lattice->get_f(ind_iq + 15) ) );

      double Nyx = 0.5*(lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 11) + lattice->get_f(ind_iq + 17)
                        -(lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 12)))
                   - oneOverThree*rho_bc_*uy_bc_;

      double Nzx = 0.5*(lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 11)
                        -(lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 17) + lattice->get_f(ind_iq + 12)))
                   - oneOverThree*rho_bc_*uz_bc_;


      lattice->set_f(ind_iq + 2, lattice->get_f(ind_iq + 1) - oneOverThree*rho_bc_*ux_bc_ );

      lattice->set_f(ind_iq + 14, lattice->get_f(ind_iq + 13) + rho_bc_*oneOverSix*(-ux_bc_ + uy_bc_) - Nyx );

      lattice->set_f(ind_iq + 8, lattice->get_f(ind_iq + 7) + rho_bc_*oneOverSix*(-ux_bc_ - uy_bc_) + Nyx );

      lattice->set_f(ind_iq + 10, lattice->get_f(ind_iq + 9) + rho_bc_*oneOverSix*(-ux_bc_ - uz_bc_) + Nzx );

      lattice->set_f(ind_iq + 16, lattice->get_f(ind_iq + 15) + rho_bc_*oneOverSix*(-ux_bc_ + uz_bc_) - Nzx );

    }
  }
}



void ZouHeBC::setZouHeVelBC3D_yn (int iy_, int ix0_, int ix1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_, int currentStep) {
  for (int i = ix0_; i <= ix1_; ++i){
    for (int k = iz0_; k <= iz1_; ++k){
      int ind_iq = lattice->index_fi(i, iy_, k, 0, currentStep);

      double rho_tmp =  1.0/(1.0-uy_bc_)*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 2) +
                                          lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 9) +
                                          lattice->get_f(ind_iq + 15) + lattice->get_f(ind_iq + 16) + lattice->get_f(ind_iq + 10) +
                                2.0*( lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 13) + lattice->get_f(ind_iq + 8) +
                                      lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 12) ) );

      double Nxy = 0.5*(lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 9) + lattice->get_f(ind_iq + 15)
                        -(lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 16) + lattice->get_f(ind_iq + 10)))
                   - oneOverThree*rho_tmp*ux_bc_;

      double Nzy = 0.5*(lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 9) + lattice->get_f(ind_iq + 16)
                        -(lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 15) + lattice->get_f(ind_iq + 10)))
                   - oneOverThree*rho_tmp*uz_bc_;


      lattice->set_f(ind_iq + 3, lattice->get_f(ind_iq + 4) + oneOverThree*rho_tmp*uy_bc_ );

      lattice->set_f(ind_iq + 7, lattice->get_f(ind_iq + 8) + rho_tmp*oneOverSix*(uy_bc_ + ux_bc_) - Nxy );

      lattice->set_f(ind_iq + 14, lattice->get_f(ind_iq + 13) + rho_tmp*oneOverSix*(uy_bc_ - ux_bc_) + Nxy );

      lattice->set_f(ind_iq + 11, lattice->get_f(ind_iq + 12) + rho_tmp*oneOverSix*(uy_bc_ + uz_bc_) - Nzy );

      lattice->set_f(ind_iq + 17, lattice->get_f(ind_iq + 18) + rho_tmp*oneOverSix*(uy_bc_ - uz_bc_) + Nzy );

    }
  }
}


void ZouHeBC::setZouHeVelBC3D_yp (int iy_, int ix0_, int ix1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_, int currentStep) {
  for (int i = ix0_; i <= ix1_; ++i){
    for (int k = iz0_; k <= iz1_; ++k){
      int ind_iq = lattice->index_fi(i, iy_, k, 0, currentStep);

      double rho_tmp =  1.0/(1.0+uy_bc_)*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 2) +
                                          lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 9) +
                                          lattice->get_f(ind_iq + 15) + lattice->get_f(ind_iq + 16) + lattice->get_f(ind_iq + 10) +
                                2.0*( lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 7) + lattice->get_f(ind_iq + 14) +
                                      lattice->get_f(ind_iq + 11) + lattice->get_f(ind_iq + 17) ) );

      double Nxy = 0.5*(lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 9) + lattice->get_f(ind_iq + 15)
                        -(lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 16) + lattice->get_f(ind_iq + 10)))
                   - oneOverThree*rho_tmp*ux_bc_;

      double Nzy = 0.5*(lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 9) + lattice->get_f(ind_iq + 16)
                        -(lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 15) + lattice->get_f(ind_iq + 10)))
                   - oneOverThree*rho_tmp*uz_bc_;


      lattice->set_f(ind_iq + 4, lattice->get_f(ind_iq + 3) - oneOverThree*rho_tmp*uy_bc_ );

      lattice->set_f(ind_iq + 8, lattice->get_f(ind_iq + 7) + rho_tmp*oneOverSix*(-uy_bc_ - ux_bc_) + Nxy );

      lattice->set_f(ind_iq + 13, lattice->get_f(ind_iq + 14) + rho_tmp*oneOverSix*(-uy_bc_ + ux_bc_) - Nxy );

      lattice->set_f(ind_iq + 12, lattice->get_f(ind_iq + 11) + rho_tmp*oneOverSix*(-uy_bc_ - uz_bc_) + Nzy );

      lattice->set_f(ind_iq + 18, lattice->get_f(ind_iq + 17) + rho_tmp*oneOverSix*(-uy_bc_ + uz_bc_) - Nzy );

    }
  }
}


void ZouHeBC::setZouHeVelBC3D_zn (int iz_, int ix0_, int ix1_, int iy0_, int iy1_, double ux_bc_, double uy_bc_, double uz_bc_, int currentStep) {
  for (int i = ix0_; i <= ix1_; ++i){
    for (int j = iy0_; j <= iy1_; ++j){
      int ind_iq = lattice->index_fi(i, j, iz_, 0, currentStep);

      double rho_tmp =  1.0/(1.0-uz_bc_)*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 2) +
                                          lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 7) +
                                          lattice->get_f(ind_iq + 13) + lattice->get_f(ind_iq + 14) + lattice->get_f(ind_iq + 8) +
                                2.0*( lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 15) + lattice->get_f(ind_iq + 10) +
                                      lattice->get_f(ind_iq + 17) + lattice->get_f(ind_iq + 12) ) );

      double Nxz = 0.5*(lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 7) + lattice->get_f(ind_iq + 13)
                        -(lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 14) + lattice->get_f(ind_iq + 8)))
                   - oneOverThree*rho_tmp*ux_bc_;

      double Nyz = 0.5*(lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 7) + lattice->get_f(ind_iq + 14)
                        -(lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 13) + lattice->get_f(ind_iq + 8)))
                   - oneOverThree*rho_tmp*uy_bc_;


      lattice->set_f(ind_iq + 5, lattice->get_f(ind_iq + 6) + oneOverThree*rho_tmp*uz_bc_ );

      lattice->set_f(ind_iq + 9, lattice->get_f(ind_iq + 10) + rho_tmp*oneOverSix*(uz_bc_ + ux_bc_) - Nxz );

      lattice->set_f(ind_iq + 16, lattice->get_f(ind_iq + 15) + rho_tmp*oneOverSix*(uz_bc_ - ux_bc_) + Nxz );

      lattice->set_f(ind_iq + 11, lattice->get_f(ind_iq + 12) + rho_tmp*oneOverSix*(uz_bc_ + uy_bc_) - Nyz );

      lattice->set_f(ind_iq + 18, lattice->get_f(ind_iq + 17) + rho_tmp*oneOverSix*(uz_bc_ - uy_bc_) + Nyz );

    }
  }
}


void ZouHeBC::setZouHeVelBC3D_zp (int iz_, int ix0_, int ix1_, int iy0_, int iy1_, double ux_bc_, double uy_bc_, double uz_bc_, int currentStep) {
  for (int i = ix0_; i <= ix1_; ++i){
    for (int j = iy0_; j <= iy1_; ++j){
      int ind_iq = lattice->index_fi(i, j, iz_, 0, currentStep);

      double rho_tmp =  1.0/(1.0+uz_bc_)*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 2) +
                                          lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 7) +
                                          lattice->get_f(ind_iq + 14) + lattice->get_f(ind_iq + 8) + lattice->get_f(ind_iq + 13) +
                                2.0*( lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 9) + lattice->get_f(ind_iq + 16) +
                                      lattice->get_f(ind_iq + 11) + lattice->get_f(ind_iq + 18) ) );

      double Nxz = 0.5*(lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 7) + lattice->get_f(ind_iq + 13)
                        -(lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 14) + lattice->get_f(ind_iq + 8)))
                   - oneOverThree*rho_tmp*ux_bc_;

      double Nyz = 0.5*(lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 7) + lattice->get_f(ind_iq + 14)
                        -(lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 13) + lattice->get_f(ind_iq + 8)))
                   - oneOverThree*rho_tmp*uy_bc_;


      lattice->set_f(ind_iq + 6, lattice->get_f(ind_iq + 5) - oneOverThree*rho_tmp*uz_bc_ );

      lattice->set_f(ind_iq + 15, lattice->get_f(ind_iq + 16) + rho_tmp*oneOverSix*(-uz_bc_ + ux_bc_) - Nxz );

      lattice->set_f(ind_iq + 10, lattice->get_f(ind_iq + 9) + rho_tmp*oneOverSix*(-uz_bc_ - ux_bc_) + Nxz );

      lattice->set_f(ind_iq + 17, lattice->get_f(ind_iq + 18) + rho_tmp*oneOverSix*(-uz_bc_ + uy_bc_) - Nyz );

      lattice->set_f(ind_iq + 12, lattice->get_f(ind_iq + 11) + rho_tmp*oneOverSix*(-uz_bc_ - uy_bc_) + Nyz );

    }
  }
}
