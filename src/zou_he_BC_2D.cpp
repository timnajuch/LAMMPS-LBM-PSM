/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#include "zou_he_BC_2D.h"

ZouHeBC2D::ZouHeBC2D (Lattice2D *lattice2D_) : lattice2D(lattice2D_) {};

ZouHeBC2D::~ZouHeBC2D () {};

void ZouHeBC2D::setZouHeVelBC2D_xn (int ix_, int iy0_, int iy1_, double ux_bc_) {
  for (int j = iy0_; j <= iy1_; ++j){
    int ind_iq = ix_ * lattice2D->ny * lattice2D->q + j * lattice2D->q; // + iq;

    double rho_tmp = 1.0/(1.0-ux_bc_)*( lattice2D->get_f(ind_iq + 0) + lattice2D->get_f(ind_iq + 2) + lattice2D->get_f(ind_iq + 4) +
                          2.0*( lattice2D->get_f(ind_iq + 3) + lattice2D->get_f(ind_iq + 6) + lattice2D->get_f(ind_iq + 7) ) );

    lattice2D->set_f(ind_iq + 5, lattice2D->get_f(ind_iq + 7) + 0.5*(lattice2D->get_f(ind_iq + 4) - lattice2D->get_f(ind_iq + 2))
                          + rho_tmp*ux_bc_/6.0 );

    lattice2D->set_f(ind_iq + 8, lattice2D->get_f(ind_iq + 6) + 0.5*(lattice2D->get_f(ind_iq + 2) - lattice2D->get_f(ind_iq + 4))
                          + rho_tmp*ux_bc_/6.0 );

    lattice2D->set_f(ind_iq + 1, lattice2D->get_f(ind_iq + 3) + 2.0/3.0*rho_tmp*ux_bc_ );
  }
};

void ZouHeBC2D::setZouHeDensBC2D_xp (int ix_, int iy0_, int iy1_, double rho_bc_) {
  for (int j = iy0_; j <= iy1_; ++j){
    int ind_iq = ix_ * lattice2D->ny * lattice2D->q + j * lattice2D->q; // + iq;

    double ux_tmp = -1.0 + 1.0/rho_bc_*(lattice2D->get_f(ind_iq + 0) + lattice2D->get_f(ind_iq + 2) + lattice2D->get_f(ind_iq + 4) 
                      + 2.0*(lattice2D->get_f(ind_iq + 1) + lattice2D->get_f(ind_iq + 5) + lattice2D->get_f(ind_iq + 8)) );

    lattice2D->set_f(ind_iq + 6, lattice2D->get_f(ind_iq + 8) + 0.5*(lattice2D->get_f(ind_iq + 4) - lattice2D->get_f(ind_iq + 2))
                          - rho_bc_*ux_tmp/6.0 );
                          
    lattice2D->set_f(ind_iq + 7, lattice2D->get_f(ind_iq + 5) + 0.5*(lattice2D->get_f(ind_iq + 2) - lattice2D->get_f(ind_iq + 4))
                          - rho_bc_*ux_tmp/6.0 );
                          
    lattice2D->set_f(ind_iq + 3, lattice2D->get_f(ind_iq + 1) - 2.0/3.0*rho_bc_*ux_tmp );
  }
};

void ZouHeBC2D::setZouHeVelBC2D_yn (int iy_, int ix0_, int ix1_, double ux_bc_) {
  for (int i = ix0_; i <= ix1_; ++i){
    int ind_iq = i * lattice2D->ny * lattice2D->q + iy_ * lattice2D->q; // + iq;

    double rho_tmp =  lattice2D->get_f(ind_iq + 0) + lattice2D->get_f(ind_iq + 3) + lattice2D->get_f(ind_iq + 1) +
                              2.0*( lattice2D->get_f(ind_iq + 8) + lattice2D->get_f(ind_iq + 4) + lattice2D->get_f(ind_iq + 7) );

    lattice2D->set_f(ind_iq + 6, lattice2D->get_f(ind_iq + 8) + 0.5*(lattice2D->get_f(ind_iq + 1) - lattice2D->get_f(ind_iq + 3))
                          - rho_tmp*ux_bc_/2.0 );

    lattice2D->set_f(ind_iq + 5, lattice2D->get_f(ind_iq + 7) + 0.5*(lattice2D->get_f(ind_iq + 3) - lattice2D->get_f(ind_iq + 1))
                          + rho_tmp*ux_bc_/2.0 );

    lattice2D->set_f(ind_iq + 2, lattice2D->get_f(ind_iq + 4) );
  }
};

void ZouHeBC2D::setZouHeVelBC2D_yp (int iy_, int ix0_, int ix1_, double ux_bc_) {
  for (int i = ix0_; i <= ix1_; ++i){
    int ind_iq = i * lattice2D->ny * lattice2D->q + iy_ * lattice2D->q; // + iq;

    double rho_tmp =  lattice2D->get_f(ind_iq + 0) + lattice2D->get_f(ind_iq + 3) + lattice2D->get_f(ind_iq + 1) +
                              2.0*( lattice2D->get_f(ind_iq + 6) + lattice2D->get_f(ind_iq + 2) + lattice2D->get_f(ind_iq + 5) );

    lattice2D->set_f(ind_iq + 7, lattice2D->get_f(ind_iq + 5) + 0.5*(lattice2D->get_f(ind_iq + 1) - lattice2D->get_f(ind_iq + 3))
                          - rho_tmp*ux_bc_/2.0 );

    lattice2D->set_f(ind_iq + 8, lattice2D->get_f(ind_iq + 6) + 0.5*(lattice2D->get_f(ind_iq + 3) - lattice2D->get_f(ind_iq + 1))
                          + rho_tmp*ux_bc_/2.0 );

    lattice2D->set_f(ind_iq + 4, lattice2D->get_f(ind_iq + 2) );
  }
};
