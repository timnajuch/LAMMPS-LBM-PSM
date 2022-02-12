/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#include "LBM_PSM_zou_he_BC.h"

ZouHeBC::ZouHeBC (LBMPSMLattice *lattice_) : lattice(lattice_) {};

ZouHeBC::~ZouHeBC () {};


void ZouHeBC::setZouHeVelBC2D_xn (int ix_, int iy0_, int iy1_, double ux_bc_) {
  for (int j = iy0_; j <= iy1_; ++j){
    int ind_iq = ix_ * lattice->ny * lattice->q + j * lattice->q; // + iq;

    double rho_tmp = 1.0/(1.0-ux_bc_)*( lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 4) +
                          2.0*( lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 7) ) );

    lattice->set_f(ind_iq + 5, lattice->get_f(ind_iq + 7) + 0.5*(lattice->get_f(ind_iq + 4) - lattice->get_f(ind_iq + 2))
                          + rho_tmp*ux_bc_/6.0 );

    lattice->set_f(ind_iq + 8, lattice->get_f(ind_iq + 6) + 0.5*(lattice->get_f(ind_iq + 2) - lattice->get_f(ind_iq + 4))
                          + rho_tmp*ux_bc_/6.0 );

    lattice->set_f(ind_iq + 1, lattice->get_f(ind_iq + 3) + 2.0/3.0*rho_tmp*ux_bc_ );
  }
};


void ZouHeBC::setZouHeDensBC2D_xp (int ix_, int iy0_, int iy1_, double rho_bc_) {
  for (int j = iy0_; j <= iy1_; ++j){
    int ind_iq = ix_ * lattice->ny * lattice->q + j * lattice->q; // + iq;

    double ux_tmp = -1.0 + 1.0/rho_bc_*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 4) 
                      + 2.0*(lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 8)) );

    lattice->set_f(ind_iq + 6, lattice->get_f(ind_iq + 8) + 0.5*(lattice->get_f(ind_iq + 4) - lattice->get_f(ind_iq + 2))
                          - rho_bc_*ux_tmp/6.0 );
                          
    lattice->set_f(ind_iq + 7, lattice->get_f(ind_iq + 5) + 0.5*(lattice->get_f(ind_iq + 2) - lattice->get_f(ind_iq + 4))
                          - rho_bc_*ux_tmp/6.0 );
                          
    lattice->set_f(ind_iq + 3, lattice->get_f(ind_iq + 1) - 2.0/3.0*rho_bc_*ux_tmp );
  }
};


void ZouHeBC::setZouHeVelBC2D_yn (int iy_, int ix0_, int ix1_, double ux_bc_) {
  for (int i = ix0_; i <= ix1_; ++i){
    int ind_iq = i * lattice->ny * lattice->q + iy_ * lattice->q; // + iq;

    double rho_tmp =  lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 1) +
                              2.0*( lattice->get_f(ind_iq + 8) + lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 7) );

    lattice->set_f(ind_iq + 6, lattice->get_f(ind_iq + 8) + 0.5*(lattice->get_f(ind_iq + 1) - lattice->get_f(ind_iq + 3))
                          - rho_tmp*ux_bc_/2.0 );

    lattice->set_f(ind_iq + 5, lattice->get_f(ind_iq + 7) + 0.5*(lattice->get_f(ind_iq + 3) - lattice->get_f(ind_iq + 1))
                          + rho_tmp*ux_bc_/2.0 );

    lattice->set_f(ind_iq + 2, lattice->get_f(ind_iq + 4) );
  }
};


void ZouHeBC::setZouHeVelBC2D_yp (int iy_, int ix0_, int ix1_, double ux_bc_) {
  for (int i = ix0_; i <= ix1_; ++i){
    int ind_iq = i * lattice->ny * lattice->q + iy_ * lattice->q; // + iq;

    double rho_tmp =  lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 1) +
                              2.0*( lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 5) );

    lattice->set_f(ind_iq + 7, lattice->get_f(ind_iq + 5) + 0.5*(lattice->get_f(ind_iq + 1) - lattice->get_f(ind_iq + 3))
                          - rho_tmp*ux_bc_/2.0 );

    lattice->set_f(ind_iq + 8, lattice->get_f(ind_iq + 6) + 0.5*(lattice->get_f(ind_iq + 3) - lattice->get_f(ind_iq + 1))
                          + rho_tmp*ux_bc_/2.0 );

    lattice->set_f(ind_iq + 4, lattice->get_f(ind_iq + 2) );
  }
};



void ZouHeBC::setZouHeVelBC3D_xn (int ix_, int iy0_, int iy1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_) {
  for (int j = iy0_; j <= iy1_; ++j){
    for (int k = iz0_; k <= iz1_; ++k){

      int ind_iq = ix_ * lattice->ny * lattice->nz * lattice->q + j * lattice->nz * lattice->q + k * lattice->q;

      double rho_tmp =  1.0/(1.0-ux_bc_)*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 4) +
                                          lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 11) +
                                          lattice->get_f(ind_iq + 17) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 12) +
                                2.0*( lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 14) + lattice->get_f(ind_iq + 8) +
                                      lattice->get_f(ind_iq + 16) + lattice->get_f(ind_iq + 10) ) );

      double Nyx = 0.5*(lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 11) + lattice->get_f(ind_iq + 17)
                        -(lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 12)))
                   - 1.0/3.0*rho_tmp*uy_bc_;

      double Nzx = 0.5*(lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 11)
                        -(lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 17) + lattice->get_f(ind_iq + 12)))
                   - 1.0/3.0*rho_tmp*uz_bc_;


      lattice->set_f(ind_iq + 1, lattice->get_f(ind_iq + 2) + 1.0/3.0*rho_tmp*ux_bc_ );

      lattice->set_f(ind_iq + 13, lattice->get_f(ind_iq + 14) + rho_tmp/6.0*(ux_bc_ - uy_bc_) - Nyx );

      lattice->set_f(ind_iq + 7, lattice->get_f(ind_iq + 8) + rho_tmp/6.0*(ux_bc_ + uy_bc_) - Nyx );

      lattice->set_f(ind_iq + 9, lattice->get_f(ind_iq + 10) + rho_tmp/6.0*(ux_bc_ + uz_bc_) - Nzx );

      lattice->set_f(ind_iq + 15, lattice->get_f(ind_iq + 16) + rho_tmp/6.0*(ux_bc_ - uz_bc_) + Nzx );

    }
  }
};



void ZouHeBC::setZouHeVelBC3D_xp (int ix_, int iy0_, int iy1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_) {
  for (int j = iy0_; j <= iy1_; ++j){
    for (int k = iz0_; k <= iz1_; ++k){

      int ind_iq = ix_ * lattice->ny * lattice->nz * lattice->q + j * lattice->nz * lattice->q + k * lattice->q;

      double rho_tmp =  1.0/(1.0+ux_bc_)*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 4) +
                                          lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 11) +
                                          lattice->get_f(ind_iq + 17) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 12) +
                                2.0*( lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 7) + lattice->get_f(ind_iq + 13) +
                                      lattice->get_f(ind_iq + 9) + lattice->get_f(ind_iq + 15) ) );

      double Nyx = 0.5*(lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 11) + lattice->get_f(ind_iq + 17)
                        -(lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 12)))
                   - 1.0/3.0*rho_tmp*uy_bc_;

      double Nzx = 0.5*(lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 11)
                        -(lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 17) + lattice->get_f(ind_iq + 12)))
                   - 1.0/3.0*rho_tmp*uz_bc_;


      lattice->set_f(ind_iq + 2, lattice->get_f(ind_iq + 1) - 1.0/3.0*rho_tmp*ux_bc_ );

      lattice->set_f(ind_iq + 14, lattice->get_f(ind_iq + 113) + rho_tmp/6.0*(-ux_bc_ + uy_bc_) - Nyx );

      lattice->set_f(ind_iq + 8, lattice->get_f(ind_iq + 7) + rho_tmp/6.0*(-ux_bc_ - uy_bc_) + Nyx );

      lattice->set_f(ind_iq + 10, lattice->get_f(ind_iq + 9) + rho_tmp/6.0*(-ux_bc_ - uz_bc_) + Nzx );

      lattice->set_f(ind_iq + 16, lattice->get_f(ind_iq + 15) + rho_tmp/6.0*(-ux_bc_ + uz_bc_) - Nzx );

    }
  }
};



void ZouHeBC::setZouHeDensBC3D_xp (int ix_, int iy0_, int iy1_, int iz0_, int iz1_, double rho_bc_, double uy_bc_, double uz_bc_) {
  for (int j = iy0_; j <= iy1_; ++j){
    for (int k = iz0_; k <= iz1_; ++k){

      int ind_iq = ix_ * lattice->ny * lattice->nz * lattice->q + j * lattice->nz * lattice->q + k * lattice->q;

      double ux_bc_ =  -1.0 + 1.0/rho_bc_*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 4) +
                                          lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 11) +
                                          lattice->get_f(ind_iq + 17) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 12) +
                                2.0*( lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 7) + lattice->get_f(ind_iq + 13) +
                                      lattice->get_f(ind_iq + 9) + lattice->get_f(ind_iq + 15) ) );

      double Nyx = 0.5*(lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 11) + lattice->get_f(ind_iq + 17)
                        -(lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 12)))
                   - 1.0/3.0*rho_bc_*uy_bc_;

      double Nzx = 0.5*(lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 11)
                        -(lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 17) + lattice->get_f(ind_iq + 12)))
                   - 1.0/3.0*rho_bc_*uz_bc_;


      lattice->set_f(ind_iq + 2, lattice->get_f(ind_iq + 1) - 1.0/3.0*rho_bc_*ux_bc_ );

      lattice->set_f(ind_iq + 14, lattice->get_f(ind_iq + 113) + rho_bc_/6.0*(-ux_bc_ + uy_bc_) - Nyx );

      lattice->set_f(ind_iq + 8, lattice->get_f(ind_iq + 7) + rho_bc_/6.0*(-ux_bc_ - uy_bc_) + Nyx );

      lattice->set_f(ind_iq + 10, lattice->get_f(ind_iq + 9) + rho_bc_/6.0*(-ux_bc_ - uz_bc_) + Nzx );

      lattice->set_f(ind_iq + 16, lattice->get_f(ind_iq + 15) + rho_bc_/6.0*(-ux_bc_ + uz_bc_) - Nzx );

    }
  }
};



void ZouHeBC::setZouHeVelBC3D_yn (int iy_, int ix0_, int ix1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_) {
  for (int i = ix0_; i <= ix1_; ++i){
    for (int k = iz0_; k <= iz1_; ++k){

      int ind_iq = i * lattice->ny * lattice->nz * lattice->q + iy_ * lattice->nz * lattice->q + k * lattice->q; // + iq;

      double rho_tmp =  1.0/(1.0-uy_bc_)*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 2) +
                                          lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 9) +
                                          lattice->get_f(ind_iq + 15) + lattice->get_f(ind_iq + 16) + lattice->get_f(ind_iq + 10) +
                                2.0*( lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 13) + lattice->get_f(ind_iq + 8) +
                                      lattice->get_f(ind_iq + 18) + lattice->get_f(ind_iq + 12) ) );

      double Nxy = 0.5*(lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 9) + lattice->get_f(ind_iq + 15)
                        -(lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 16) + lattice->get_f(ind_iq + 10)))
                   - 1.0/3.0*rho_tmp*ux_bc_;

      double Nzy = 0.5*(lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 9) + lattice->get_f(ind_iq + 16)
                        -(lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 15) + lattice->get_f(ind_iq + 10)))
                   - 1.0/3.0*rho_tmp*uz_bc_;


      lattice->set_f(ind_iq + 3, lattice->get_f(ind_iq + 4) + 1.0/3.0*rho_tmp*uy_bc_ );

      lattice->set_f(ind_iq + 7, lattice->get_f(ind_iq + 8) + rho_tmp/6.0*(uy_bc_ + ux_bc_) - Nxy );

      lattice->set_f(ind_iq + 14, lattice->get_f(ind_iq + 13) + rho_tmp/6.0*(uy_bc_ - ux_bc_) + Nxy );

      lattice->set_f(ind_iq + 11, lattice->get_f(ind_iq + 12) + rho_tmp/6.0*(uy_bc_ + uz_bc_) - Nzy );

      lattice->set_f(ind_iq + 17, lattice->get_f(ind_iq + 18) + rho_tmp/6.0*(uy_bc_ - uz_bc_) + Nzy );

    }
  }
};


void ZouHeBC::setZouHeVelBC3D_yp (int iy_, int ix0_, int ix1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_) {
  for (int i = ix0_; i <= ix1_; ++i){
    for (int k = iz0_; k <= iz1_; ++k){

      int ind_iq = i * lattice->ny * lattice->nz * lattice->q + iy_ * lattice->nz * lattice->q + k * lattice->q; // + iq;

      double rho_tmp =  1.0/(1.0+uy_bc_)*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 2) +
                                          lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 9) +
                                          lattice->get_f(ind_iq + 15) + lattice->get_f(ind_iq + 16) + lattice->get_f(ind_iq + 10) +
                                2.0*( lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 7) + lattice->get_f(ind_iq + 14) +
                                      lattice->get_f(ind_iq + 11) + lattice->get_f(ind_iq + 17) ) );

      double Nxy = 0.5*(lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 9) + lattice->get_f(ind_iq + 15)
                        -(lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 16) + lattice->get_f(ind_iq + 10)))
                   - 1.0/3.0*rho_tmp*ux_bc_;

      double Nzy = 0.5*(lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 9) + lattice->get_f(ind_iq + 16)
                        -(lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 15) + lattice->get_f(ind_iq + 10)))
                   - 1.0/3.0*rho_tmp*uz_bc_;


      lattice->set_f(ind_iq + 4, lattice->get_f(ind_iq + 3) - 1.0/3.0*rho_tmp*uy_bc_ );

      lattice->set_f(ind_iq + 8, lattice->get_f(ind_iq + 7) + rho_tmp/6.0*(-uy_bc_ - ux_bc_) + Nxy );

      lattice->set_f(ind_iq + 13, lattice->get_f(ind_iq + 14) + rho_tmp/6.0*(-uy_bc_ + ux_bc_) - Nxy );

      lattice->set_f(ind_iq + 12, lattice->get_f(ind_iq + 11) + rho_tmp/6.0*(-uy_bc_ - uz_bc_) + Nzy );

      lattice->set_f(ind_iq + 18, lattice->get_f(ind_iq + 17) + rho_tmp/6.0*(-uy_bc_ + uz_bc_) - Nzy );

    }
  }
};


void ZouHeBC::setZouHeVelBC3D_zn (int iz_, int ix0_, int ix1_, int iy0_, int iy1_, double ux_bc_, double uy_bc_, double uz_bc_) {
  for (int i = ix0_; i <= ix1_; ++i){
    for (int j = iy0_; j <= iy1_; ++j){

      int ind_iq = i * lattice->ny * lattice->nz * lattice->q + j * lattice->nz * lattice->q + iz_ * lattice->q; // + iq;

      double rho_tmp =  1.0/(1.0-uz_bc_)*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 2) +
                                          lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 7) +
                                          lattice->get_f(ind_iq + 13) + lattice->get_f(ind_iq + 14) + lattice->get_f(ind_iq + 8) +
                                2.0*( lattice->get_f(ind_iq + 6) + lattice->get_f(ind_iq + 15) + lattice->get_f(ind_iq + 10) +
                                      lattice->get_f(ind_iq + 17) + lattice->get_f(ind_iq + 12) ) );

      double Nxz = 0.5*(lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 7) + lattice->get_f(ind_iq + 13)
                        -(lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 14) + lattice->get_f(ind_iq + 8)))
                   - 1.0/3.0*rho_tmp*ux_bc_;

      double Nyz = 0.5*(lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 7) + lattice->get_f(ind_iq + 14)
                        -(lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 13) + lattice->get_f(ind_iq + 8)))
                   - 1.0/3.0*rho_tmp*uy_bc_;


      lattice->set_f(ind_iq + 5, lattice->get_f(ind_iq + 6) + 1.0/3.0*rho_tmp*uz_bc_ );

      lattice->set_f(ind_iq + 9, lattice->get_f(ind_iq + 10) + rho_tmp/6.0*(uz_bc_ + ux_bc_) - Nxz );

      lattice->set_f(ind_iq + 16, lattice->get_f(ind_iq + 15) + rho_tmp/6.0*(uz_bc_ - ux_bc_) + Nxz );

      lattice->set_f(ind_iq + 11, lattice->get_f(ind_iq + 12) + rho_tmp/6.0*(uz_bc_ + uy_bc_) - Nyz );

      lattice->set_f(ind_iq + 18, lattice->get_f(ind_iq + 17) + rho_tmp/6.0*(uz_bc_ - uy_bc_) + Nyz );

    }
  }
};


void ZouHeBC::setZouHeVelBC3D_zp (int iz_, int ix0_, int ix1_, int iy0_, int iy1_, double ux_bc_, double uy_bc_, double uz_bc_) {
  for (int i = ix0_; i <= ix1_; ++i){
    for (int j = iy0_; j <= iy1_; ++j){

      int ind_iq = i * lattice->ny * lattice->nz * lattice->q + j * lattice->nz * lattice->q + iz_ * lattice->q; // + iq;

      double rho_tmp =  1.0/(1.0+uz_bc_)*(lattice->get_f(ind_iq + 0) + lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 2) +
                                          lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 7) +
                                          lattice->get_f(ind_iq + 14) + lattice->get_f(ind_iq + 8) + lattice->get_f(ind_iq + 13) +
                                2.0*( lattice->get_f(ind_iq + 5) + lattice->get_f(ind_iq + 9) + lattice->get_f(ind_iq + 16) +
                                      lattice->get_f(ind_iq + 11) + lattice->get_f(ind_iq + 18) ) );

      double Nxz = 0.5*(lattice->get_f(ind_iq + 1) + lattice->get_f(ind_iq + 7) + lattice->get_f(ind_iq + 13)
                        -(lattice->get_f(ind_iq + 2) + lattice->get_f(ind_iq + 14) + lattice->get_f(ind_iq + 8)))
                   - 1.0/3.0*rho_tmp*ux_bc_;

      double Nyz = 0.5*(lattice->get_f(ind_iq + 3) + lattice->get_f(ind_iq + 7) + lattice->get_f(ind_iq + 14)
                        -(lattice->get_f(ind_iq + 4) + lattice->get_f(ind_iq + 13) + lattice->get_f(ind_iq + 8)))
                   - 1.0/3.0*rho_tmp*uy_bc_;


      lattice->set_f(ind_iq + 6, lattice->get_f(ind_iq + 5) - 1.0/3.0*rho_tmp*uz_bc_ );

      lattice->set_f(ind_iq + 10, lattice->get_f(ind_iq + 16) + rho_tmp/6.0*(-uz_bc_ + ux_bc_) - Nxz );

      lattice->set_f(ind_iq + 10, lattice->get_f(ind_iq + 9) + rho_tmp/6.0*(-uz_bc_ - ux_bc_) + Nxz );

      lattice->set_f(ind_iq + 17, lattice->get_f(ind_iq + 18) + rho_tmp/6.0*(-uz_bc_ + uy_bc_) - Nyz );

      lattice->set_f(ind_iq + 12, lattice->get_f(ind_iq + 11) + rho_tmp/6.0*(-uz_bc_ - uy_bc_) + Nyz );

    }
  }
};
