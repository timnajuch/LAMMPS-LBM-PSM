/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#ifndef ZOU_HE_BC_H
#define ZOU_HE_BC_H

#include <iostream>
#include <vector>
#include <math.h>
#include "LBM_PSM_lattice.h"

using namespace std;

class ZouHeBC {
  private:
    LBMPSMLattice *lattice;
    
  public:
    ZouHeBC(LBMPSMLattice *lattice_);
    ~ZouHeBC();
    
    void setZouHeVelBC2D_xn(int ix_, int iy0_, int iy1_, double ux_bc_);
    void setZouHeDensBC2D_xp(int ix_, int iy0_, int iy1_, double rho_bc_);
    void setZouHeVelBC2D_yn(int iy_, int ix0_, int ix1_, double ux_bc_);
    void setZouHeVelBC2D_yp(int iy_, int ix0_, int ix1_, double ux_bc_);


    void setZouHeVelBC3D_xn(int ix_, int iy0_, int iy1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_);
    void setZouHeVelBC3D_xp(int ix_, int iy0_, int iy1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_);
    void setZouHeDensBC3D_xp(int ix_, int iy0_, int iy1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_);
    void setZouHeVelBC3D_yn(int iy_, int ix0_, int ix1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_);
    void setZouHeVelBC3D_yp(int iy_, int ix0_, int ix1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_);
    void setZouHeVelBC3D_zn(int iz_, int ix0_, int ix1_, int iy0_, int iy1_, double ux_bc_, double uy_bc_, double uz_bc_);
    void setZouHeVelBC3D_zp(int iz_, int ix0_, int ix1_, int iy0_, int iy1_, double ux_bc_, double uy_bc_, double uz_bc_);

};

#endif
