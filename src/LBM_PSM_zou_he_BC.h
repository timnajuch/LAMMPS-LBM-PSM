/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
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
    double oneOverThree;
    double twoOverThree;
    double oneOverSix;
    
  public:
    ZouHeBC(LBMPSMLattice *lattice_);
    ~ZouHeBC();

    void setZouHeVelBC2D_xn(int ix_, int iy0_, int iy1_, double ux_bc_, int currentStep);
    void setZouHeVelBC2D_xp(int ix_, int iy0_, int iy1_, double ux_bc_, int currentStep);
    void setZouHeDensBC2D_xp(int ix_, int iy0_, int iy1_, double rho_bc_, int currentStep);
    void setZouHeVelBC2D_yn(int iy_, int ix0_, int ix1_, double ux_bc_, int currentStep);
    void setZouHeVelBC2D_yp(int iy_, int ix0_, int ix1_, double ux_bc_, int currentStep);
    void setZouHeNeumannVelBC2D_yn (int iy_, int ix0_, int ix1_, int currentStep);
    void setZouHeNeumannVelBC2D_yp (int iy_, int ix0_, int ix1_, int currentStep);


    void setZouHeVelBC3D_xn(int ix_, int iy0_, int iy1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_, int currentStep);
    void setZouHeVelBC3D_xp(int ix_, int iy0_, int iy1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_, int currentStep);
    void setZouHeDensBC3D_xp(int ix_, int iy0_, int iy1_, int iz0_, int iz1_, double rho_bc_, double uy_bc_, double uz_bc_, int currentStep);
    void setZouHeVelBC3D_yn(int iy_, int ix0_, int ix1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_, int currentStep);
    void setZouHeVelBC3D_yp(int iy_, int ix0_, int ix1_, int iz0_, int iz1_, double ux_bc_, double uy_bc_, double uz_bc_, int currentStep);
    void setZouHeVelBC3D_zn(int iz_, int ix0_, int ix1_, int iy0_, int iy1_, double ux_bc_, double uy_bc_, double uz_bc_, int currentStep);
    void setZouHeVelBC3D_zp(int iz_, int ix0_, int ix1_, int iy0_, int iy1_, double ux_bc_, double uy_bc_, double uz_bc_, int currentStep);

};

#endif
