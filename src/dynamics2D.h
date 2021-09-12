/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#ifndef DYNAMICS_2D
#define DYNAMICS_2D

#include <iostream>
#include <math.h>
#include <vector>

#include "lattice2D.h"

using namespace std;

class Dynamics2D : public Lattice2D {
  
  private:
  //  double tau;

  public:
    Dynamics2D(int nx_, int ny_, int q_, int decomposition_[3], int procCoordinates_[3]);
    ~Dynamics2D();

    double feq(int iq_, int ind_phys_1D_, int ind_phys_2D_, vector<double> &rho_, vector<double> &u_);
    double feq(int iq_, int ind_phys_1D_, int ind_phys_2D_, double rho, vector<double> u);

    virtual void compute_macro_values() = 0;
    void streaming_periodic_2D();
    void streaming();
    void streamBulk(int i_, int j_, int iq_);
    void streamBC_xn(int i_, int j_, int iq_, int corner_);
    void streamBC_xp(int i_, int j_, int iq_, int corner_);
    void streamBC_yn(int i_, int j_, int iq_, int corner_);
    void streamBC_yp(int i_, int j_, int iq_, int corner_);
    virtual void collision() = 0;

};

#endif
