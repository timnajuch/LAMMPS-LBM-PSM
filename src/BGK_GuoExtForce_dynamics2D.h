/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#ifndef BGK_GUOEXTFORCE_DYNAMICS_2D
#define BGK_GUOEXTFORCE_DYNAMICS_2D

#include "lattice2D.h"
#include "dynamics2D.h"
#include "domain.h"

class BGK_GuoExtForce_Dynamics2D : public Dynamics2D {
  
  private:
    double tau;
    vector<double> F_lbm;
    int dimension;

  public:
    BGK_GuoExtForce_Dynamics2D(double tau_, int nx_, int ny_, int nz_, int q_, vector<double> F_lbm_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_, double dx_);
    ~BGK_GuoExtForce_Dynamics2D();

    void compute_macro_values();
 //   void collision();

    void initialise_dynamics(double rho_, double ux_, double uy_, double uz_);

    void macroCollideStream();
    void compute_macro_values(int i_, int j_, int k_);
    void collision(int i_, int j_, int k_, int iq_);
    
};

#endif
