/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details. 

Tim Najuch, 2022
------------------------------------------------------*/

#ifndef LBMPSM_BGK_DYNAMICS
#define LBMPSM_BGK_DYNAMICS

#include "LBM_PSM_lattice.h"
#include "LBM_PSM_dynamics.h"
#include "domain.h"

class LBMPSMBGKDynamics : public LBMPSMDynamics {
  
  private:
    double tau, omega;      // Single relaxation parameter tau and its inverse omega
    vector<double> F_lbm;
    double F_lbm_mag_pow2;

  public:
    LBMPSMBGKDynamics(double tau_, int nx_, int ny_, int nz_, vector<double> F_lbm_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_);
    ~LBMPSMBGKDynamics();

    void initialise_dynamics(double rho_, double ux_, double uy_, double uz_);
    void compute_macro_values(int i_, int j_, int k_, int currentStep_);
    inline void collisionAndStream(int i_, int j_, int k_, int iq_, int iShift_, int jShift_, int kShift_, int currentStep_, int nextStep_);
    void macroCollideStream();
    
};



inline void LBMPSMBGKDynamics::collisionAndStream(int i_, int j_, int k_, int iq_, int iShift_, int jShift_, int kShift_, int currentStep_, int nextStep_){
  int ind_iq = index_fi(i_, j_, k_, iq_, currentStep_); // Index for populations
  int ind_iq0 = index_fi(i_, j_, k_, iq_, 0);           // Index for equilibrium populations
  int ind_phys_1D = index_1D(i_, j_, k_);
  int ind_phys_2D = index_2D(i_, j_, k_, 0);

  double rho = get_rho(ind_phys_1D);
  double fi = f[ind_iq];
  double fi_eq = feq(iq_, ind_phys_1D, ind_phys_2D);;
  f0[ind_iq0] = fi_eq;

  auto& nodeData = getReferenceParticleDataOnLatticeNode(ind_phys_1D);

  // Solid phase / particle data
  LAMMPS_NS::tagint pID1 = nodeData.particleID[0];
  LAMMPS_NS::tagint pID2 = nodeData.particleID[1];

  // External forcing according to Guo et al. (2002)
  double F_lbm_iq = 0.0;
  if (F_lbm_mag_pow2 > 0.0)
    { F_lbm_iq = (1.0-0.5*omega) * F_iq(iq_, ind_phys_2D, F_lbm); }


  // Collision and streaming for pure fluid node
  if (pID1 <= 0 && pID2 <= 0){    
    set_f(iShift_, jShift_, kShift_, iq_, nextStep_, fi  + (fi_eq - fi)*omega  + F_lbm_iq);
    return;
  }

  // Solid phase on node
  double B1 = 0.0;
  double B2 = 0.0;
  double solid_coll1 = 0.0;
  double solid_coll2 = 0.0;
  double oneMinusOmega = 1.0 - omega;

  // Solid phase collision terms
  if (pID1 > 0){
    B1 = nodeData.solidFraction[0];
    std::array<double, 3> uSolid1{nodeData.particleVelocity[0], nodeData.particleVelocity[1], nodeData.particleVelocity[2]};
    double f0_solid1 = feq(iq_, rho, uSolid1);
    solid_coll1 = f0_solid1 - fi + oneMinusOmega * (fi - fi_eq);
  }
  if (pID2 > 0){
    B2 = nodeData.solidFraction[1];
    std::array<double, 3> uSolid2{nodeData.particleVelocity[3], nodeData.particleVelocity[4], nodeData.particleVelocity[5]};
    double f0_solid2 = feq(iq_, rho, uSolid2);
    solid_coll2 = f0_solid2 - fi + oneMinusOmega * (fi - fi_eq);
  }

  double Btot = B1 + B2;
  if(Btot > 1.0){
    B1 = B1/Btot;
    B2 = B2/Btot;
    Btot = 1.0;
  }

  // Collision and streaming
  set_f(iShift_, jShift_, kShift_, iq_, nextStep_, fi  + (1.0 - Btot)*(fi_eq - fi)*omega  + B1*solid_coll1 + B2*solid_coll2 + F_lbm_iq);

  // Add hydrodynamic interaction force
  double ex = e[3*iq_];
  double ey = e[3*iq_+1];
  double ez = e[3*iq_+2];
  if (pID1 > 0){
    add_Fhyd(ind_phys_1D, pID1, -B1 * solid_coll1 * ex, 0);
    add_Fhyd(ind_phys_1D, pID1, -B1 * solid_coll1 * ey, 1);
    add_Fhyd(ind_phys_1D, pID1, -B1 * solid_coll1 * ez, 2);
  }
  if (pID2 > 0){
    add_Fhyd(ind_phys_1D, pID2, -B2 * solid_coll2 * ex, 0);
    add_Fhyd(ind_phys_1D, pID2, -B2 * solid_coll2 * ey, 1);
    add_Fhyd(ind_phys_1D, pID2, -B2 * solid_coll2 * ez, 2);
  }
}



#endif
