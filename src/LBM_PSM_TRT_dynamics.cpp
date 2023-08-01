/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#include "LBM_PSM_TRT_dynamics.h"

LBMPSMTRTDynamics::LBMPSMTRTDynamics(double tau_, double tau_m_, int nx_, int ny_, int nz_, vector<double> F_lbm_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_) :
  LBMPSMDynamics(nx_, ny_, nz_, decomposition_, procCoordinates_, origin_, boxLength_, dimension_), tau_p(tau_), tau_m(tau_m_), omega_p(1.0/tau_), omega_m(1.0/tau_m_), F_lbm(F_lbm_)
  {
    F_hyd_1.resize(3, 0.0);
    F_hyd_2.resize(3, 0.0);
  }


LBMPSMTRTDynamics::~LBMPSMTRTDynamics(){}


void  LBMPSMTRTDynamics::initialise_dynamics(double rho_, double ux_, double uy_, double uz_){
  for(int i = 0; i < nx; ++i){
    for(int j = 0; j < ny; ++j){
      for(int k = 0; k < nz; ++k){
        int ind_phys_1D = index_1D(i, j, k);
        int ind_phys_2D = index_2D(i, j, k, 0);
        rho[ind_phys_1D] = rho_;
        u[ind_phys_2D] = ux_;
        u[ind_phys_2D+1] = uy_;
        u[ind_phys_2D+2] = uz_;

        for(int iq = 0; iq < q; ++iq){
            set_f0(i, j, k, iq, feq(iq, ind_phys_1D, ind_phys_2D) );
            set_f(i, j, k, iq, 0, feq(iq, ind_phys_1D, ind_phys_2D) );
            set_f(i, j, k, iq, 1, feq(iq, ind_phys_1D, ind_phys_2D) );
        }
      }
    }
  }
  F_lbm_mag_pow2 = F_lbm[0]*F_lbm[0] + F_lbm[1]*F_lbm[1] + F_lbm[2]*F_lbm[2];
}



void LBMPSMTRTDynamics::compute_macro_values(int i_, int j_, int k_, int currentStep_){
  int ind_phys_1D = index_1D(i_, j_, k_);
  int ind_phys_2D = index_2D(i_, j_, k_, 0);

  double rho_tmp = 0.0;
  double jx = 0.0;
  double jy = 0.0;
  double jz = 0.0;
  
  int ind_iq = 0;
  for(int iq = 0; iq < q; ++iq){
    ind_iq = index_fi(i_, j_, k_, iq, currentStep_);
    rho_tmp += f[ind_iq];
    jx += f[ind_iq] * e[3*iq];
    jy += f[ind_iq] * e[3*iq+1];
    jz += f[ind_iq] * e[3*iq+2];
  }
  rho[ind_phys_1D] = rho_tmp;

  double halfPI = 3.14159265358979323846/2.0;
  // Determine and add hydrodynamic interaction force
  LAMMPS_NS::tagint pID1 = getParticleDataOnLatticeNode(ind_phys_1D).particleID[0];
  LAMMPS_NS::tagint pID2 = getParticleDataOnLatticeNode(ind_phys_1D).particleID[1];
  F_hyd_1[0] = F_hyd_1[1] = F_hyd_1[2] = 0.0;
  F_hyd_2[0] = F_hyd_2[1] = F_hyd_2[2] = 0.0;
  if (pID1 > 0){
    //double B1 = getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[0];
    //double B1 = sin(getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[0]*halfPI);
    double B1 = pow(getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[0], 0.5);
    vector<double> uSolid1{0.0, 0.0, 0.0};
    uSolid1[0] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[0];
    uSolid1[1] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[1];
    uSolid1[2] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[2];

    // Velocity corrections
    vector<double> dU{0.0, 0.0, 0.0};
    dU[0] = (uSolid1[0] - jx/rho_tmp);// /B1*2.0;
    dU[1] = (uSolid1[1] - jy/rho_tmp);// /B1*2.0;
    dU[2] = (uSolid1[2] - jz/rho_tmp);// /B1*2.0;

    // Hydrodynamic interaction force based on corrections
    F_hyd_1[0] = dU[0]*B1;
    F_hyd_1[1] = dU[1]*B1;
    F_hyd_1[2] = dU[2]*B1;

    //F_hyd_1[0] = rho[ind_phys_1D]*(uSolid1[0] - u[ind_phys_1D][0]);
    //F_hyd_1[1] = rho[ind_phys_1D]*(uSolid1[1] - u[ind_phys_1D][1]);
    //F_hyd_1[2] = rho[ind_phys_1D]*(uSolid1[2] - u[ind_phys_1D][2]);
  //F_hyd_1_iq = (1.0-0.5*omega) * F_iq(iq_, ind_phys_2D, F_hyd_1);
    add_Fhyd(ind_phys_1D, pID1, -rho_tmp*F_hyd_1[0], 0);
    add_Fhyd(ind_phys_1D, pID1, -rho_tmp*F_hyd_1[1], 1);
    add_Fhyd(ind_phys_1D, pID1, -rho_tmp*F_hyd_1[2], 2);
  }

  if (pID2 > 0){
    //double B2 = getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[1];
    //double B2 = sin(getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[1]*halfPI);
    double B2 = pow(getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[1], 0.5);
    vector<double> uSolid2{0.0, 0.0, 0.0};
    uSolid2[0] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[3];
    uSolid2[1] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[4];
    uSolid2[2] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[5];

    // Velocity corrections
    vector<double> dU{0.0, 0.0, 0.0};
    dU[0] = (uSolid2[0] - jx/rho_tmp);// /B2*2.0;
    dU[1] = (uSolid2[1] - jy/rho_tmp);// /B2*2.0;
    dU[2] = (uSolid2[2] - jz/rho_tmp);// /B2*2.0;

    // Hydrodynamic interaction force based on corrections
    F_hyd_2[0] = dU[0]*B2;
    F_hyd_2[1] = dU[1]*B2;
    F_hyd_2[2] = dU[2]*B2;

    //F_hyd_2[0] = rho[ind_phys_1D]*(uSolid2[0] - u[ind_phys_1D][0]);
    //F_hyd_2[1] = rho[ind_phys_1D]*(uSolid2[1] - u[ind_phys_1D][1]);
    //F_hyd_2[2] = rho[ind_phys_1D]*(uSolid2[2] - u[ind_phys_1D][2]);
    //F_hyd_2_iq = (1.0-0.5*omega) * F_iq(iq_, ind_phys_2D, F_hyd_2);
    add_Fhyd(ind_phys_1D, pID2, -rho_tmp*F_hyd_2[0], 0);
    add_Fhyd(ind_phys_1D, pID2, -rho_tmp*F_hyd_2[1], 1);
    add_Fhyd(ind_phys_1D, pID2, -rho_tmp*F_hyd_2[2], 2);
  }

  // External force according to Guo et al. (2002). Hydrodynamic force added analogously
  //jx += (F_lbm[0] + F_hyd_1[0] + F_hyd_2[0])*0.5;
  //jy += (F_lbm[1] + F_hyd_1[1] + F_hyd_2[1])*0.5;
  //jz += (F_lbm[2] + F_hyd_1[2] + F_hyd_2[2])*0.5;
  jx += F_lbm[0]*0.5;
  jy += F_lbm[1]*0.5;
  jz += F_lbm[2]*0.5;

  u[ind_phys_2D] = jx/rho_tmp + (F_hyd_1[0] + F_hyd_2[0])*0.5;
  u[ind_phys_2D+1] = jy/rho_tmp + (F_hyd_1[1] + F_hyd_2[1])*0.5;
  u[ind_phys_2D+2] = jz/rho_tmp + (F_hyd_1[2] + F_hyd_2[2])*0.5;
/*
  u[ind_phys_2D] = jx/rho_tmp;
  u[ind_phys_2D+1] = jy/rho_tmp;
  u[ind_phys_2D+2] = jz/rho_tmp;
*/
}



void LBMPSMTRTDynamics::collisionAndStream(int i_, int j_, int k_, int iq_, int iShift_, int jShift_, int kShift_, int currentStep_, int nextStep_){
  int ind_iq = index_fi(i_, j_, k_, iq_, currentStep_); // Index for populations
  int ind_iq0 = index_fi(i_, j_, k_, iq_, 0);           // Index for equilibrium populations
  int ind_phys_1D = index_1D(i_, j_, k_);
  int ind_phys_2D = index_2D(i_, j_, k_, 0);

  // Equilibrium function
  f0[ind_iq0] = feq(iq_, ind_phys_1D, ind_phys_2D);

  // Solid phase / particle data
  LAMMPS_NS::tagint pID1 = getParticleDataOnLatticeNode(ind_phys_1D).particleID[0];
  LAMMPS_NS::tagint pID2 = getParticleDataOnLatticeNode(ind_phys_1D).particleID[1];
  double B1 = 0.0;
  double B2 = 0.0;
  vector<double> uSolid1{0.0, 0.0, 0.0};
  vector<double> uSolid2{0.0, 0.0, 0.0};
  double f0_solid1 = 0.0;
  double f0_solid2 = 0.0;
  double solid_coll1 = 0.0;
  double solid_coll2 = 0.0;
  double F_hyd_1_iq = 0.0;
  double F_hyd_2_iq = 0.0;

  double F_hyd_1_iq_p = 0.0;
  double F_hyd_1_iq_m = 0.0;
  double F_hyd_2_iq_p = 0.0;
  double F_hyd_2_iq_m = 0.0;

  vector<double> Fhydro1{0.0, 0.0, 0.0};
  vector<double> Fhydro2{0.0, 0.0, 0.0};

  double eps0 = 0.2;
  double halfPI = 3.14159265358979323846/2.0;
  // Solid phase collision terms
  if (pID1 > 0){
    //B1 = pow(getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[0], 0.5);
    //B1 = eps0 + (1.0-eps0)*getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[0];
    B1 = sin(getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[0]*halfPI);
    uSolid1[0] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[0];
    uSolid1[1] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[1];
    uSolid1[2] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[2];
    f0_solid1 = feq(iq_, get_rho(ind_phys_1D), uSolid1);
    //solid_coll1 = f0_solid1 - f[ind_iq] + ( 1.0 - omega_p) * (f[ind_iq] - f0[ind_iq0]);
    solid_coll1 = f0_solid1 - f0[ind_iq0];
    //solid_coll1 = (1.0-0.5*omega) * F_iq(iq_, ind_phys_2D, F_lbm);
    /*F_hyd_1[0] = rho[ind_phys_1D]*(uSolid1[0] - u[ind_phys_1D][0]);
    F_hyd_1[1] = rho[ind_phys_1D]*(uSolid1[1] - u[ind_phys_1D][1]);
    F_hyd_1[2] = rho[ind_phys_1D]*(uSolid1[2] - u[ind_phys_1D][2]);*/
    //F_hyd_1_iq = (1.0-0.5*omega_p) * rho[ind_phys_1D] * F_iq(iq_, ind_phys_2D, F_hyd_1);
    //Fhydro1[0] = getParticleDataOnLatticeNode(ind_phys_1D).hydrodynamicForce[0];
    //Fhydro1[1] = getParticleDataOnLatticeNode(ind_phys_1D).hydrodynamicForce[1];
    //Fhydro1[2] = getParticleDataOnLatticeNode(ind_phys_1D).hydrodynamicForce[2];
  }
  if (pID2 > 0){
    //B2 = pow(getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[1], 0.5);
    //B2 = eps0 + (1.0-eps0)*getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[1];
    B2 = sin(getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[1]*halfPI);
    uSolid2[0] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[3];
    uSolid2[1] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[4];
    uSolid2[2] = getParticleDataOnLatticeNode(ind_phys_1D).particleVelocity[5];
    f0_solid2 = feq(iq_, get_rho(ind_phys_1D), uSolid2);
    //solid_coll2 = f0_solid2 - f[ind_iq] + ( 1.0 - omega_p) * (f[ind_iq] - f0[ind_iq0]);
    solid_coll2 = f0_solid2 - f0[ind_iq0];
    /*F_hyd_2[0] = rho[ind_phys_1D]*(uSolid2[0] - u[ind_phys_1D][0]);
    F_hyd_2[1] = rho[ind_phys_1D]*(uSolid2[1] - u[ind_phys_1D][1]);
    F_hyd_2[2] = rho[ind_phys_1D]*(uSolid2[2] - u[ind_phys_1D][2]);*/
    //F_hyd_2_iq = (1.0-0.5*omega_p) * rho[ind_phys_1D] * F_iq(iq_, ind_phys_2D, F_hyd_2);
    //Fhydro2[0] = getParticleDataOnLatticeNode(ind_phys_1D).hydrodynamicForce[3];
    //Fhydro2[1] = getParticleDataOnLatticeNode(ind_phys_1D).hydrodynamicForce[4];
    //Fhydro2[2] = getParticleDataOnLatticeNode(ind_phys_1D).hydrodynamicForce[5];
  }

  double Btot = B1 + B2;
  if(Btot > 1.0){
    B1 = B1/Btot;
    B2 = B2/Btot;
    Btot = 1.0;
  }

  int iq_m = 0;
  if(iq_ > 0)
  {
    if(dimension == 2){
      if(iq_ < 5)
      {
        iq_m = iq_ < 3 ? iq_ + 2 : iq_ - 2;
      }
      else
      {
        iq_m = iq_ < 7 ? iq_ + 2 : iq_ - 2;
      }
    } else { //3D
      if (iq_ % 2 == 0){
        iq_m = iq_ - 1;
      } else {
        iq_m = iq_ + 1;
      }
    }
  }
  int ind_iq_m = index_fi(i_, j_, k_, iq_m, currentStep_);
  int ind_iq0_m = index_fi(i_, j_, k_, iq_m, 0);           // Index for equilibrium populations

  // External forcing according to Guo et al. (2002)
  double F_lbm_iq = 0.0;
  if (F_lbm_mag_pow2 > 0.0)
    { F_lbm_iq = (1.0-0.5*omega_p) * F_iq(iq_, ind_phys_2D, F_lbm); }

  double f_p = 0.5*(f[ind_iq] + f[ind_iq_m]);
  double f_m = 0.5*(f[ind_iq] - f[ind_iq_m]);
  f0[ind_iq0_m] = feq(iq_m, ind_phys_1D, ind_phys_2D);
  double f0_p = 0.5*(f0[ind_iq0] + f0[ind_iq0_m]);
  double f0_m = 0.5*(f0[ind_iq0] - f0[ind_iq0_m]);

  F_hyd_1_iq_p = (1.0-0.5*omega_p) * rho[ind_phys_1D] * 0.5 * (F_iq(iq_, ind_phys_2D, F_hyd_1) + F_iq(iq_m, ind_phys_2D, F_hyd_1));
  F_hyd_1_iq_m = (1.0-0.5*omega_m) * rho[ind_phys_1D] * 0.5 * (F_iq(iq_, ind_phys_2D, F_hyd_1) - F_iq(iq_m, ind_phys_2D, F_hyd_1));
  F_hyd_2_iq_p = (1.0-0.5*omega_p) * rho[ind_phys_1D] * 0.5 * (F_iq(iq_, ind_phys_2D, F_hyd_2) + F_iq(iq_m, ind_phys_2D, F_hyd_2));
  F_hyd_2_iq_m = (1.0-0.5*omega_m) * rho[ind_phys_1D] * 0.5 * (F_iq(iq_, ind_phys_2D, F_hyd_2) - F_iq(iq_m, ind_phys_2D, F_hyd_2));

  // Collision and streaming
  set_f(iShift_, jShift_, kShift_, iq_, nextStep_, f[ind_iq] + (f0_p - f_p)*omega_p + (f0_m - f_m)*omega_m + F_hyd_1_iq_p + F_hyd_1_iq_m + F_hyd_2_iq_p + F_hyd_2_iq_m + F_lbm_iq);

/*
if (pID1 > 0){
    double f0_solid1_iqp = feq(iq_, get_rho(ind_phys_1D), uSolid1);
    double f0_solid1_iqm = feq(iq_m, get_rho(ind_phys_1D), uSolid1);
    double f0_solid1_p = 0.5*(f0_solid1_iqp + f0_solid1_iqm);
    double f0_solid1_m = 0.5*(f0_solid1_iqp - f0_solid1_iqm);
    solid_coll1 = omega_p*(f0_solid1_p - f0_p) + omega_m*(f0_solid1_m - f0_m);
}
if (pID2 > 0){
    double f0_solid2_iqp = feq(iq_, get_rho(ind_phys_1D), uSolid2);
    double f0_solid2_iqm = feq(iq_m, get_rho(ind_phys_1D), uSolid2);
    double f0_solid2_p = 0.5*(f0_solid2_iqp + f0_solid2_iqm);
    double f0_solid2_m = 0.5*(f0_solid2_iqp - f0_solid2_iqm);
    solid_coll2 = omega_p*(f0_solid2_p - f0_p) + omega_m*(f0_solid2_m - f0_m);
}
*/
  //set_f(iShift_, jShift_, kShift_, iq_, nextStep_, f[ind_iq]  + (1.0 - Btot)*( (f0_p - f_p)*omega_p + (f0_m - f_m)*omega_m )  + B1*solid_coll1 + B2*solid_coll2 + F_lbm_iq);
//  set_f(iShift_, jShift_, kShift_, iq_, nextStep_, f[ind_iq]  + (f0_p - f_p)*omega_p + (f0_m - f_m)*omega_m  + B1*solid_coll1 + B2*solid_coll2 + F_lbm_iq);

/*
  // Add hydrodynamic interaction force
  if (pID1 > 0){
    add_Fhyd(ind_phys_1D, pID1, -B1 * solid_coll1 * e[3*iq_], 0);
    add_Fhyd(ind_phys_1D, pID1, -B1 * solid_coll1 * e[3*iq_+1], 1);
    add_Fhyd(ind_phys_1D, pID1, -B1 * solid_coll1 * e[3*iq_+2], 2);
  }
  if (pID2 > 0){
    add_Fhyd(ind_phys_1D, pID2, -B2 * solid_coll2 * e[3*iq_], 0);
    add_Fhyd(ind_phys_1D, pID2, -B2 * solid_coll2 * e[3*iq_+1], 1);
    add_Fhyd(ind_phys_1D, pID2, -B2 * solid_coll2 * e[3*iq_+2], 2);
  }
*/
// Add hydrodynamic interaction force
/*
  if (pID1 > 0){
    add_Fhyd(ind_phys_1D, pID1, -(F_hyd_1_iq_p + F_hyd_1_iq_m) * e[3*iq_], 0);
    add_Fhyd(ind_phys_1D, pID1, -(F_hyd_1_iq_p + F_hyd_1_iq_m) * e[3*iq_+1], 1);
    add_Fhyd(ind_phys_1D, pID1, -(F_hyd_1_iq_p + F_hyd_1_iq_m) * e[3*iq_+2], 2);
  }
  if (pID2 > 0){
    add_Fhyd(ind_phys_1D, pID2, -(F_hyd_2_iq_p + F_hyd_2_iq_m) * e[3*iq_], 0);
    add_Fhyd(ind_phys_1D, pID2, -(F_hyd_2_iq_p + F_hyd_2_iq_m) * e[3*iq_+1], 1);
    add_Fhyd(ind_phys_1D, pID2, -(F_hyd_2_iq_p + F_hyd_2_iq_m) * e[3*iq_+2], 2);
  }
  */
}



void LBMPSMTRTDynamics::macroCollideStream(){
  int iShift = 0;
  int jShift = 0;
  int kShift = 0;

  for(int i = 0; i < nx; ++i){
    for(int j = 0; j < ny; ++j){
      for(int k = 0; k < nz; ++k){

        compute_macro_values(i, j, k, currentStep);

        for(int iq = 0; iq < q; ++iq){

          iShift = i + e[3*iq];
          if (iShift < 0 || iShift > nx-1) { continue; }
          jShift = j + e[3*iq+1];
          if (jShift < 0 || jShift > ny-1) { continue; }
          kShift = k + e[3*iq+2];
          if (kShift < 0 || kShift > nz-1) { continue; }

          collisionAndStream(i, j, k, iq, iShift, jShift, kShift, currentStep, nextStep);
        }

      }
    }
  }

  currentStep = nextStep;
  nextStep = 1 - currentStep;
}