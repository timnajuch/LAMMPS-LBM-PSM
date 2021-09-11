/* -*- c++ -*- -----------------------------------------------------------
   LAMMPS 2003 (July 31) - Molecular Dynamics Simulator
   Sandia National Laboratories, www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   For more info, see the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------ */

#ifdef FIX_CLASS

FixStyle(lbm-psm,fix_PSM_LBM)

#else

#ifndef LMP_FIX_PSM_LBM_FLUID_H
#define LMP_FIX_PSM_LBM_FLUID_H


#include <cmath>
#include <cstring>
#include <algorithm>
#include <utility>
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "atom.h"
#include "group.h"
#include "random_mars.h"
#include "update.h"
#include "force.h"
#include "modify.h"

#include "BGK_GuoExtForce_dynamics2D.h"
//#include "mpiCommunication.h"
//#include "fix_PSM_LBM_MPICOMM.h"
#include "PSM_LBM_MPICOMM.h"
#include "fix_PSM_LBM_BC.h"
//#include "fix_exchangeParticleData.h"
#include "ExchangeParticleData.h"
#include "unit_conversion.h"

#include "fix.h"
#include "comm.h"

//#if defined(MPI_STUBS)
//#error "The USER-LB package cannot be compiled in serial with MPI STUBS"
//#endif

//namespace LAMMPS_NS {

using namespace LAMMPS_NS;
using namespace FixConst;

class fix_PSM_LBM : public Fix {
  //friend class BGK_GuoExtForce_Dynamics2D;
  //friend class fix_PSM_LBM_MPI;
 // friend class fix_PSM_LBM_BC;
  //friend class fix_ExchangeParticleDate;
//  friend class WriteVTK;

public:
  fix_PSM_LBM(class LAMMPS *, int, char **);
  ~fix_PSM_LBM();
  int setmask();
  void init();
//    void initial_integrate(int);
//    void setup(int);
//    void post_force(int);
//    void end_of_step();
  void pre_force(int);

  BGK_GuoExtForce_Dynamics2D *dynamics;
//  BGK_GuoExtForce_Dynamics2D dynamics;
  int get_nx();
  int get_ny();
  vector<double> get_x();
  vector<double> get_y();
  vector<double> get_rho();
  vector<double> get_u();
  vector<double> get_B();

  Unit_Conversion *unitConversion;
  ExchangeParticleData *exchangeParticleData;
  PSM_LBM_MPI *lbmmpicomm;

};

//}
#endif
#endif

