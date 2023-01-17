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

FixStyle(lbm-psm,fix_LBM_PSM)

#else

#ifndef LMP_FIX_LBM_PSM_FLUID_H
#define LMP_FIX_LBM_PSM_FLUID_H

#include <algorithm>
#include <cmath>
#include <cstring>
#include <utility>

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_property_atom.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "update.h"

#include "LBM_PSM_BGK_dynamics.h"
#include "LBM_PSM_MPICOMM.h"
#include "LBM_PSM_exchangeParticleData.h"
#include "LBM_PSM_unitConversion.h"

using namespace LAMMPS_NS;
using namespace FixConst;

class fix_LBM_PSM : public Fix {

  public:
    fix_LBM_PSM(class LAMMPS *, int, char **);
    ~fix_LBM_PSM();
    int setmask();
    void init();
    void post_force(int);

    LBMPSMBGKDynamics *dynamics;

    UnitConversion *unitConversion;
    ExchangeParticleData *exchangeParticleData;
    LBMPSMMPI *lbmmpicomm;

  private:
    int Nlc;            // Number of lattice grid nodes discretising the characteristic length lc
    double lc;          // characteristic length
    double rho;         // fluid density
    double nu;          // fluid viscosity
    double Re;          // Reynolds number of system (based on characteristic velocity, characteristic length, and fluid viscosity)
    double tau;         // BGK relaxation parameter (optional, default is tau = 0.7)

  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  int pack_reverse_comm_size(int, int);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

  double **hydrodynamicInteractions;

};

#endif
#endif