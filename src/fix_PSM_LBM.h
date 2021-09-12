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

#include "BGK_GuoExtForce_dynamics2D.h"
#include "PSM_LBM_MPICOMM.h"
#include "fix_PSM_LBM_BC.h"
#include "ExchangeParticleData.h"
#include "unit_conversion.h"

using namespace LAMMPS_NS;
using namespace FixConst;

class fix_PSM_LBM : public Fix {

  public:
    fix_PSM_LBM(class LAMMPS *, int, char **);
    ~fix_PSM_LBM();
    int setmask();
    void init();
    void pre_force(int);

    BGK_GuoExtForce_Dynamics2D *dynamics;
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

    double **get_force_ptr();
    double **get_torque_ptr();
    double **get_stresslet_ptr();

  private:
    int Nlc;            // Number of lattice grid nodes discretising the characteristic length lc
    double lc;          // characteristic length
    double rho;         // fluid density
    double nu;          // fluid viscosity
    double Re;          // Reynolds number of system (based on characteristic velocity, characteristic length, and fluid viscosity)
    double tau;         // BGK relaxation parameter (optional, default is tau = 0.7)

    class FixPropertyAtom *fix_hydroForce_;
    class FixPropertyAtom *fix_hydroTorque_;
    class FixPropertyAtom *fix_stresslet_; 
};

#endif
#endif

// TODO add comminucation of hyd force etc to ghost atoms every step (otherwise only communicated when neighborl list newly build)
