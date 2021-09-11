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

FixStyle(lbm-psm-bc,fix_PSM_LBM_BC)

#else

#ifndef LMP_FIX_PSM_LBM_BC_H
#define LMP_FIX_PSM_LBM_BC_H


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
#include "fix.h"
#include "comm.h"

#include "fix_PSM_LBM.h"
#include "zou_he_BC_2D.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

class fix_PSM_LBM_BC : public Fix {

public:
  fix_PSM_LBM_BC(class LAMMPS *, int, char **);
  ~fix_PSM_LBM_BC();
  int setmask();
  void init();
//    void initial_integrate(int);
//    void setup(int);
  void pre_force(int);
//    void end_of_step();

  ZouHeBC2D *zouHe2D;

  class fix_PSM_LBM *fixPSMLBM;

};

//}
#endif
#endif

