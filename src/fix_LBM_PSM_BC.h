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

FixStyle(lbm-psm-bc,fix_LBM_PSM_BC)

#else

#ifndef LMP_FIX_LBM_PSM_BC_H
#define LMP_FIX_LBM_PSM_BC_H

#include <algorithm>
#include <cmath>
#include <cstring>
#include <utility>

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "update.h"

#include "fix_LBM_PSM.h"
#include "LBM_PSM_zou_he_BC.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

class fix_LBM_PSM_BC : public Fix {

  public:
    fix_LBM_PSM_BC(class LAMMPS *, int, char **);
    ~fix_LBM_PSM_BC();
    int setmask();
    void init();
    void post_force(int);

    ZouHeBC *zouHe;

    class fix_LBM_PSM *fixLBMPSM;

    int typeBC; // 0: periodic set-up, 1: shear set-up, 2: imposed flow in x-direction set-up

/* TODO Use when genreral BC are implemented
    int lowerBC[3];  // 0: periodic, 1: velocity, 2: density
    int upperBC[3];  // 0: periodic, 1: velocity, 2: density
//    int lowerBC_y;  // 0: periodic, 1: velocity, 2: density
//    int upperBC_y;  // 0: periodic, 1: velocity, 2: density
//    int lowerBC_z;  // 0: periodic, 1: velocity, 2: density
//    int upperBC_z;  // 0: periodic, 1: velocity, 2: density

    double lowerRhoBC[3];
    double lowerVelBC[3];
    double upperRhoBC[3];
    double upperVelBC[3];
*/
};

#endif
#endif
