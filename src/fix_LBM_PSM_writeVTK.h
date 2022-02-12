/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(lbm-psm-vtk,WriteVTK)

#else

#ifndef WRITEVTK_H
#define WRITEVTK_H

#include <algorithm>
#include <fstream> 
#include <functional>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vector>
#include <sstream>
#include <string>

#include "comm.h"
#include "fix.h"
#include "modify.h"

#include "fix_LBM_PSM.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

class WriteVTK : public Fix{
  private:
    int nx, ny, nz;
    int decomposition[3];

  public:
    WriteVTK(class LAMMPS *, int, char **);
    ~WriteVTK();
    int setmask();
    void init();
    void pre_force(int);

  void write_vtk(string name_, vector<double> &x_, double x0_, vector<double> &y_, double y0_, vector<double> &z_, double z0_, vector<double> &B_, double B0_, vector<double> &rho_, double rho0_, vector<double> &u_, double u0_);

  void write_profile(string name_, int ix_, vector<double> &y_, double y0_, vector<double> &B_, double B0_, vector<double> &rho_, double rho0_, vector<double> &u_, double u0_);

  void scale_vector(vector<double> &vec_, double scaling_);

  class fix_LBM_PSM *fixLBMPSM;

};

#endif
#endif
