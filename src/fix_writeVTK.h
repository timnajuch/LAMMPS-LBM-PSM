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
#include <iostream>
#include <math.h>
#include <vector>

#include "comm.h"
#include "fix.h"
#include "modify.h"

#include "fix_PSM_LBM.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

class WriteVTK : public Fix{
  private:
    int nx, ny;
    int decomposition[3];
    //vector<double> point_id;

  public:
    //WriteVTK(class LAMMPS *, int, char **, int nx_, int ny_, int decomposition_[3]);
    WriteVTK(class LAMMPS *, int, char **);
    ~WriteVTK();
    int setmask();
    void init();
    void pre_force(int);

//  void write_vtk(string name_, vector<double> &x_, vector<double> &y_, vector<double> &B_, vector<double> &rho_, vector<double> &u_, MPICommunication &mpicomm);
//  void write_vtk(string name_, vector<double> &x_, double x0_, vector<double> &y_, double y0_, vector<double> &B_, double B0_, vector<double> &rho_, double rho0_, vector<double> &u_, double u0_, MPICommunication &mpicomm);
  void write_vtk(string name_, vector<double> &x_, vector<double> &y_, vector<double> &B_, vector<double> &rho_, vector<double> &u_);
  void write_vtk(string name_, vector<double> &x_, double x0_, vector<double> &y_, double y0_, vector<double> &B_, double B0_, vector<double> &rho_, double rho0_, vector<double> &u_, double u0_);

  void write_profile(string name_, int ix_, vector<double> &y_, double y0_, vector<double> &B_, double B0_, vector<double> &rho_, double rho0_, vector<double> &u_, double u0_);

  //vector<double> scale_vector(vector<double> &vec_, double scaling_);
  void scale_vector(vector<double> &vec_, double scaling_);

  class fix_PSM_LBM *fixPSMLBM;

};

#endif
#endif
