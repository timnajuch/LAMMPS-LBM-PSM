/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(lbm-psm-restart,LBMPSMWriteReadRestart)

#else

#ifndef LBMPWMWRITEREADRESTART_H
#define LBMPSMWRITEREADRESTART_H

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

class LBMPSMWriteReadRestart : public Fix{
  private:
    int iWrite, iRead; // Write or read when set to 1.

  public:
    LBMPSMWriteReadRestart(class LAMMPS *, int, char **);
    ~LBMPSMWriteReadRestart();
    int setmask();
    void init();
    void pre_force(int);

  void write_restart(string name_, vector<double> &f_, int currentStep);
  void read_restart(string name_, vector<double> &f_);

  class fix_LBM_PSM *fixLBMPSM;

};

#endif
#endif
