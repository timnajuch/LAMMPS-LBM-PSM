/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#include "fix_LBM_PSM_writeReadRestart.h"

LBMPSMWriteReadRestart::LBMPSMWriteReadRestart(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {

  if (narg < 9) error->all(FLERR,"Illegal fix lbm-psm-restart command");

  for(int ifix=0; ifix<modify->nfix; ifix++)
    if(strcmp(modify->fix[ifix]->style,"lbm-psm")==0)
      fixLBMPSM = (fix_LBM_PSM *)modify->fix[ifix];

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm-restart command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix lbm-psm-restart command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"write") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm-restart command");
      iWrite = atoi(arg[iarg+1]);
      if (iWrite < 0 || iWrite > 1) error->all(FLERR,"Illegal fix lbm-psm-restart command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"read") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm-restart command");
      iRead = atoi(arg[iarg+1]);
      if (iRead < 0 || iRead > 1) error->all(FLERR,"Illegal fix lbm-psm-restart command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix lbm-psm-restart command");
  }
};


LBMPSMWriteReadRestart::~LBMPSMWriteReadRestart() {};


int LBMPSMWriteReadRestart::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}


void LBMPSMWriteReadRestart::init()
{
  if (iRead == 1){
    ostringstream restartFileStringTmp;
    restartFileStringTmp << "LBM-PSM-restartFile-processor-" << comm->me << ".bin";
    string timeString(restartFileStringTmp.str());
    read_restart(timeString, fixLBMPSM->dynamics->getVector_f());
  }
}


void LBMPSMWriteReadRestart::pre_force(int)
{
  if (update->ntimestep % nevery) return;
  if (iWrite == 1){
    ostringstream restartFileStringTmp;
    restartFileStringTmp << "LBM-PSM-restartFile-processor-" << comm->me << ".bin";
    string timeString(restartFileStringTmp.str());

    write_restart(timeString, fixLBMPSM->dynamics->getVector_f(), fixLBMPSM->dynamics->get_currentStep());
  }
}


void LBMPSMWriteReadRestart::write_restart(string name_, vector<double> &f_, int currentStep){
  ofstream out;
  out.open(name_, ios::binary | ios::out);

  out.write ( reinterpret_cast<char *>(&currentStep), sizeof(currentStep) );

  unsigned int f_size = f_.size();
  out.write ( reinterpret_cast<char *>(&f_size), sizeof(unsigned) );
  out.write ( reinterpret_cast<char *>(&f_[0]), f_.size()*sizeof(double) );

  out.close();
}


void LBMPSMWriteReadRestart::read_restart(string name_, vector<double> &f_){
  ifstream in(name_, ios_base::binary | ios::in);

  int currentStep = 0;
  in.read( reinterpret_cast<char *>(&currentStep), sizeof(int) );

  unsigned vsize;
  in.read( reinterpret_cast<char *>(&vsize), sizeof(unsigned) );
  vector<double> f_in(vsize);
  in.read( reinterpret_cast<char *>(&f_in[0]), vsize*sizeof(double) );

  in.close();

  fixLBMPSM->dynamics->set_currentStep(currentStep);
  f_ = f_in;
}