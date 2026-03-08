/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch
------------------------------------------------------*/

#include "fix_LBM_PSM_writeVTI.h"

WriteVTI::WriteVTI(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {

  if (narg < 5) error->all(FLERR,"Illegal fix lbm-psm-vti command");

  for(int ifix=0; ifix<modify->nfix; ifix++)
    if(strcmp(modify->fix[ifix]->style,"lbm-psm")==0)
      fixLBMPSM = (fix_LBM_PSM *)modify->fix[ifix];

  fileName = "flow_field";
  binary_ = true;
  useDouble_ = false;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm-vti command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix lbm-psm-vti command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"binary") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm-vti command");
      int readBinaryOption = atoi(arg[iarg+1]);
      if (readBinaryOption < 0 || readBinaryOption > 1) error->all(FLERR,"Illegal fix lbm-psm-vti command");
      binary_ = (bool)readBinaryOption;
      iarg += 2;
    } else if (strcmp(arg[iarg],"useDouble") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm-vti command");
      int readDoubleOption = atoi(arg[iarg+1]);
      if (readDoubleOption < 0 || readDoubleOption > 1) error->all(FLERR,"Illegal fix lbm-psm-vti command");
      useDouble_ = (bool)readDoubleOption;
      iarg += 2;
    } else error->all(FLERR,"Illegal fix lbm-psm-vti command");
  }
};


WriteVTI::~WriteVTI() {};


int WriteVTI::setmask()
{
  int mask =0;
  mask |= PRE_FORCE;
  return mask;
}


void WriteVTI::init()
{
  decomposition[0] = comm->procgrid[0];
  decomposition[1] = comm->procgrid[1];
  decomposition[2] = comm->procgrid[2];

  nx = fixLBMPSM->dynamics->get_nx();
  ny = fixLBMPSM->dynamics->get_ny();
  nz = fixLBMPSM->dynamics->get_nz();

  write_vti_wrapper(fileName, 0, binary_, useDouble_,
            fixLBMPSM->dynamics->get_x_reference(),
            fixLBMPSM->dynamics->get_y_reference(),
            fixLBMPSM->dynamics->get_z_reference(),
            fixLBMPSM->dynamics->get_B_reference(), 1.0,
            fixLBMPSM->dynamics->get_rho_reference(), fixLBMPSM->get_rho(),
            fixLBMPSM->dynamics->get_u_reference(), 0.0);
}


void WriteVTI::pre_force(int)
{
  if (update->ntimestep % nevery) return;
  double u_infty = fixLBMPSM->unitConversion->get_u_lb();
  double Uc = fixLBMPSM->unitConversion->get_Uc();

  write_vti_wrapper(fileName, update->ntimestep, binary_, useDouble_,
            fixLBMPSM->dynamics->get_x_reference(),
            fixLBMPSM->dynamics->get_y_reference(),
            fixLBMPSM->dynamics->get_z_reference(),
            fixLBMPSM->dynamics->get_B_reference(), 1.0,
            fixLBMPSM->dynamics->get_rho_reference(), fixLBMPSM->get_rho(),
            fixLBMPSM->dynamics->get_u_reference(), Uc/u_infty);

    
}



void WriteVTI::write_vti_wrapper(std::string name_, int timestep, bool binary, bool useDouble,
                         vector<double> &x_,
                         vector<double> &y_,
                         vector<double> &z_,
                         std::vector<double> &rho_, double rho0_,
                         std::vector<double> &u_, double u0_,
                         std::vector<double> &B_, double B0_)
{
    if (useDouble) {
        execute_write_vti<double>(name_, timestep, binary, x_, y_, z_, rho_, rho0_, u_, u0_, B_, B0_);
    } else {
        execute_write_vti<float>(name_, timestep, binary, x_, y_, z_, rho_, rho0_, u_, u0_, B_, B0_);
    }
}
