/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#include "fix_LBM_PSM_writeVTK.h"

WriteVTK::WriteVTK(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {

  if (narg < 5) error->all(FLERR,"Illegal fix lbm-psm-vtk command");

  for(int ifix=0; ifix<modify->nfix; ifix++)
    if(strcmp(modify->fix[ifix]->style,"lbm-psm")==0)
      fixLBMPSM = (fix_LBM_PSM *)modify->fix[ifix];

  fileName = "flow_field";
  binary_ = true;
  useDouble_ = false;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm-vtk command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix lbm-psm-vtk command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"binary") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm-vtk command");
      int readBinaryOption = atoi(arg[iarg+1]);
      if (readBinaryOption < 0 || readBinaryOption > 1) error->all(FLERR,"Illegal fix lbm-psm-vtk command");
      binary_ = (bool)readBinaryOption;
      iarg += 2;
    } else if (strcmp(arg[iarg],"useDouble") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm-vtk command");
      int readDoubleOption = atoi(arg[iarg+1]);
      if (readDoubleOption < 0 || readDoubleOption > 1) error->all(FLERR,"Illegal fix lbm-psm-vtk command");
      useDouble_ = (bool)readDoubleOption;
      iarg += 2;
    } else error->all(FLERR,"Illegal fix lbm-psm-vtk command");
  }
};


WriteVTK::~WriteVTK() {};


int WriteVTK::setmask()
{
  int mask =0;
  mask |= PRE_FORCE;
  return mask;
}


void WriteVTK::init()
{
  decomposition[0] = comm->procgrid[0];
  decomposition[1] = comm->procgrid[1];
  decomposition[2] = comm->procgrid[2];

  nx = fixLBMPSM->dynamics->get_nx();
  ny = fixLBMPSM->dynamics->get_ny();
  nz = fixLBMPSM->dynamics->get_nz();

  write_vtk_wrapper(fileName, 0, binary_, useDouble_,
            fixLBMPSM->dynamics->get_x_reference(), 1.0,
            fixLBMPSM->dynamics->get_y_reference(), 1.0,
            fixLBMPSM->dynamics->get_z_reference(), 1.0,
            fixLBMPSM->dynamics->get_B_reference(), 1.0,
            fixLBMPSM->dynamics->get_rho_reference(), fixLBMPSM->get_rho(),
            fixLBMPSM->dynamics->get_u_reference(), 0.0);
}


void WriteVTK::pre_force(int)
{
  if (update->ntimestep % nevery) return;
  double u_infty = fixLBMPSM->unitConversion->get_u_lb();
  double Uc = fixLBMPSM->unitConversion->get_Uc();

  write_vtk_wrapper(fileName, update->ntimestep, binary_, useDouble_,
            fixLBMPSM->dynamics->get_x_reference(), 1.0,
            fixLBMPSM->dynamics->get_y_reference(), 1.0,
            fixLBMPSM->dynamics->get_z_reference(), 1.0,
            fixLBMPSM->dynamics->get_B_reference(), 1.0,
            fixLBMPSM->dynamics->get_rho_reference(), fixLBMPSM->get_rho(),
            fixLBMPSM->dynamics->get_u_reference(), Uc/u_infty);

    
}



void WriteVTK::write_vtk_wrapper(std::string name_, int timestep, bool binary, bool useDouble,
                         vector<double> &x_, double x0_,
                         vector<double> &y_, double y0_,
                         vector<double> &z_, double z0_,
                         std::vector<double> &rho_, double rho0_,
                         std::vector<double> &u_, double u0_,
                         std::vector<double> &B_, double B0_) 
{
    if (useDouble) {
        execute_write_vtk<double>(name_, timestep, binary, x_, x0_,
                         y_, y0_,
                         z_, z0_,rho_, rho0_, u_, u0_, B_, B0_);
    } else {
        execute_write_vtk<float>(name_, timestep, binary, x_, x0_,
                         y_, y0_,
                         z_, z0_,rho_, rho0_, u_, u0_, B_, B0_);
    }
}
