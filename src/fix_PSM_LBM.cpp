/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */


#include "fix_PSM_LBM.h"
/*
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

#include "BGK_GuoExtForce_dynamics2D.h"
//#include "mpiCommunication.h"
#include "fix_PSM_LBM_MPICOMM.h"

using namespace LAMMPS_NS;
using namespace FixConst;
*/

fix_PSM_LBM::fix_PSM_LBM(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
//  int decomposition[3] = {comm->procgrid[0], comm->procgrid[1], comm->procgrid[2]};
//  int procCoordinates[3] = {0, 0, 0};
//  vector<double> F_lbm(2,0.0);
}

fix_PSM_LBM::~fix_PSM_LBM()
{
  delete dynamics;
}

int fix_PSM_LBM::setmask()
{
  int mask =0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= END_OF_STEP;
  return mask;
}


void fix_PSM_LBM::init()
{

//  MPICommunication mpicomm(&argc, &argv);
  int decomposition[3] = {comm->procgrid[0], comm->procgrid[1], comm->procgrid[2]};
//  mpicomm.domainDecomposition(decomposition);

  int procCoordinates[3] = {0, 0, 0}; // todo: get proc coordinate of current proc from lammps mpi class
//  vector<int> procCoordinatesTmp;
//  procCoordinatesTmp.resize(3);
//  mpicomm.returnProcCoordinatesArray(procCoordinatesTmp);
//  procCoordinates[0] = procCoordinatesTmp[0];
//  procCoordinates[1] = procCoordinatesTmp[1];
//  procCoordinates[2] = procCoordinatesTmp[2];
/*
  int procNeigh[3][2];
    = comm->procneigh[0][0];
    = comm->procneigh[0][1];
   = comm->procneigh[1][0];
   = comm->procneigh[1][1];
    = comm->procneigh[2][0];
  procNeigh[2][1]    = comm->procneigh[2][1];
*/
  int procNeigh[6];
  procNeigh[0] = comm->procneigh[0][0];
  procNeigh[1] = comm->procneigh[0][1];
  procNeigh[2] = comm->procneigh[1][0];
  procNeigh[3] = comm->procneigh[1][1];
  procNeigh[4] = comm->procneigh[2][0];
  procNeigh[5] = comm->procneigh[2][1];
//std::cout << "DEBUG A" << std::endl;
  lbmmpicomm = new PSM_LBM_MPI(world, decomposition, procNeigh);
//std::cout << "DEBUG B" << std::endl;

//  BGK_GuoExtForce_Dynamics2D dynamics(tau, nx, ny, q, F_lbm, decomposition, procCoordinates);
  vector<double> F_lbm(2,0.0);
//  BGK_GuoExtForce_Dynamics2D dynamics(0.7, 100, 50, 9, F_lbm, decomposition, procCoordinates);
  dynamics = new BGK_GuoExtForce_Dynamics2D(0.7, 100, 50, 9, F_lbm, decomposition, procCoordinates);


  double rhof = 1000.0;
  double Re = 0.1;
  double nu = 0.0001;
  double tau = 0.7;
  int N_particle = 10;
  double dp = 0.01;   // particle diameter
  unitConversion = new Unit_Conversion(rhof, nu, dp, Re, N_particle, tau);

  dynamics->initialise_domain(unitConversion->get_dx(), unitConversion->get_dx());

  dynamics->initialise_dynamics(1.0, 0.0, 0.0);
/*
  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix addforce command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix addforce command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
    } else error->all(FLERR,"Illegal fix addforce command");
  }
*/


  int nPart = atom->nlocal + atom->nghost;
  std::cout << "nPart: " << nPart << std::endl;

  //double dp = 10.0; /a/ 10 lattice cells over particle diameter
  vector<double> xp;
  vector<double> us;

  xp.resize(2);
  //xp[0] = 0.5*(double)dynamics->get_nx();
  //xp[1] = 0.5*(double)dynamics->get_ny();
  xp[0] = atom->x[0][0]/domain->xprd*(double)dynamics->get_nx();
  xp[1] = atom->x[0][1]/domain->yprd*(double)dynamics->get_ny();

  us.resize(2);
  us[0] = 0.0;
  us[1] = 0.0;
  
  exchangeParticleData = new ExchangeParticleData((double)N_particle, xp, us);
  //exchangeParticleData->setParticlesOnLattice(dynamics);
  exchangeParticleData->setParticlesOnLattice(dynamics, nPart, atom->x);

}

void fix_PSM_LBM::pre_force(int)
{
//  if (update->ntimestep % nevery) return;

  int envelopeWidth = 1;
  lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), false, 0, dynamics->get_nx(), dynamics->get_ny(), 1, envelopeWidth, true);

  dynamics->macroCollideStream();
}

int fix_PSM_LBM::get_nx()
{
  return dynamics->get_nx();
}

int fix_PSM_LBM::get_ny()
{
  return dynamics->get_ny();
}

vector<double> fix_PSM_LBM::get_x()
{
  return dynamics->get_x();
}

vector<double> fix_PSM_LBM::get_y()
{
  return dynamics->get_y();
}

vector<double> fix_PSM_LBM::get_rho()
{
  return dynamics->get_rho();
}

vector<double> fix_PSM_LBM::get_u()
{
  return dynamics->get_u();
}

vector<double> fix_PSM_LBM::get_B()
{
  return dynamics->get_B();
}
