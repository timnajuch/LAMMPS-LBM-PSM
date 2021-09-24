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

fix_PSM_LBM::fix_PSM_LBM(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 15) error->all(FLERR,"Illegal fix lbm-psm command");

  fix_hydroForce_ = nullptr;
  fix_hydroTorque_ = nullptr;
  fix_stresslet_ = nullptr;

  tau = 0.7; // default value for BGK relaxation parameter

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix lbm-psm command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"Nlc") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm command");
      Nlc = atoi(arg[iarg+1]);
      if (Nlc <= 0) error->all(FLERR,"Illegal fix lbm-psm command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"lc") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm command");
      lc = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (lc <= 0.0) error->all(FLERR,"Illegal fix lbm-psm command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"rho") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm command");
      rho = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (rho <= 0.0) error->all(FLERR,"Illegal fix lbm-psm command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"nu") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm command");
      nu = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (nu <= 0.0) error->all(FLERR,"Illegal fix lbm-psm command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"Re") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm command");
      Re = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (Re <= 0.0) error->all(FLERR,"Illegal fix lbm-psm command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"tau") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm command");
      tau = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (tau <= 0.5) error->all(FLERR,"Illegal fix lbm-psm command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix lbm-psm command");
  }
}


fix_PSM_LBM::~fix_PSM_LBM()
{
  delete dynamics;
}


int fix_PSM_LBM::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= END_OF_STEP;
  return mask;
}


void fix_PSM_LBM::init()
{
  int decomposition[3] = {comm->procgrid[0], comm->procgrid[1], comm->procgrid[2]};
  int procCoordinates[3] = {comm->myloc[0], comm->myloc[1], comm->myloc[2]};
  int procNeigh[6];
  procNeigh[0] = comm->procneigh[0][0];
  procNeigh[1] = comm->procneigh[0][1];
  procNeigh[2] = comm->procneigh[1][0];
  procNeigh[3] = comm->procneigh[1][1];
  procNeigh[4] = comm->procneigh[2][0];
  procNeigh[5] = comm->procneigh[2][1];
  lbmmpicomm = new PSM_LBM_MPI(world, decomposition, procNeigh, procCoordinates);

  unitConversion = new Unit_Conversion(rho, nu, lc, Re, Nlc, tau);

  vector<double> F_lbm(2,0.0); // external forcing such as gravity. not incorporated yet.
  int nx = domain->xprd/unitConversion->get_dx()+1;
  int ny = domain->yprd/unitConversion->get_dx()+1;
  vector<double> boxLength{domain->xprd, domain->yprd, domain->zprd};
  vector<double> origin{domain->boxlo[0], domain->boxlo[1], domain->boxlo[2]};

  dynamics = new BGK_GuoExtForce_Dynamics2D(tau, nx, ny, 9, F_lbm, decomposition, procCoordinates, origin, boxLength);

  dynamics->initialise_domain(unitConversion->get_dx(), unitConversion->get_dx());

  dynamics->initialise_dynamics(1.0, 0.0, 0.0);
  exchangeParticleData = new ExchangeParticleData();

  if(!fix_hydroForce_)
  {
    char *fixarg[] = {
      (char *)"hydroForce",     // fix id
      (char *)"all",            // fix group
      (char *)"property/atom",  // fix style: property/atom
      (char *)"d_hydroForce",   // name
      (char *)"ghost",          // communicate ghost atom
      (char *)"yes"
    };
 
    modify->add_fix(6, fixarg, 1);
    fix_hydroForce_ = (FixPropertyAtom *) modify->fix[modify->nfix-1];
  }

  if(!fix_hydroTorque_)
  {
    char *fixarg[] = {
      (char *)"hydroTorque",    // fix id
      (char *)"all",            // fix group
      (char *)"property/atom",  // fix style: property/atom
      (char *)"d_hydroTorque",  // name
      (char *)"ghost",          // communicate ghost atom
      (char *)"yes"
    };
    modify->add_fix(6, fixarg, 1);
    fix_hydroTorque_ = (FixPropertyAtom *) modify->fix[modify->nfix-1];
  }

  if(!fix_stresslet_)
  {
    char *fixarg[] = {
      (char *)"stresslet",      // fix id
      (char *)"all",            // fix group
      (char *)"property/atom",  // fix style: property/atom
      (char *)"d2_stresslet",   // name
      (char *)"3",              // number of array columns
      (char *)"ghost",          // communicate ghost atom
      (char *)"yes"
    };
    modify->add_fix(7, fixarg, 1);
    fix_stresslet_ = (FixPropertyAtom *) modify->fix[modify->nfix-1];
  }
}


void fix_PSM_LBM::pre_force(int)
{
  if (update->ntimestep % nevery) return;

  int nPart = atom->nlocal + atom->nghost;
  vector<double> boxLength{domain->xprd, domain->yprd, domain->zprd};
  vector<double> origin{domain->boxlo[0], domain->boxlo[1], domain->boxlo[2]};
  exchangeParticleData->setParticlesOnLattice(dynamics, unitConversion, nPart, atom->tag, atom->x, atom->v, atom->radius, boxLength, origin);

  //lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), false, 0, dynamics->get_nx(), dynamics->get_ny(), 1, dynamics->get_envelopeWidth(), false);
  lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), false, 0, dynamics->get_nx(), dynamics->get_ny(), 1, dynamics->get_envelopeWidth(), domain->xperiodic);
  lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), false, 1, dynamics->get_nx(), dynamics->get_ny(), 1, dynamics->get_envelopeWidth(), domain->yperiodic);
  //lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), true, 2, dynamics->get_nx(), dynamics->get_ny(), dynamics->get_nz(), dynamics->get_envelopeWidth(), domain->zperiodic);

  dynamics->macroCollideStream();


  double **f = atom->f;
  double **t = atom->torque;

  int flagHF, ncolumnsHF;
  int hydroForceFixID = atom->find_custom((char *)"hydroForce", flagHF, ncolumnsHF);
  int flagHT, ncolumnsHT;
  int hydroTorqueFixID = atom->find_custom((char *)"hydroTorque", flagHT, ncolumnsHT);
  int flagS, ncolumnsS;
  int stressletFixID = atom->find_custom((char *)"stresslet", flagS, ncolumnsS);


  //exchangeParticleData->calculateHydrodynamicInteractions(dynamics, unitConversion, nPart, atom->tag, atom->x, hydroForceFixID, hydroTorqueFixID, stressletFixID, atom, origin);


  for(int i=0;i<atom->nlocal;i++){
    f[i][0] += atom->dvector[hydroForceFixID][0];
    f[i][1] += atom->dvector[hydroForceFixID][1];
    f[i][2] += atom->dvector[hydroForceFixID][2];

    t[i][0] += atom->dvector[hydroTorqueFixID][0];
    t[i][1] += atom->dvector[hydroTorqueFixID][1];
    t[i][2] += atom->dvector[hydroTorqueFixID][2];

    // Torque added because it is the antisymmetric part 
    // of the first moment of the stress over the surface
    // Stress defined as negative, hence the flipped signs
    // TODO: Check what middle indice of darray stands for. [][1][] gives seg fault
    virial[0] += atom->darray[stressletFixID][0][0];
    virial[1] += atom->darray[stressletFixID][0][1];
    virial[2] += atom->darray[stressletFixID][0][2];     
    virial[3] += atom->darray[stressletFixID][0][3] - 0.5*atom->dvector[hydroTorqueFixID][2];
    virial[4] += atom->darray[stressletFixID][0][4] + 0.5*atom->dvector[hydroTorqueFixID][1];
    virial[5] += atom->darray[stressletFixID][0][5] - 0.5*atom->dvector[hydroTorqueFixID][0];
  }

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


double **fix_PSM_LBM::get_force_ptr()
{
  return fix_hydroForce_->array_atom;
}

double **fix_PSM_LBM::get_torque_ptr()
{
  return fix_hydroTorque_->array_atom;
}

double **fix_PSM_LBM::get_stresslet_ptr()
{
  return fix_stresslet_->array_atom;
}


/*
  void FixLbCouplingOnetoone::post_run()
  {
    // need one very last forward_comm to make sure 
    // that velocities on owned and ghost particles match
    timer->stamp();
    comm->forward_comm();
    timer->stamp(TIME_COMM);
  }
  void FixLbCouplingOnetoone::comm_force_torque()
  {
    fix_dragforce_->do_reverse_comm();
    fix_hdtorque_->do_reverse_comm();
    fix_stresslet_->do_reverse_comm();
  }
*/
