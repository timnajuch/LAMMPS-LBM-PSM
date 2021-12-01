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

  virial_global_flag = virial_peratom_flag = 1;
  thermo_virial = 1;

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
  lbmmpicomm = new PSM_LBM_MPI(world, decomposition, procNeigh, procCoordinates, domain->dimension);

  unitConversion = new Unit_Conversion(rho, nu, lc, Re, Nlc, tau, domain->dimension);

  update->dt = unitConversion->get_phys_time(1.0)/((double)nevery);

  vector<double> F_lbm(3,0.0); // external forcing such as gravity. not incorporated yet.
/*
  if ( fmod(domain->xprd, unitConversion->get_dx()) != 0
    || fmod(domain->xprd, unitConversion->get_dx()) != 0
    || fmod(domain->xprd, unitConversion->get_dx()) != 0 ){
      error->all(FLERR, "Illegal cell width. Division of domain length and cell width has to give an integer.");
  }
*/
double SMALL = 1e-15;
  if(  pow(((double)((int)(domain->xprd/unitConversion->get_dx()+0.5)) - domain->xprd/unitConversion->get_dx()), 2.0) > SMALL
    || pow(((double)((int)(domain->yprd/unitConversion->get_dx()+0.5)) - domain->yprd/unitConversion->get_dx()), 2.0) > SMALL
    || pow(((double)((int)(domain->zprd/unitConversion->get_dx()+0.5)) - domain->zprd/unitConversion->get_dx()), 2.0) > SMALL ){
      error->all(FLERR, "Illegal cell width. Division of domain length and cell width has to give an integer.");
  }

  if ( ((int)(domain->xprd/unitConversion->get_dx()+1.5) % decomposition[0] != 0)
    || ((int)(domain->yprd/unitConversion->get_dx()+1.5) % decomposition[1] != 0)
    || ((int)(domain->zprd/unitConversion->get_dx()+1.5) % decomposition[2] != 0) ){
      error->all(FLERR, "Illegal combination of decomposition and lattice cell number. Division of total lattice cell number by domain decomposition has to be even (i.e. modulo == 0) in each direction.");
  }
  

  int nx = domain->xprd/unitConversion->get_dx()+1;
  int ny = domain->yprd/unitConversion->get_dx()+1;
  int nz = 0;
  if (domain->dimension == 3)
    {nz = domain->zprd/unitConversion->get_dx()+1; }

  vector<double> boxLength{domain->xprd, domain->yprd, domain->zprd};
  vector<double> origin{domain->boxlo[0], domain->boxlo[1], domain->boxlo[2]};

  dynamics = new BGK_GuoExtForce_Dynamics2D(tau, nx, ny, nz, 9, F_lbm, decomposition, procCoordinates, origin, boxLength, domain->dimension, unitConversion->get_dx());

  dynamics->initialise_domain(unitConversion->get_dx(), unitConversion->get_dx(), unitConversion->get_dx());

  dynamics->initialise_dynamics(1.0, 0.0, 0.0, 0.0);

  exchangeParticleData = new ExchangeParticleData(domain->dimension);

  if(!fix_hydroForce_)
  {
    char *fixarg[] = {
      (char *)"hydroForce",     // fix id
      (char *)"all",            // fix group
      (char *)"property/atom",  // fix style: property/atom
      (char *)"d2_hydroForce",  // name
      (char *)"3",              // number of array columns
      (char *)"ghost",          // communicate ghost atom
      (char *)"yes"
    };
 
    modify->add_fix(7, fixarg, 1);
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
      (char *)"6",              // number of array columns
      (char *)"ghost",          // communicate ghost atom
      (char *)"yes"
    };
    modify->add_fix(7, fixarg, 1);
    fix_stresslet_ = (FixPropertyAtom *) modify->fix[modify->nfix-1];
  }
}


void fix_PSM_LBM::pre_force(int vflag)
{
  if (update->ntimestep % nevery) return;
  int nPart = atom->nlocal + atom->nghost;
  vector<double> boxLength{domain->xprd, domain->yprd, domain->zprd};
  vector<double> origin{domain->boxlo[0], domain->boxlo[1], domain->boxlo[2]};
  exchangeParticleData->setParticlesOnLattice(dynamics, unitConversion, nPart, atom->tag, atom->x, atom->v, atom->omega, atom->radius, boxLength, origin);

  if(domain->dimension == 2){
    lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), false, 0, dynamics->get_nx(), dynamics->get_ny(), 1, dynamics->get_envelopeWidth(), domain->xperiodic);
    lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), false, 1, dynamics->get_nx(), dynamics->get_ny(), 1, dynamics->get_envelopeWidth(), domain->yperiodic);
  }else{
    lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), false, 0, dynamics->get_nx(), dynamics->get_ny(), dynamics->get_nz(), dynamics->get_envelopeWidth(), domain->xperiodic);
    lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), false, 1, dynamics->get_nx(), dynamics->get_ny(), dynamics->get_nz(), dynamics->get_envelopeWidth(), domain->yperiodic);
    lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), false, 2, dynamics->get_nx(), dynamics->get_ny(), dynamics->get_nz(), dynamics->get_envelopeWidth(), domain->zperiodic);
  }

  dynamics->macroCollideStream();

  double **f = atom->f;
  double **t = atom->torque;

  int flagHF, ncolumnsHF;
  int hydroForceFixID = atom->find_custom((char *)"hydroForce", flagHF, ncolumnsHF);
  int flagHT, ncolumnsHT;
  int hydroTorqueFixID = atom->find_custom((char *)"hydroTorque", flagHT, ncolumnsHT);
  int flagS, ncolumnsS;
  int stressletFixID = atom->find_custom((char *)"stresslet", flagS, ncolumnsS);

  v_init(vflag);


  for(int i=0;i<nPart;i++){

    vector<double> fh;
    fh.resize(3);
    fh[0] = 0.0;
    fh[1] = 0.0;
    fh[2] = 0.0;

    vector<double> th;
    th.resize(3);
    th[0] = 0.0;
    th[1] = 0.0;
    th[2] = 0.0;

    vector<double> stresslet;
    stresslet.resize(6);
    stresslet[0] = 0.0;
    stresslet[1] = 0.0;
    stresslet[2] = 0.0;
    stresslet[3] = 0.0;
    stresslet[4] = 0.0;
    stresslet[5] = 0.0;

    exchangeParticleData->calculateHydrodynamicInteractions(dynamics, unitConversion, atom->x[i], atom->radius[i], fh, th, stresslet);

    if (i < atom->nlocal){
      f[i][0] += fh[0];
      f[i][1] += fh[1];
      f[i][2] += fh[2];
    }else{ // ghost atoms
      f[i][0] = fh[0];
      f[i][1] = fh[1];
      f[i][2] = fh[2];
    }


    if (i < atom->nlocal){
      t[i][0] += th[0];
      t[i][1] += th[1];
      t[i][2] += th[2];
    }else{ // ghost atoms
      t[i][0] = th[0];
      t[i][1] = th[1];
      t[i][2] = th[2];
    }

    // Torque added because it is the antisymmetric part 
    // of the first moment of the stress over the surface
    // Stress defined as negative, hence the flipped signs
    // TODO: Check what middle indice of darray stands for. [][1][] gives seg fault
//    if (i < atom->nlocal){
      double stresslet_arr[6] = {stresslet[0], stresslet[1], stresslet[2], stresslet[3], stresslet[4], stresslet[5]};
      v_tally(i, stresslet_arr);
//    }

/*
    if (i < atom->nlocal){
      virial[0] += stresslet[0];
      virial[1] += stresslet[1];
      virial[2] += stresslet[2];
      virial[3] += stresslet[3];// - 0.5*th[2];
      //virial[4] += stresslet[4];// + 0.5*th[1];
      virial[4] += th[2];
      virial[5] += stresslet[5];// - 0.5*th[0];
    }else{
      virial[0] = stresslet[0];
      virial[1] = stresslet[1];
      virial[2] = stresslet[2];
      virial[3] = stresslet[3];// - 0.5*th[2];
      //virial[4] = stresslet[4];// + 0.5*th[1];
      virial[4] = th[2];
      virial[5] = stresslet[5];// - 0.5*th[0];
    }
*/
  }
  comm->reverse_comm();// todo check if sufficient or if i need to define some virtual functions etc

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
