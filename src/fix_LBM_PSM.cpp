/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#include "fix_LBM_PSM.h"

fix_LBM_PSM::fix_LBM_PSM(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 15) error->all(FLERR,"Illegal fix lbm-psm command");

  virial_global_flag = virial_peratom_flag = 1;
  thermo_virial = 1;
  comm_reverse = 6;

  tau = 0.7; // default value for BGK relaxation parameter
  F_ext = vector<double>(3, 0.0);
  dynamicsScheme = 1;
  tau_m = 0.7;

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
    } else if (strcmp(arg[iarg],"dynamicsScheme") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm command");
      dynamicsScheme = atoi(arg[iarg+1]);
      if (dynamicsScheme < 1 || dynamicsScheme > 2) error->all(FLERR,"Illegal fix lbm-psm command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"tau_m") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm-trt command");
      tau_m = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (tau_m <= 0.5) error->all(FLERR,"Illegal fix lbm-psm-trt command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"Fext") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix lbm-psm command");
      F_ext[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      F_ext[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      F_ext[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;
    } else error->all(FLERR,"Illegal fix lbm-psm command");
  }

  // Allocate memory for storage of forces so that hydrodynamic interaction forces can be still imposed for timesteps between LBM calculations
  hydrodynamicInteractions = nullptr;
  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  for (int i=0; i<atom->nmax; i++){
    for (int j=0; j<6; j++){
      hydrodynamicInteractions[i][j] = 0.0;
    }
  }

}


fix_LBM_PSM::~fix_LBM_PSM()
{
  delete dynamics;
  delete unitConversion;
  delete exchangeParticleData;
  delete lbmmpicomm;
  memory->destroy(hydrodynamicInteractions);
}


int fix_LBM_PSM::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}


void fix_LBM_PSM::init()
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

  lbmmpicomm = new LBMPSMMPI(world, decomposition, procNeigh, procCoordinates, domain->dimension);

  unitConversion = new UnitConversion(rho, nu, lc, Re, Nlc, tau, domain->dimension);

  update->dt = unitConversion->get_phys_time(1.0)/((double)nevery);

  double Ma = unitConversion->get_u_lb()/(1.0/sqrt(3.0));
  if (Ma > 0.3){
      error->all(FLERR,"Mach number is Ma > 0.3. Aborting simulation. Choose parameters which give Ma < 0.3 for stability reasons.");
  }
  else{
    if (comm->me == 0){
      std::string mesg = fmt::format("Mach number is Ma = {:.2}\n", Ma);
      utils::logmesg(lmp,mesg);
      mesg = fmt::format("DEM timestep is dt = {:.6}\n", update->dt);
      utils::logmesg(lmp,mesg);
    }
  }

  double SMALL = 1e-15;
  if(  pow(((double)((int)(domain->xprd/unitConversion->get_dx()+0.5)) - domain->xprd/unitConversion->get_dx()), 2.0) > SMALL
    || pow(((double)((int)(domain->yprd/unitConversion->get_dx()+0.5)) - domain->yprd/unitConversion->get_dx()), 2.0) > SMALL
    || pow(((double)((int)(domain->zprd/unitConversion->get_dx()+0.5)) - domain->zprd/unitConversion->get_dx()), 2.0) > SMALL ){
      error->all(FLERR, "Illegal cell width. Division of domain length and cell width has to give an integer.");
  }

  int nx = domain->xprd/unitConversion->get_dx()+1.5;
  if (domain->xperiodic == true)
    { nx -= 1; }
  int ny = domain->yprd/unitConversion->get_dx()+1.5;
  if (domain->yperiodic == true)
    { ny -= 1; }
  int nz = 0;
  if (domain->dimension == 3)
  {
    nz = domain->zprd/unitConversion->get_dx()+1.5;
    if (domain->zperiodic == true)
      { nz -= 1; }
  }

  vector<double> boxLength{domain->xprd, domain->yprd, domain->zprd};
  vector<double> origin{domain->boxlo[0], domain->boxlo[1], domain->boxlo[2]};
  if (domain->dimension == 2)
    { origin[2] = 0.0; }

  vector<double> F_lbm;
  F_lbm = unitConversion->get_volume_force_lb(F_ext);

  if (dynamicsScheme == 1){
    dynamics = new LBMPSMBGKDynamics(tau, nx, ny, nz, F_lbm, decomposition, procCoordinates, origin, boxLength, domain->dimension);
  }else if (dynamicsScheme == 2){
    dynamics = new LBMPSMTRTDynamics(tau, tau_m, nx, ny, nz, F_lbm, decomposition, procCoordinates, origin, boxLength, domain->dimension);
  }

  dynamics->initialise_domain(unitConversion->get_dx(), unitConversion->get_dx(), unitConversion->get_dx());

  dynamics->initialise_dynamics(1.0, 0.0, 0.0, 0.0);

  exchangeParticleData = new ExchangeParticleData(domain->dimension, origin);

  for (int i=0; i<atom->nmax; i++){
    for (int j=0; j<6; j++){
      hydrodynamicInteractions[i][j] = 0.0;
    }
  }
}



void fix_LBM_PSM::post_force(int vflag)
{
  if (update->ntimestep % nevery)
  {
    double **f = atom->f;
    double **t = atom->torque;

    int nPart = atom->nlocal;
    for(int i = 0; i < nPart; i++){
        f[i][0] += hydrodynamicInteractions[i][0];
        f[i][1] += hydrodynamicInteractions[i][1];
        f[i][2] += hydrodynamicInteractions[i][2];
        t[i][0] += hydrodynamicInteractions[i][3];
        t[i][1] += hydrodynamicInteractions[i][4];
        t[i][2] += hydrodynamicInteractions[i][5];
    }
  }else{
    comm->forward_comm();

    int nPart = atom->nlocal + atom->nghost;
    vector<double> boxLength{domain->xprd, domain->yprd, domain->zprd};
    vector<double> origin{domain->boxlo[0], domain->boxlo[1], domain->boxlo[2]};

    exchangeParticleData->setParticlesOnLattice(dynamics, unitConversion, nPart, atom->tag, atom->x, atom->v, atom->omega, atom->radius);

    if(domain->dimension == 2){
      lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), 0, dynamics->get_nx(), dynamics->get_ny(), 1, dynamics->get_envelopeWidth(), domain->xperiodic, dynamics->get_currentStep());
      lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), 1, dynamics->get_nx(), dynamics->get_ny(), 1, dynamics->get_envelopeWidth(), domain->yperiodic, dynamics->get_currentStep());
    }else{
      lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), 0, dynamics->get_nx(), dynamics->get_ny(), dynamics->get_nz(), dynamics->get_envelopeWidth(), domain->xperiodic, dynamics->get_currentStep());
      lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), 1, dynamics->get_nx(), dynamics->get_ny(), dynamics->get_nz(), dynamics->get_envelopeWidth(), domain->yperiodic, dynamics->get_currentStep());
      lbmmpicomm->sendRecvData<double>(dynamics->getVector_f(), 2, dynamics->get_nx(), dynamics->get_ny(), dynamics->get_nz(), dynamics->get_envelopeWidth(), domain->zperiodic, dynamics->get_currentStep());
    }

    dynamics->macroCollideStream();

    double **f = atom->f;
    double **t = atom->torque;
    v_init(vflag);

    vector<double> fh;
    fh.resize(3);
    vector<double> th;
    th.resize(3);
    vector<double> stresslet;
    stresslet.resize(6);

    for(int i = 0; i < nPart; i++){
      fh[0] = 0.0;
      fh[1] = 0.0;
      fh[2] = 0.0;
      
      th[0] = 0.0;
      th[1] = 0.0;
      th[2] = 0.0;

      stresslet[0] = 0.0;
      stresslet[1] = 0.0;
      stresslet[2] = 0.0;
      stresslet[3] = 0.0;
      stresslet[4] = 0.0;
      stresslet[5] = 0.0;

      exchangeParticleData->calculateHydrodynamicInteractions(dynamics, unitConversion, atom->tag[i], atom->x[i], atom->radius[i], fh, th, stresslet);

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

      // Store forces and torque for timesteps when LBM is not called
      hydrodynamicInteractions[i][0] = fh[0];
      hydrodynamicInteractions[i][1] = fh[1];
      hydrodynamicInteractions[i][2] = fh[2];
      hydrodynamicInteractions[i][3] = th[0];
      hydrodynamicInteractions[i][4] = th[1];
      hydrodynamicInteractions[i][5] = th[2];

      double stresslet_arr[6] = {stresslet[0], stresslet[1], stresslet[2], stresslet[3], stresslet[4], stresslet[5]};
      v_tally(i, stresslet_arr);
    }

    // Communicate (and sum) forces from ghost particle back to the original particle
    comm->reverse_comm();
    comm->reverse_comm_fix(this,0);
  }
}


double fix_LBM_PSM::get_rho()
{
    return rho;
}


//==========================================================================
//   allocate atom-based array
//==========================================================================
void fix_LBM_PSM::grow_arrays(int nmax)
{
  memory->grow(hydrodynamicInteractions,nmax,6,"fix_LBM_PSM:hydrodynamicInteractions");
}

//==========================================================================
//   copy values within local atom-based array
//==========================================================================
void fix_LBM_PSM::copy_arrays(int i, int j, int /*delflag*/)
{
  hydrodynamicInteractions[j][0] = hydrodynamicInteractions[i][0];
  hydrodynamicInteractions[j][1] = hydrodynamicInteractions[i][1];
  hydrodynamicInteractions[j][2] = hydrodynamicInteractions[i][2];
  hydrodynamicInteractions[j][3] = hydrodynamicInteractions[i][3];
  hydrodynamicInteractions[j][4] = hydrodynamicInteractions[i][4];
  hydrodynamicInteractions[j][5] = hydrodynamicInteractions[i][5];
}

//==========================================================================
//   pack values in local atom-based array for exchange with another proc
//==========================================================================
int fix_LBM_PSM::pack_exchange(int i, double *buf)
{
  buf[0] = hydrodynamicInteractions[i][0];
  buf[1] = hydrodynamicInteractions[i][1];
  buf[2] = hydrodynamicInteractions[i][2];
  buf[3] = hydrodynamicInteractions[i][3];
  buf[4] = hydrodynamicInteractions[i][4];
  buf[5] = hydrodynamicInteractions[i][5];

  return 6;
}


int fix_LBM_PSM::pack_reverse_comm_size(int n, int first)
{
  int i,last;

  int m = 0;
  last = first + n;

  for (i = first; i < last; i++)
    m += 6;

  return m;
}


//==========================================================================
//   unpack values in local atom-based array from exchange with another proc
//==========================================================================
int fix_LBM_PSM::unpack_exchange(int nlocal, double *buf)
{
  hydrodynamicInteractions[nlocal][0] = buf[0];
  hydrodynamicInteractions[nlocal][1] = buf[1];
  hydrodynamicInteractions[nlocal][2] = buf[2];
  hydrodynamicInteractions[nlocal][3] = buf[3];
  hydrodynamicInteractions[nlocal][4] = buf[4];
  hydrodynamicInteractions[nlocal][5] = buf[5];

  return 6;
}


int fix_LBM_PSM::pack_reverse_comm(int n, int first, double *buf)
{
  int i,k,last;

  int m = 0;
  last = first + n;

  for (i = first; i < last; i++) {
    for (k = 0; k < 6; k++){
      buf[m++] = hydrodynamicInteractions[i][k];
    }
  }

  return m;
}


void fix_LBM_PSM::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,k,kk,ncount;

  int m = 0;

  for (i = 0; i < n; i++) {
    for (k = 0; k < 6; k++){
      j = list[i];
      hydrodynamicInteractions[j][k] += buf[m++];
    }
  }
}
