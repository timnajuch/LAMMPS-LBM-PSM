//#include "mpiCommunication.h"
#include "fix_PSM_LBM_MPICOMM.h"

fix_PSM_LBM_MPI::fix_PSM_LBM_MPI(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
//MPICommunication::MPICommunication(int *argc, char ***argv){
//    MPI_Init(argc, argv);
    //MPI_Comm_size(MPI_COMM_WORLD, &size);
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(world, &size);
    MPI_Comm_rank(world, &rank);

//    cartComm3D = (MPI_Comm)NULL;
/*
    procCoordinates[0] = 0;
    procCoordinates[1] = 0;
    procCoordinates[2] = 0;
    northRank = 0;
    eastRank = 0;
    southRank = 0;
    westRank = 0;
    upRank = 0;
    downRank = 0;
*/

  int decomposition[3] = {comm->procgrid[0], comm->procgrid[1], comm->procgrid[2]};
  domainDecomposition(decomposition);

  for(int ifix=0; ifix<modify->nfix; ifix++)
    if(strcmp(modify->fix[ifix]->style,"lbm-psm")==0)
      fixPSMLBM = (fix_PSM_LBM *)modify->fix[ifix];
}

//MPICommunication::~MPICommunication(){}
fix_PSM_LBM_MPI::~fix_PSM_LBM_MPI(){}

int fix_PSM_LBM_MPI::setmask()
{
  int mask =0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

void fix_PSM_LBM_MPI::init()
{
}

//void MPICommunication::domainDecomposition(int decomposition[3]){
void fix_PSM_LBM_MPI::domainDecomposition(int decomposition[3]){
    int periodicity[3] = {1, 1, 1};
    dimensions[0] = decomposition[0];
    dimensions[1] = decomposition[1];
    dimensions[2] = decomposition[2];
/*
    MPI_Cart_create(MPI_COMM_WORLD, 3, decomposition, periodicity, false, &cartComm3D);

    MPI_Comm_rank(cartComm3D, &rank);

    MPI_Cart_coords(cartComm3D, rank, 3, procCoordinates);
    MPI_Cart_shift(cartComm3D, 0, 1, &westRank, &eastRank);
    MPI_Cart_shift(cartComm3D, 1, 1, &southRank, &northRank);
    MPI_Cart_shift(cartComm3D, 2, 1, &downRank, &upRank);
*/
    westRank  = comm->procneigh[0][0];
    eastRank  = comm->procneigh[0][1];
    southRank = comm->procneigh[1][0];
    northRank = comm->procneigh[1][1];
    downRank  = comm->procneigh[2][0];
    upRank    = comm->procneigh[2][1];
}


void fix_PSM_LBM_MPI::pre_force(int)
{
  int envelopeWidth = 1;
  //sendRecvData<double>(fixPSMLBM->dynamics->getVector_f(), false, 0, fixPSMLBM->dynamics->get_nx(), fixPSMLBM->dynamics->get_ny(), 1, envelopeWidth, true);
  sendRecvData<double>(fixPSMLBM->dynamics->getVector_f(), false, 0, fixPSMLBM->dynamics->get_nx(), fixPSMLBM->dynamics->get_ny(), 1, envelopeWidth, false);
}
