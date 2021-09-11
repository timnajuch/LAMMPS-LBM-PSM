#include "PSM_LBM_MPICOMM.h"

PSM_LBM_MPI::PSM_LBM_MPI(MPI_Comm world_, int decomposition[3], int procNeigh[6])
{
//MPICommunication::MPICommunication(int *argc, char ***argv){
//    MPI_Init(argc, argv);
    //MPI_Comm_size(MPI_COMM_WORLD, &size);
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    world = world_;
    MPI_Comm_size(world, &size);
    MPI_Comm_rank(world, &rank);
    int periodicity[3] = {1, 1, 1};
    dimensions[0] = decomposition[0];
    dimensions[1] = decomposition[1];
    dimensions[2] = decomposition[2];
    westRank  = procNeigh[0];
    eastRank  = procNeigh[1];
    southRank = procNeigh[2];
    northRank = procNeigh[3];
    downRank  = procNeigh[4];
    upRank    = procNeigh[5];
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

//  domainDecomposition(decomposition, procNeigh);

}

//MPICommunication::~MPICommunication(){}
PSM_LBM_MPI::~PSM_LBM_MPI(){}

/*
//void MPICommunication::domainDecomposition(int decomposition[3]){
void PSM_LBM_MPI::domainDecomposition(int decomposition[3], int **procNeigh){
    int periodicity[3] = {1, 1, 1};
    dimensions[0] = decomposition[0];
    dimensions[1] = decomposition[1];
    dimensions[2] = decomposition[2];
    westRank  = procNeigh[0][0];
    eastRank  = procNeigh[0][1];
    southRank = procNeigh[1][0];
    northRank = procNeigh[1][1];
    downRank  = procNeigh[2][0];
    upRank    = procNeigh[2][1];
}
*/

/*
void fix_PSM_LBM_MPI::pre_force(int)
{
  int envelopeWidth = 1;
  //sendRecvData<double>(fixPSMLBM->dynamics->getVector_f(), false, 0, fixPSMLBM->dynamics->get_nx(), fixPSMLBM->dynamics->get_ny(), 1, envelopeWidth, true);
  sendRecvData<double>(fixPSMLBM->dynamics->getVector_f(), false, 0, fixPSMLBM->dynamics->get_nx(), fixPSMLBM->dynamics->get_ny(), 1, envelopeWidth, false);
}
*/
