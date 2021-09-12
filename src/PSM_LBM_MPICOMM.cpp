#include "PSM_LBM_MPICOMM.h"

PSM_LBM_MPI::PSM_LBM_MPI(MPI_Comm world_, int decomposition[3], int procNeigh[6])
{
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
}

PSM_LBM_MPI::~PSM_LBM_MPI() {}
