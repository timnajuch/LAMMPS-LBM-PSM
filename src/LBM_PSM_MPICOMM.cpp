/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#include "LBM_PSM_MPICOMM.h"

LBMPSMMPI::LBMPSMMPI(MPI_Comm world_, int decomposition[3], int procNeigh[6], int procCoordinates_[3], int dimension_)
{
    world = world_;
    MPI_Comm_size(world, &size);
    MPI_Comm_rank(world, &rank);

    dimension = dimension_;
    if(dimension == 2){
      q = 9;
    }else{
      q = 19;
    }

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

    procCoordinates[0] = procCoordinates_[0];
    procCoordinates[1] = procCoordinates_[1];
    procCoordinates[2] = procCoordinates_[2];

}

LBMPSMMPI::~LBMPSMMPI() {}
