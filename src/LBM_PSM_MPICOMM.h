/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#ifndef LBM_PSM_MPICOMM_H
#define LBM_PSM_MPICOMM_H

#include <iostream>
#include <vector>
#include "mpi.h"
#include <typeinfo>

using namespace std;

class LBMPSMMPI{

    private:
        int dimension;
        int q;

    public:
        LBMPSMMPI(MPI_Comm world_, int decomposition[3], int procNeigh[6], int procCoordinates_[3], int dimension_);
        ~LBMPSMMPI();
    
        int size;
        int rank;

        int dimensions[3];
        int procCoordinates[3];
        int northRank, eastRank, southRank, westRank, upRank, downRank;

        MPI_Comm world;
        MPI_Status status;

        template<typename T> MPI_Datatype get_type();
        template<typename T> void sendRecvData(vector<T> &data_, int commDirection, int nx, int ny, int nz, int envelopeWidth, bool isPeriodic, int currentStep);
        template<typename T> void packData(vector<T> &sendBuf, vector<T> &data, int direction[3], int envelopeIterSend[3],
                                           int envlopeStart, int dataSize, int nx, int ny, int nz, int envelopeWidth, int currentStep);
        template<typename T> void unpackData(vector<T> &recvBuf, vector<T> &data, int direction[3], int envelopeIterRecv[3],
                                           int envelopeStart, int dataSize, int nx, int ny, int nz, int envelopeWidth, int currentStep);

};



template<typename T> MPI_Datatype LBMPSMMPI::get_type()
{
  char name = typeid(T).name()[0];
  MPI_Datatype returnMPIDataType = MPI_DOUBLE; // Default MPI_Datatype (surpresses compiler warnings)
  switch (name)
  {
    case 'i':
      returnMPIDataType = MPI_INT;
      break;
    case 'f':
      returnMPIDataType = MPI_FLOAT;
      break;
    case 'j':
      returnMPIDataType = MPI_UNSIGNED;
      break;
    case 'd':
      returnMPIDataType = MPI_DOUBLE;
      break;
    case 'c':
      returnMPIDataType = MPI_CHAR;
      break;
    case 's':
      returnMPIDataType = MPI_SHORT;
      break;
    case 'l':
      returnMPIDataType = MPI_LONG;
      break;
    case 'm':
      returnMPIDataType = MPI_UNSIGNED_LONG;
      break;
    case 'b':
      returnMPIDataType = MPI_BYTE;
      break;
    default:
      std::cout << "MPI_DataTyppe is not defined. Check the data which is sent in the LBM lattice class." << std::endl;
      // todo die function
      break;
  }
  return returnMPIDataType;
}


template<typename T> void LBMPSMMPI::sendRecvData(vector<T> &data_, int commDirection, int nx, int ny, int nz, int envelopeWidth, bool isPeriodic, int currentStep)
{
  int dataSize = q;
  int commDataSize = 0;
  int envelopeStart;
  int direction[3] = {0, 0, 0};
  int sendRank[3][2] = { {westRank, eastRank},
                         {southRank, northRank},
                         {downRank, upRank} };
  int recvRank[3][2] = { {eastRank, westRank},
                         {northRank, southRank},
                         {upRank, downRank} };

  MPI_Datatype commDataType = get_type<T>();
  switch(commDirection)
  {
    case 0: // x-direction
      if(dimension == 2){
        commDataSize = dataSize*ny*envelopeWidth;
      }else{
        commDataSize = dataSize*ny*nz*envelopeWidth;
      }
      direction[0] = 1;
      break;
    case 1: // y-direction
      if(dimension == 2){
        commDataSize = dataSize*nx*envelopeWidth;
      }else{
        commDataSize = dataSize*nx*nz*envelopeWidth;
      }
      direction[1] = 1;
      break;
    case 2: // z-direction // Only for 3D. In 2D this direction is not possible. 2D in x-y-plane
      commDataSize = dataSize*nx*ny*envelopeWidth;
      direction[2] = 1;
      break;
  }

  int envelopeIterSend[3];
  int envelopeIterRecv[3];
  envelopeIterSend[0] = direction[0]*(nx-envelopeWidth*2);
  envelopeIterSend[1] = direction[1]*(ny-envelopeWidth*2);
  envelopeIterSend[2] = direction[2]*(nz-envelopeWidth*2);
  envelopeIterRecv[0] = 0;
  envelopeIterRecv[1] = 0;
  envelopeIterRecv[2] = 0;
  vector<T> sendBuf1;
  vector<T> recvBuf1;
  sendBuf1.resize(commDataSize);
  recvBuf1.resize(commDataSize);

  envelopeStart = direction[0]*envelopeIterSend[0] +
                  direction[1]*envelopeIterSend[1] +
                  direction[2]*envelopeIterSend[2];

  packData(sendBuf1, data_, direction, envelopeIterSend, envelopeStart, dataSize, nx, ny, nz, envelopeWidth, currentStep);

  MPI_Sendrecv(&sendBuf1[0], commDataSize, commDataType, recvRank[commDirection][0], 0,
               &recvBuf1[0], commDataSize, commDataType, sendRank[commDirection][0], 0,
               world, &status);

  if(procCoordinates[commDirection] != 0 || (isPeriodic == true))
  {
    envelopeStart = direction[0]*envelopeIterRecv[0] +
                    direction[1]*envelopeIterRecv[1] +
                    direction[2]*envelopeIterRecv[2];

    unpackData(recvBuf1, data_, direction, envelopeIterRecv, envelopeStart, dataSize, nx, ny, nz, envelopeWidth, currentStep);
  }

  envelopeIterSend[0] = direction[0]*envelopeWidth;
  envelopeIterSend[1] = direction[1]*envelopeWidth;
  envelopeIterSend[2] = direction[2]*envelopeWidth;
  envelopeIterRecv[0] = direction[0]*(nx-envelopeWidth);
  envelopeIterRecv[1] = direction[1]*(ny-envelopeWidth);
  envelopeIterRecv[2] = direction[2]*(nz-envelopeWidth);

  vector<T> sendBuf2;
  vector<T> recvBuf2;
  sendBuf2.resize(commDataSize);
  recvBuf2.resize(commDataSize);

  envelopeStart = direction[0]*envelopeIterSend[0] +
                  direction[1]*envelopeIterSend[1] +
                  direction[2]*envelopeIterSend[2];

  packData(sendBuf2, data_, direction, envelopeIterSend, envelopeStart, dataSize, nx, ny, nz, envelopeWidth, currentStep);

  MPI_Sendrecv(&sendBuf2[0], commDataSize, commDataType, recvRank[commDirection][1], 0,
               &recvBuf2[0], commDataSize, commDataType, sendRank[commDirection][1], 0,
               world, &status);

  if(procCoordinates[commDirection] != dimensions[commDirection]-1 || (isPeriodic == true))
  {
    envelopeStart = direction[0]*envelopeIterRecv[0] +
                    direction[1]*envelopeIterRecv[1] +
                    direction[2]*envelopeIterRecv[2];

    unpackData(recvBuf2, data_, direction, envelopeIterRecv, envelopeStart, dataSize, nx, ny, nz, envelopeWidth, currentStep); 
  }

}


template<typename T> void LBMPSMMPI::packData(vector<T> &sendBuf, vector<T> &data, int direction[3], int envelopeIterSend[3], int envelopeStart, int dataSize, int nx, int ny, int nz, int envelopeWidth, int currentStep)
{
  int iMax = nx*(1-direction[0])+direction[0];
  int jMax = ny*(1-direction[1])+direction[1];
  int kMax = nz*(1-direction[2])+direction[2];  

  int envelopeIterSendIndex = 0;
  int sendBufIndex = 0;

  for(int iEnvelope = 0; iEnvelope < envelopeWidth; ++iEnvelope){
    for(int i = 0; i < iMax; ++i){
      for(int j = 0; j < jMax; ++j){
        for(int k = 0; k < kMax; ++k){
          envelopeIterSend[0] = direction[0]*(envelopeStart + iEnvelope) + i;
          envelopeIterSend[1] = direction[1]*(envelopeStart + iEnvelope) + j;
          envelopeIterSend[2] = direction[2]*(envelopeStart + iEnvelope) + k;
          for(int iq = 0; iq < dataSize; ++iq){
            envelopeIterSendIndex = nx*ny*nz*dataSize*currentStep + (envelopeIterSend[0] * ny * nz + envelopeIterSend[1] * nz + envelopeIterSend[2])*dataSize + iq;
            sendBufIndex = (k + j*kMax + i*jMax*kMax + iEnvelope*iMax*jMax*kMax)*dataSize + iq;

            sendBuf[sendBufIndex] = data[envelopeIterSendIndex];
          }
        }
      }
    }
  }

}


template<typename T> void LBMPSMMPI::unpackData(vector<T> &recvBuf, vector<T> &data, int direction[3], int envelopeIterRecv[3], int envelopeStart, int dataSize, int nx, int ny, int nz, int envelopeWidth, int currentStep)
{
  int iMax = nx*(1-direction[0])+direction[0];
  int jMax = ny*(1-direction[1])+direction[1];
  int kMax = nz*(1-direction[2])+direction[2];

  int envelopeIterRecvIndex = 0;
  int recvBufIndex = 0;

  for(int iEnvelope = 0; iEnvelope < envelopeWidth; ++iEnvelope){
    for(int i = 0; i < iMax; ++i){
      for(int j = 0; j < jMax; ++j){
        for(int k = 0; k < kMax; ++k){
          envelopeIterRecv[0] = direction[0]*(envelopeStart + iEnvelope) + i;
          envelopeIterRecv[1] = direction[1]*(envelopeStart + iEnvelope) + j;
          envelopeIterRecv[2] = direction[2]*(envelopeStart + iEnvelope) + k;
          for(int iq = 0; iq < dataSize; ++iq){
            envelopeIterRecvIndex = nx*ny*nz*dataSize*currentStep + (envelopeIterRecv[0] * ny * nz + envelopeIterRecv[1] * nz + envelopeIterRecv[2])*dataSize + iq;
            recvBufIndex = (k + j*kMax + i*jMax*kMax + iEnvelope*iMax*jMax*kMax)*dataSize + iq;

            data[envelopeIterRecvIndex] = recvBuf[recvBufIndex];
          }
        }
      }
    }
  }

}

#endif
