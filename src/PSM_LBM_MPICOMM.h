/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#ifndef PSM_LBM_MPICOMM_H
#define PSM_LBM_MPICOMM_H

#include <iostream>
#include <vector>
#include "mpi.h"

using namespace std;

class PSM_LBM_MPI{

    private:
        int dimension;
        int q;

    public:
        PSM_LBM_MPI(MPI_Comm world_, int decomposition[3], int procNeigh[6], int procCoordinates_[3], int dimension_);
        ~PSM_LBM_MPI();
    
        int size;
        int rank;

        int dimensions[3];
        int procCoordinates[3];
        int northRank, eastRank, southRank, westRank, upRank, downRank;

        MPI_Comm world;
        MPI_Status status;

//        void returnProcCoordinatesArray(vector<int>& procCoordinates_){ procCoordinates_[0] = procCoordinates[0]; procCoordinates_[1] = procCoordinates[1]; procCoordinates_[2] = procCoordinates[2]; }

        template<typename T> MPI_Datatype get_type();
        template<typename T> void sendRecvData(vector<T> &data_, bool isVector3D, int commDirection, int nx, int ny, int nz, int envelopeWidth, bool periodicInX);
        template<typename T> void packData(vector<T> &sendBuf, vector<T> &data, int direction[3], int envelopeIterSend[3],
                                           int envlopeStart, int dataSize, int nx, int ny, int nz, int envelopeWidth);
        template<typename T> void unpackData(vector<T> &recvBuf, vector<T> &data, int direction[3], int envelopeIterRecv[3],
                                           int envelopeStart, int dataSize, int nx, int ny, int nz, int envelopeWidth);

};



template<typename T> MPI_Datatype PSM_LBM_MPI::get_type()
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


template<typename T> void PSM_LBM_MPI::sendRecvData(vector<T> &data_, bool isVector3D, int commDirection, int nx, int ny, int nz, int envelopeWidth, bool periodicInX)
{
  //int dataSize = 1+2*(int)isVector3D;
  int dataSize = q; //9; // TODO: Extend to 3D
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
  {/*
    case 0: // x-direction
      commDataSize = dataSize*ny*envelopeWidth;
      direction[0] = 1;
      break;
    case 1: // y-direction
      commDataSize = dataSize*nx*envelopeWidth;
      direction[1] = 1;
      break;*/
    //case 2: // z-direction // not 3D yet
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

  packData(sendBuf1, data_, direction, envelopeIterSend, envelopeStart, dataSize, nx, ny, nz, envelopeWidth);

  MPI_Sendrecv(&sendBuf1[0], commDataSize, commDataType, recvRank[commDirection][0], 0,
               &recvBuf1[0], commDataSize, commDataType, sendRank[commDirection][0], 0,
               world, &status);

  //if(procCoordinates[commDirection] != 0 || (periodicInX == true && commDirection == 0))
  if(procCoordinates[commDirection] != 0 || (periodicInX == true))
  {
    envelopeStart = direction[0]*envelopeIterRecv[0] +
                    direction[1]*envelopeIterRecv[1] +
                    direction[2]*envelopeIterRecv[2];

    unpackData(recvBuf1, data_, direction, envelopeIterRecv, envelopeStart, dataSize, nx, ny, nz, envelopeWidth);
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

  packData(sendBuf2, data_, direction, envelopeIterSend, envelopeStart, dataSize, nx, ny, nz, envelopeWidth);

  MPI_Sendrecv(&sendBuf2[0], commDataSize, commDataType, recvRank[commDirection][1], 0,
               &recvBuf2[0], commDataSize, commDataType, sendRank[commDirection][1], 0,
               world, &status);

  //if(procCoordinates[commDirection] != dimensions[commDirection]-1 || (periodicInX == true && commDirection == 0))
  if(procCoordinates[commDirection] != dimensions[commDirection]-1 || (periodicInX == true))
  {
    envelopeStart = direction[0]*envelopeIterRecv[0] +
                    direction[1]*envelopeIterRecv[1] +
                    direction[2]*envelopeIterRecv[2];

    unpackData(recvBuf2, data_, direction, envelopeIterRecv, envelopeStart, dataSize, nx, ny, nz, envelopeWidth); 
  }

}


template<typename T> void PSM_LBM_MPI::packData(vector<T> &sendBuf, vector<T> &data, int direction[3], int envelopeIterSend[3], int envelopeStart, int dataSize, int nx, int ny, int nz, int envelopeWidth)
{
  int iMax = nx*(1-direction[0])+direction[0];
  int jMax = ny*(1-direction[1])+direction[1];
  int kMax = nz*(1-direction[2])+direction[2];  

  int envelopeIterSendIndex = 0;
  int sendBufIndex = 0;

//  int q = 9;
  for(int iEnvelope = 0; iEnvelope < envelopeWidth; ++iEnvelope){
    for(int i = 0; i < iMax; ++i){
      for(int j = 0; j < jMax; ++j){
        for(int k = 0; k < kMax; ++k){
          for(int iq = 0; iq < q; ++iq){
            envelopeIterSend[0] = direction[0]*(envelopeStart + iEnvelope) + i;
            envelopeIterSend[1] = direction[1]*(envelopeStart + iEnvelope) + j;
            envelopeIterSend[2] = direction[2]*(envelopeStart + iEnvelope) + k;

            //envelopeIterSendIndex = envelopeIterSend[0] * ny * q + envelopeIterSend[1] * q + iq;
            envelopeIterSendIndex = envelopeIterSend[0] * ny * nz * q + envelopeIterSend[1] * nz * q + envelopeIterSend[2] * q + iq;
            sendBufIndex = (k + j*kMax + i*jMax*kMax + iEnvelope*iMax*jMax*kMax)*dataSize + iq;

            sendBuf[sendBufIndex] = data[envelopeIterSendIndex];
            //if(dataSize > 1)
            //{
            //  sendBuf[sendBufIndex + 1] = data[envelopIterSendIndex + 1];
            //  sendBuf[sendBufIndex + 2] = data[envelopIterSendIndex + 2];
            //  todo
            //}
          }
        }
      }
    }
  }

}


template<typename T> void PSM_LBM_MPI::unpackData(vector<T> &recvBuf, vector<T> &data, int direction[3], int envelopeIterRecv[3], int envelopeStart, int dataSize, int nx, int ny, int nz, int envelopeWidth)
{
  int iMax = nx*(1-direction[0])+direction[0];
  int jMax = ny*(1-direction[1])+direction[1];
  int kMax = nz*(1-direction[2])+direction[2];

  int envelopeIterRecvIndex = 0;
  int recvBufIndex = 0;

//  int q = 9;

  for(int iEnvelope = 0; iEnvelope < envelopeWidth; ++iEnvelope){
    for(int i = 0; i < iMax; ++i){
      for(int j = 0; j < jMax; ++j){
        for(int k = 0; k < kMax; ++k){
          for(int iq = 0; iq < q; ++iq){
            envelopeIterRecv[0] = direction[0]*(envelopeStart + iEnvelope) + i;
            envelopeIterRecv[1] = direction[1]*(envelopeStart + iEnvelope) + j;
            envelopeIterRecv[2] = direction[2]*(envelopeStart + iEnvelope) + k;

            //envelopeIterRecvIndex = i * ny * q + j * q + iq;
            //envelopeIterRecvIndex = envelopeIterRecv[0] * ny * q + envelopeIterRecv[1] * q + iq;
            envelopeIterRecvIndex = envelopeIterRecv[0] * ny * nz * q + envelopeIterRecv[1] * nz * q + envelopeIterRecv[2] * q + iq;
            recvBufIndex = (k + j*kMax + i*jMax*kMax + iEnvelope*iMax*jMax*kMax)*dataSize + iq;

            data[envelopeIterRecvIndex] = recvBuf[recvBufIndex];
            //if(dataSize > 1)
            //{
            //  data[envelopeIterRecvIndex + 1] = recvBuf[recvBufIndex + 1];
            //  data[envelopeIterRecvIndex + 2] = recvBuf[recvBufIndex + 2];
            //  todo
            //}
          }
        }
      }
    }
  }

}

#endif
