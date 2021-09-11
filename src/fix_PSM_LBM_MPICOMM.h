/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

//#ifndef MPI_COMMUNICATION_H
//#define MPI_COMMUNICATION_H
#ifdef FIX_CLASS

FixStyle(lbm-psm-mpi,fix_PSM_LBM_MPI)

#else

#ifndef LMP_FIX_PSM_LBM_MPICOMM_H
#define LMP_FIX_PSM_LBM_MPICOMM_H

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
#include "fix.h"
#include "comm.h"
#include "mpi.h"

#include "fix_PSM_LBM.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;


//namespace LAMMPS_NS {

class fix_PSM_LBM_MPI : public Fix {
//MPICommunication{

    public:
        //MPICommunication(int *argc, char ***argv);
        //~MPICommunication();
        fix_PSM_LBM_MPI(class LAMMPS *, int, char **);
        ~fix_PSM_LBM_MPI();
    
        int setmask();
        void init();
        void pre_force(int);

        int size;
        int rank;
        int dimensions[3];
        int procCoordinates[3];
        int northRank, eastRank, southRank, westRank, upRank, downRank;

//        MPI_Comm cartComm3D;
        MPI_Status status;

        void domainDecomposition(int decomposition[3]);

        void returnProcCoordinatesArray(vector<int>& procCoordinates_){ procCoordinates_[0] = procCoordinates[0]; procCoordinates_[1] = procCoordinates[1]; procCoordinates_[2] = procCoordinates[2]; }

        template<typename T> MPI_Datatype get_type();
        template<typename T> void sendRecvData(vector<T> &data_, bool isVector3D, int commDirection, int nx, int ny, int nz, int envelopeWidth, bool periodicInX);
        template<typename T> void packData(vector<T> &sendBuf, vector<T> &data, int direction[3], int envelopeIterSend[3],
                                           int envlopeStart, int dataSize, int nx, int ny, int nz, int envelopeWidth);
        template<typename T> void unpackData(vector<T> &recvBuf, vector<T> &data, int direction[3], int envelopeIterRecv[3],
                                           int envelopeStart, int dataSize, int nx, int ny, int nz, int envelopeWidth);


        class fix_PSM_LBM *fixPSMLBM;
};



template<typename T> MPI_Datatype fix_PSM_LBM_MPI::get_type()
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


template<typename T> void fix_PSM_LBM_MPI::sendRecvData(vector<T> &data_, bool isVector3D, int commDirection, int nx, int ny, int nz, int envelopeWidth, bool periodicInX)
{
  //int dataSize = 1+2*(int)isVector3D;
  int dataSize = 9;
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
//cout << "DEBUG SENDRECV. A" << endl;
  switch(commDirection)
  {
    case 0: // x-direction
      commDataSize = dataSize*ny*envelopeWidth;
      direction[0] = 1;
      break;
    case 1: // y-direction
      commDataSize = dataSize*nx*envelopeWidth;
      direction[1] = 1;
      break;
    //case 2: // z-direction // not 3D yet
  }

  int envelopeIterSend[3];
  int envelopeIterRecv[3];
  envelopeIterSend[0] = direction[0]*(nx-envelopeWidth*2);
  envelopeIterSend[1] = direction[1]*(ny-envelopeWidth*2);
  envelopeIterSend[2] = direction[2]*(nz-envelopeWidth*2);
  envelopeIterRecv[0] = 0;
  envelopeIterRecv[1] = 0;
  envelopeIterRecv[2] = 0;
//int a = 1;
//int b = 2;
  vector<T> sendBuf1;
  vector<T> recvBuf1;
  //commDataSize = 600;
  sendBuf1.resize(commDataSize);
  recvBuf1.resize(commDataSize);

  envelopeStart = direction[0]*envelopeIterSend[0] +
                  direction[1]*envelopeIterSend[1] +
                  direction[2]*envelopeIterSend[2];
//cout << "DEBUG SENDRECV. B" << endl;
  packData(sendBuf1, data_, direction, envelopeIterSend, envelopeStart, dataSize, nx, ny, nz, envelopeWidth);
//cout << "DEBUG SENDRECV. C: " << commDataSize << " / " << commDirection << " / " << rank << " / " << recvRank[0][0] << " / " << sendRank[0][0] << " / " << westRank << " / " << eastRank << endl;
//for(int i =0; i<commDataSize; ++i){
//    sendBuf1[i] = 0.0;
//    recvBuf1[i] = 0.0;
//    //cout << "DEBUG  sendbuf and rcbuf: " << i << " / " << sendBuf1[i] << " / " << recvBuf1[i] << endl;
//}
//cout << "DEBUG SENDRECV. Cb" << endl;
  MPI_Sendrecv(&sendBuf1[0], commDataSize, commDataType, recvRank[commDirection][0], 0,
               &recvBuf1[0], commDataSize, commDataType, sendRank[commDirection][0], 0,
               world, &status);
  //MPI_Sendrecv(&sendBuf1[0], commDataSize, MPI_DOUBLE, recvRank[commDirection][0], 0,
   //            &recvBuf1[0], commDataSize, MPI_DOUBLE, sendRank[commDirection][0], 0,
   //            world, &status);
  //MPI_Sendrecv(&a, 1, MPI_INT, recvRank[commDirection][0], 0,
  //            &b, 1, MPI_INT, sendRank[commDirection][0], 0,
  //            world, &status);
  //MPI_Barrier(world);
//cout << "DEBUG SENDRECV. D" << endl;
  if(procCoordinates[commDirection] != 0 || (periodicInX == true && commDirection == 0))
  {
    envelopeStart = direction[0]*envelopeIterRecv[0] +
                    direction[1]*envelopeIterRecv[1] +
                    direction[2]*envelopeIterRecv[2];
//cout << "DEBUG SENDRECV. E" << endl;
    unpackData(recvBuf1, data_, direction, envelopeIterRecv, envelopeStart, dataSize, nx, ny, nz, envelopeWidth);
  }
//cout << "DEBUG SENDRECV. F" << endl;
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

  if(procCoordinates[commDirection] != dimensions[commDirection]-1 || (periodicInX == true && commDirection == 0))
  {
    envelopeStart = direction[0]*envelopeIterRecv[0] +
                    direction[1]*envelopeIterRecv[1] +
                    direction[2]*envelopeIterRecv[2];

    unpackData(recvBuf2, data_, direction, envelopeIterRecv, envelopeStart, dataSize, nx, ny, nz, envelopeWidth); 
  }

}


template<typename T> void fix_PSM_LBM_MPI::packData(vector<T> &sendBuf, vector<T> &data, int direction[3], int envelopeIterSend[3], int envelopeStart, int dataSize, int nx, int ny, int nz, int envelopeWidth)
{
  int iMax = nx*(1-direction[0])+direction[0];
  int jMax = ny*(1-direction[1])+direction[1];
  int kMax = nz*(1-direction[2])+direction[2];  

  int envelopeIterSendIndex = 0;
  int sendBufIndex = 0;

  int q = 9;

  for(int iEnvelope = 0; iEnvelope < envelopeWidth; ++iEnvelope){
    for(int i = 0; i < iMax; ++i){
      for(int j = 0; j < jMax; ++j){
        for(int k = 0; k < kMax; ++k){
          for(int iq = 0; iq < q; ++iq){
            envelopeIterSend[0] = direction[0]*(envelopeStart + iEnvelope) + i;
            envelopeIterSend[1] = direction[1]*(envelopeStart + iEnvelope) + j;
            envelopeIterSend[2] = direction[2]*(envelopeStart + iEnvelope) + k;

            envelopeIterSendIndex = envelopeIterSend[0] * ny * q + envelopeIterSend[1] * q + iq;
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


template<typename T> void fix_PSM_LBM_MPI::unpackData(vector<T> &recvBuf, vector<T> &data, int direction[3], int envelopeIterRecv[3], int envelopeStart, int dataSize, int nx, int ny, int nz, int envelopeWidth)
{
  int iMax = nx*(1-direction[0])+direction[0];
  int jMax = ny*(1-direction[1])+direction[1];
  int kMax = nz*(1-direction[2])+direction[2];

  int envelopeIterRecvIndex = 0;
  int recvBufIndex = 0;

  int q = 9;

  for(int iEnvelope = 0; iEnvelope < envelopeWidth; ++iEnvelope){
    for(int i = 0; i < iMax; ++i){
      for(int j = 0; j < jMax; ++j){
        for(int k = 0; k < kMax; ++k){
          for(int iq = 0; iq < q; ++iq){
            envelopeIterRecv[0] = direction[0]*(envelopeStart + iEnvelope) + i;
            envelopeIterRecv[1] = direction[1]*(envelopeStart + iEnvelope) + j;
            envelopeIterRecv[2] = direction[2]*(envelopeStart + iEnvelope) + k;

            //envelopeIterRecvIndex = i * ny * q + j * q + iq;
            envelopeIterRecvIndex = envelopeIterRecv[0] * ny * q + envelopeIterRecv[1] * q + iq;
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

//}
#endif
#endif
