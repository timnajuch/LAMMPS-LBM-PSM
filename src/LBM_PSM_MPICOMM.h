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

#include <algorithm>
#include <iostream>
#include <vector>
#include "mpi.h"
#include <typeinfo>

using namespace std;

template<typename T> inline MPI_Datatype get_mpi_type();
template<> inline MPI_Datatype get_mpi_type<double>() { return MPI_DOUBLE; }
template<> inline MPI_Datatype get_mpi_type<float>()  { return MPI_FLOAT; }
template<> inline MPI_Datatype get_mpi_type<int>()    { return MPI_INT; }


class LBMPSMMPI{

    private:
        int dimension;
        int q;

        // Pre-allocated raw byte buffers — avoids heap allocation in the hot path.
        // Memory is reserved in setupBuffers() for the largest face × envelopeWidth × q.
        // In sendRecvData we resize (O(1) within capacity) and reinterpret as T*.
        std::vector<char> sendBuf_m1, recvBuf_m1;
        std::vector<char> sendBuf_m2, recvBuf_m2;

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

        void setupBuffers(int nx, int ny, int nz, int envelopeWidth, int q);

        template<typename T> void sendRecvData(vector<T> &data_, int commDirection, int nx, int ny, int nz, int envelopeWidth, bool isPeriodic);
        template<typename T> void packData(T* sendBuf, const T* data, const int direction[3],
                                           int envelopeStart, int dataSize, int nx, int ny, int nz, int envelopeWidth);
        template<typename T> void unpackData(const T* recvBuf, T* data, const int direction[3],
                                           int envelopeStart, int dataSize, int nx, int ny, int nz, int envelopeWidth);

};



inline void LBMPSMMPI::setupBuffers(int nx, int ny, int nz, int envelopeWidth, int q) {
  size_t maxSide  = std::max({nx*ny, nx*nz, ny*nz});
  // Reserve raw bytes; sizeof(double) is the widest type ever stored here.
  size_t maxBytes = maxSide * envelopeWidth * q * sizeof(double);
  sendBuf_m1.reserve(maxBytes);
  recvBuf_m1.reserve(maxBytes);
  sendBuf_m2.reserve(maxBytes);
  recvBuf_m2.reserve(maxBytes);
}


template<typename T> void LBMPSMMPI::sendRecvData(vector<T> &data_, int commDirection, int nx, int ny, int nz, int envelopeWidth, bool isPeriodic)
{
  const int dataSize = q;
  int commDataSize = 0;
  int direction[3] = {0, 0, 0};
  const int sendRank[3][2] = { {westRank, eastRank},
                               {southRank, northRank},
                               {downRank,  upRank} };
  const int recvRank[3][2] = { {eastRank,  westRank},
                               {northRank, southRank},
                               {upRank,    downRank} };

  MPI_Datatype commDataType = get_mpi_type<T>();
  switch(commDirection)
  {
    case 0: // x-direction
      commDataSize = dataSize * envelopeWidth * (dimension == 2 ? ny : ny*nz);
      direction[0] = 1;
      break;
    case 1: // y-direction
      commDataSize = dataSize * envelopeWidth * (dimension == 2 ? nx : nx*nz);
      direction[1] = 1;
      break;
    case 2: // z-direction (3D only)
      commDataSize = dataSize * nx * ny * envelopeWidth;
      direction[2] = 1;
      break;
  }

  // Resize pre-reserved member buffers — O(1), no heap allocation.
  sendBuf_m1.resize(commDataSize * sizeof(T));
  recvBuf_m1.resize(commDataSize * sizeof(T));
  sendBuf_m2.resize(commDataSize * sizeof(T));
  recvBuf_m2.resize(commDataSize * sizeof(T));
  T* sBuf1 = reinterpret_cast<T*>(sendBuf_m1.data());
  T* rBuf1 = reinterpret_cast<T*>(recvBuf_m1.data());
  T* sBuf2 = reinterpret_cast<T*>(sendBuf_m2.data());
  T* rBuf2 = reinterpret_cast<T*>(recvBuf_m2.data());

  // First exchange: inner real cells → negative neighbour, ghost layer from positive neighbour
  int envelopeStart = direction[0]*(nx - envelopeWidth*2)
                    + direction[1]*(ny - envelopeWidth*2)
                    + direction[2]*(nz - envelopeWidth*2);

  packData(sBuf1, data_.data(), direction, envelopeStart, dataSize, nx, ny, nz, envelopeWidth);

  MPI_Sendrecv(sBuf1, commDataSize, commDataType, recvRank[commDirection][0], 0,
               rBuf1, commDataSize, commDataType, sendRank[commDirection][0], 0,
               world, &status);

  if(procCoordinates[commDirection] != 0 || isPeriodic)
    unpackData(rBuf1, data_.data(), direction, 0, dataSize, nx, ny, nz, envelopeWidth);

  // Second exchange: inner real cells → positive neighbour, ghost layer from negative neighbour
  // envelopeStart = direction[d]*envelopeWidth; since exactly one direction[d]==1 this equals envelopeWidth.
  packData(sBuf2, data_.data(), direction, envelopeWidth, dataSize, nx, ny, nz, envelopeWidth);

  MPI_Sendrecv(sBuf2, commDataSize, commDataType, recvRank[commDirection][1], 0,
               rBuf2, commDataSize, commDataType, sendRank[commDirection][1], 0,
               world, &status);

  if(procCoordinates[commDirection] != dimensions[commDirection]-1 || isPeriodic)
  {
    envelopeStart = direction[0]*(nx - envelopeWidth)
                  + direction[1]*(ny - envelopeWidth)
                  + direction[2]*(nz - envelopeWidth);
    unpackData(rBuf2, data_.data(), direction, envelopeStart, dataSize, nx, ny, nz, envelopeWidth);
  }
}


template<typename T> void LBMPSMMPI::packData(T* sendBuf, const T* data, const int direction[3],
                                               int envelopeStart, int dataSize,
                                               int nx, int ny, int nz, int envelopeWidth)
{
  const int iMax = nx*(1-direction[0]) + direction[0];
  const int jMax = ny*(1-direction[1]) + direction[1];
  const int kMax = nz*(1-direction[2]) + direction[2];
  const int jMaxkMax     = jMax * kMax;
  const int iMaxjMaxkMax = iMax * jMaxkMax;
  const int nynz         = ny * nz;

  // Per-direction strides into the lattice data array (in units of T elements).
  // Exactly one direction[d] == 1; the others are 0, so these collapse cleanly.
  const int envDataStride = (direction[0]*nynz + direction[1]*nz + direction[2]) * dataSize;
  const int iDataStride   = (1-direction[0]) * nynz * dataSize;
  const int jDataStride   = (1-direction[1]) * nz   * dataSize;
  const int kDataStride   = (1-direction[2])         * dataSize;

  for (int iEnvelope = 0; iEnvelope < envelopeWidth; ++iEnvelope) {
    const int dataEnvBase = (envelopeStart + iEnvelope) * envDataStride;
    const int bufEnvBase  = iEnvelope * iMaxjMaxkMax * dataSize;
    for (int i = 0; i < iMax; ++i) {
      const int dataIBase = dataEnvBase + i * iDataStride;
      const int bufIBase  = bufEnvBase  + i * jMaxkMax * dataSize;
      for (int j = 0; j < jMax; ++j) {
        const int dataJBase = dataIBase + j * jDataStride;
        const int bufJBase  = bufIBase  + j * kMax * dataSize;
        for (int k = 0; k < kMax; ++k) {
          std::copy_n(data   + dataJBase + k * kDataStride,
                      dataSize,
                      sendBuf + bufJBase  + k * dataSize);
        }
      }
    }
  }
}


template<typename T> void LBMPSMMPI::unpackData(const T* recvBuf, T* data, const int direction[3],
                                                 int envelopeStart, int dataSize,
                                                 int nx, int ny, int nz, int envelopeWidth)
{
  const int iMax = nx*(1-direction[0]) + direction[0];
  const int jMax = ny*(1-direction[1]) + direction[1];
  const int kMax = nz*(1-direction[2]) + direction[2];
  const int jMaxkMax     = jMax * kMax;
  const int iMaxjMaxkMax = iMax * jMaxkMax;
  const int nynz         = ny * nz;

  const int envDataStride = (direction[0]*nynz + direction[1]*nz + direction[2]) * dataSize;
  const int iDataStride   = (1-direction[0]) * nynz * dataSize;
  const int jDataStride   = (1-direction[1]) * nz   * dataSize;
  const int kDataStride   = (1-direction[2])         * dataSize;

  for (int iEnvelope = 0; iEnvelope < envelopeWidth; ++iEnvelope) {
    const int dataEnvBase = (envelopeStart + iEnvelope) * envDataStride;
    const int bufEnvBase  = iEnvelope * iMaxjMaxkMax * dataSize;
    for (int i = 0; i < iMax; ++i) {
      const int dataIBase = dataEnvBase + i * iDataStride;
      const int bufIBase  = bufEnvBase  + i * jMaxkMax * dataSize;
      for (int j = 0; j < jMax; ++j) {
        const int dataJBase = dataIBase + j * jDataStride;
        const int bufJBase  = bufIBase  + j * kMax * dataSize;
        for (int k = 0; k < kMax; ++k) {
          std::copy_n(recvBuf + bufJBase  + k * dataSize,
                      dataSize,
                      data    + dataJBase + k * kDataStride);
        }
      }
    }
  }
}

#endif
