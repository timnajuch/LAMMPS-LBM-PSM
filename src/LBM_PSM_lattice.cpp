/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#include "LBM_PSM_lattice.h"


LBMPSMLattice::LBMPSMLattice(int nx_, int ny_, int nz_, int q_, int decomposition[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_, double dx){
  envelopeWidth = 1;
  dimension = dimension_;

  procCoordinates[0] = procCoordinates_[0];
  procCoordinates[1] = procCoordinates_[1];
  procCoordinates[2] = procCoordinates_[2];

  procLength.resize(3);
  procOrigin.resize(3);

  c = 1.0;
  cs = 1.0/sqrt(3.0);
  csPow2 = 1.0/3.0;
  nx = nx_/decomposition[0] + envelopeWidth*2;
  ny = ny_/decomposition[1] + envelopeWidth*2;
  nz = 1;
  q = 9;
  if (dimension == 3)
  { 
    nz = nz_/decomposition[2] + envelopeWidth*2;
    q = 19;
  }

  procOrigin[0] = origin_[0] + procCoordinates[0]*(nx-2*envelopeWidth)*dx;
  procOrigin[1] = origin_[1] + procCoordinates[1]*(ny-2*envelopeWidth)*dx;
  procOrigin[2] = origin_[2] + procCoordinates[2]*(nz-2*envelopeWidth)*dx;

  f = vector<double>(nx*ny*nz*q,0.0);
  f0 = vector<double>(nx*ny*nz*q,0.0);
  fcoll = vector<double>(nx*ny*nz*q,0.0);

  B = vector<double>(nx*ny*nz,0.0);
  rho = vector<double>(nx*ny*nz,0.0);
  x = vector<double>(nx*ny*nz,0.0);
  y = vector<double>(nx*ny*nz,0.0);
  z = vector<double>(nx*ny*nz,0.0);
  
  u = vector<double>(nx*ny*nz*3,0.0);
  us = vector<double>(nx*ny*nz*3,0.0);

  Fhydx = vector<double>(nx*ny*nz,0.0);
  Fhydy = vector<double>(nx*ny*nz,0.0);
  Fhydz = vector<double>(nx*ny*nz,0.0);

  pData.resize(nx*ny*nz);

  if(q == 9){
    e = { 0.0 ,  0.0,  0.0,
          1.0 ,  0.0,  0.0,
          0.0 ,  1.0,  0.0,
          -1.0,  0.0,  0.0,
          0.0 , -1.0,  0.0,
          1.0 ,  1.0,  0.0,
          -1.0,  1.0,  0.0,
          -1.0, -1.0,  0.0,
          1.0 ,  -1.0, 0.0
        };

    w = { 4.0/9.0,
          1.0/9.0,
          1.0/9.0,
          1.0/9.0,
          1.0/9.0,
          1.0/36.0,
          1.0/36.0,
          1.0/36.0,
          1.0/36.0
        };
    
    g = { 0.0,
          1.0/3.0,
          1.0/3.0,
          1.0/3.0,
          1.0/3.0,
          1.0/12.0,
          1.0/12.0,
          1.0/12.0,
          1.0/12.0
        };
    }

  if(q == 19){
    e = { 0.0,  0.0,  0.0,
          1.0,  0.0,  0.0,
         -1.0,  0.0,  0.0,
          0.0,  1.0,  0.0,
          0.0, -1.0,  0.0,
          0.0,  0.0,  1.0,
          0.0,  0.0, -1.0,
          1.0,  1.0,  0.0,
         -1.0, -1.0,  0.0,
          1.0,  0.0,  1.0, 
         -1.0,  0.0, -1.0,
          0.0,  1.0,  1.0,
          0.0, -1.0, -1.0,
          1.0, -1.0,  0.0,
         -1.0,  1.0,  0.0,
          1.0,  0.0, -1.0,
         -1.0,  0.0,  1.0,
          0.0,  1.0, -1.0,
          0.0, -1.0,  1.0      
        };

    w = { 1.0/3.0,
          1.0/18.0,
          1.0/18.0,
          1.0/18.0,
          1.0/18.0,
          1.0/18.0,
          1.0/18.0,
          1.0/36.0,
          1.0/36.0,
          1.0/36.0,
          1.0/36.0,
          1.0/36.0,
          1.0/36.0,
          1.0/36.0,
          1.0/36.0,
          1.0/36.0,
          1.0/36.0,
          1.0/36.0,
          1.0/36.0
        };
  }

};


LBMPSMLattice::~LBMPSMLattice(){};


void LBMPSMLattice::initialise_channel_geometry(double wallHeightHalf_, double eps_, double nychannel_, double dx_, double dy_, double dz_){
  dx = dx_;
  dy = dy_;
  dz = dz_;
  
  for(int i = 0; i < nx; ++i){
    for(int j = 0; j < ny; ++j){
      for(int k = 0; k < nz; ++k){
// todo exten the stuf within the loops to 3D 
        int ind_2D = i*ny + j;

        x[ind_2D] = i*dx - dx*envelopeWidth + procCoordinates[0]*(nx-2*envelopeWidth)*dx;
        y[ind_2D] = -(wallHeightHalf_-1.0+0.5+1.0*eps_)*dy+j*dy + procCoordinates[1]*(ny-2*envelopeWidth)*dy - dy*envelopeWidth;

        // Walls placed on horizontal boundaries
        if((j < (int)(0.0+wallHeightHalf_)) || (j > (int)(ny-1.0-wallHeightHalf_))){
          B[ind_2D] = 1.0;
        }else if((j == (int)(0.0+wallHeightHalf_)) || (j == (int)(ny-1.0-wallHeightHalf_))){
          B[ind_2D] = eps_;
        }else{
          B[ind_2D] = 0.0;
        }

        Fhydx[ind_2D] = 0.0;
        Fhydy[ind_2D] = 0.0;
        //rho[ind_2D] = 1.0;
     
      }
    }
  }
};


void LBMPSMLattice::initialise_domain(double dx_, double dy_, double dz_){
  dx = dx_;
  dy = dy_;
  dz = dz_;
  
  for(int i = 0; i < nx; ++i){
    for(int j = 0; j < ny; ++j){
      for(int k = 0; k < nz; ++k){

        int index = i*ny*nz + j*nz + k;

        x[index] = i*dx - dx*envelopeWidth + procCoordinates[0]*(nx-2*envelopeWidth)*dx;
        y[index] = j*dy - dy*envelopeWidth + procCoordinates[1]*(ny-2*envelopeWidth)*dy;
        z[index] = k*dz - dz*envelopeWidth + procCoordinates[2]*(nz-2*envelopeWidth)*dz;


        Fhydx[index] = 0.0;
        Fhydy[index] = 0.0;
        Fhydz[index] = 0.0;
      }
    }
  }

};


void LBMPSMLattice::set_f(int i_, int j_, int k_, int iq_, double value_){
  f[i_*ny*nz*q + j_*nz*q + k_*q + iq_] = value_;
};

void LBMPSMLattice::set_f(int ind_iq_, double value_){
  f[ind_iq_] = value_;
};

void LBMPSMLattice::cp_fcoll_f(){
  f = fcoll;
};

double LBMPSMLattice::get_f(int i_, int j_, int k_, int iq_){
  return f[i_*ny*nz*q + j_*nz*q + k_*q + iq_];
};

double LBMPSMLattice::get_f(int ind_iq_){
  return f[ind_iq_];
};

 void LBMPSMLattice::set_f0(int i_, int j_, int k_, int iq_, double value_){
  f0[i_*ny*nz*q + j_*nz*q + k_*q + iq_] = value_;
};

double LBMPSMLattice::get_f0(int i_, int j_, int k_, int iq_){
  return f0[i_*ny*nz*q + j_*nz*q + k_*q + iq_];
};

double LBMPSMLattice::get_f0(int ind_iq_){
  return f0[ind_iq_];
};

void LBMPSMLattice::set_fcoll(int i_, int j_, int k_, int iq_, double value_){
  fcoll[i_*ny*nz*q + j_*nz*q + k_*q + iq_] = value_;
};

double LBMPSMLattice::get_fcoll(int i_, int j_, int k_, int iq_){
  return fcoll[i_*ny*nz*q + j_*nz*q + k_*q + iq_];
};

double LBMPSMLattice::get_fcoll(int ind_iq_){
  return fcoll[ind_iq_];
};

 vector<double> LBMPSMLattice::get_B(){
  return B;
};

 vector<double> LBMPSMLattice::get_rho(){
  return rho;
};

 vector<double> LBMPSMLattice::get_x(){
  return x;
};

 vector<double> LBMPSMLattice::get_y(){
  return y;
};

 vector<double> LBMPSMLattice::get_z(){
  return z;
};

 vector<double> LBMPSMLattice::get_u(){
  return u;
};

int LBMPSMLattice::get_nx(){
  return nx;
};

int LBMPSMLattice::get_ny(){
  return ny;
};

int LBMPSMLattice::get_nz(){
  return nz;
};

void LBMPSMLattice::set_B(int index, double B_){
  B[index] = B_;
};

double LBMPSMLattice::get_B(int index){
  return B[index];
};

double LBMPSMLattice::get_rho(int index){
  return rho[index];
};

double LBMPSMLattice::get_Fhydx(int index){
  return Fhydx[index];
};

double LBMPSMLattice::get_Fhydy(int index){
  return Fhydy[index];
};

double LBMPSMLattice::get_Fhydz(int index){
  return Fhydz[index];
};

void LBMPSMLattice::set_Fhydx(int index, int pID, double Fhydx_)
{
  Fhydx[index] = Fhydx_;
  pData[index].hydrodynamicForce[3*pID+0] = Fhydx_;
}

void LBMPSMLattice::set_Fhydy(int index, int pID, double Fhydy_)
{
  Fhydy[index] = Fhydy_;
  pData[index].hydrodynamicForce[3*pID+1] = Fhydy_;
}

void LBMPSMLattice::set_Fhydz(int index, int pID, double Fhydz_)
{
  Fhydz[index] = Fhydz_;
  pData[index].hydrodynamicForce[3*pID+2] = Fhydz_;
}

void LBMPSMLattice::add_Fhydx(int index, int pID, double Fhydx_)
{
  Fhydx[index] += Fhydx_;
  pData[index].hydrodynamicForce[3*pID+0] += Fhydx_;
}

void LBMPSMLattice::add_Fhydy(int index, int pID, double Fhydy_)
{
  Fhydy[index] += Fhydy_;
  pData[index].hydrodynamicForce[3*pID+1] += Fhydy_;
}

void LBMPSMLattice::add_Fhydz(int index, int pID, double Fhydz_)
{
  Fhydz[index] += Fhydz_;
  pData[index].hydrodynamicForce[3*pID+2] += Fhydz_;
}


// Extend to two or morge particles
void LBMPSMLattice::setParticleOnLattice(int index, LAMMPS_NS::tagint pID, double uP[3], double eps)
{
  if(pData[index].particleID[0] == pID){
    pData[index].particleID[0] = pID;
    pData[index].solidFraction[0] = eps;
    pData[index].particleVelocity[0] = uP[0];
    pData[index].particleVelocity[1] = uP[1];
    pData[index].particleVelocity[2] = uP[2];
    pData[index].hydrodynamicForce[0] = 0.0;
    pData[index].hydrodynamicForce[1] = 0.0;
    pData[index].hydrodynamicForce[2] = 0.0;
  }
  if(pData[index].particleID[1] == pID){
    pData[index].particleID[1] = pID;
    pData[index].solidFraction[1] = eps;
    pData[index].particleVelocity[3] = uP[0];
    pData[index].particleVelocity[4] = uP[1];
    pData[index].particleVelocity[5] = uP[2];
    pData[index].hydrodynamicForce[3] = 0.0;
    pData[index].hydrodynamicForce[4] = 0.0;
    pData[index].hydrodynamicForce[5] = 0.0;
  }

  if(pID != pData[index].particleID[0] && pID != pData[index].particleID[1] && pData[index].particleID[0] == 0){
    pData[index].particleID[0] = pID;
    pData[index].solidFraction[0] = eps;
    pData[index].particleVelocity[0] = uP[0];
    pData[index].particleVelocity[1] = uP[1];
    pData[index].particleVelocity[2] = uP[2];
    pData[index].hydrodynamicForce[0] = 0.0;
    pData[index].hydrodynamicForce[1] = 0.0;
    pData[index].hydrodynamicForce[2] = 0.0;
  }
  if(pID != pData[index].particleID[0] && pID != pData[index].particleID[1] && pData[index].particleID[1] == 0){
    pData[index].particleID[1] = pID;
    pData[index].solidFraction[1] = eps;
    pData[index].particleVelocity[3] = uP[0];
    pData[index].particleVelocity[4] = uP[1];
    pData[index].particleVelocity[5] = uP[2];
    pData[index].hydrodynamicForce[3] = 0.0;
    pData[index].hydrodynamicForce[4] = 0.0;
    pData[index].hydrodynamicForce[5] = 0.0;
  }
}


void LBMPSMLattice::setToZero(int index, LAMMPS_NS::tagint pID)
{
  if(pData[index].particleID[0] == pID){
    pData[index].particleID[0] = 0;
    pData[index].solidFraction[0] = 0.0;
    pData[index].particleVelocity[0] = 0.0;
    pData[index].particleVelocity[1] = 0.0;
    pData[index].particleVelocity[2] = 0.0;
    pData[index].hydrodynamicForce[0] = 0.0;
    pData[index].hydrodynamicForce[1] = 0.0;
    pData[index].hydrodynamicForce[2] = 0.0;
  }
  if(pData[index].particleID[1] == pID){
    pData[index].particleID[1] = 0;
    pData[index].solidFraction[1] = 0.0;
    pData[index].particleVelocity[3] = 0.0;
    pData[index].particleVelocity[4] = 0.0;
    pData[index].particleVelocity[5] = 0.0;
    pData[index].hydrodynamicForce[3] = 0.0;
    pData[index].hydrodynamicForce[4] = 0.0;
    pData[index].hydrodynamicForce[5] = 0.0;
  }
}


double LBMPSMLattice::getSolidFractionOnLattice(int index, int pID)
{
  return pData[index].particleID[pID];
}

vector<double> LBMPSMLattice::getSolidVelocityOnLattice(int index, int pID)
{
  vector<double> returnVelVector;
  returnVelVector.resize(3);
  if(pData[index].particleID[0] == pID){
    returnVelVector[0] = pData[index].particleVelocity[0];
    returnVelVector[1] = pData[index].particleVelocity[1];
    returnVelVector[2] = pData[index].particleVelocity[2];
  }
  else if(pData[index].particleID[1] == pID){
    returnVelVector[0] = pData[index].particleVelocity[3];
    returnVelVector[1] = pData[index].particleVelocity[4];
    returnVelVector[2] = pData[index].particleVelocity[5];
  }
  else
  {
    returnVelVector[0] = 0.0;
    returnVelVector[1] = 0.0;
    returnVelVector[2] = 0.0;
  }
  return returnVelVector;
}


ParticleDataOnLattice LBMPSMLattice::getParticleDataOnLatticeNode(int index)
{
  return pData[index];
}


vector<int> LBMPSMLattice::get_procCoordinates()
{
  vector<int> returnProcCoordinates {procCoordinates[0], procCoordinates[1], procCoordinates[2] };
  return returnProcCoordinates;
}
