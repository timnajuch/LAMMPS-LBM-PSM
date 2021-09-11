/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#include "lattice2D.h"


Lattice2D::Lattice2D(int nx_, int ny_, int q_, int decomposition[3], int procCoordinates_[3]){
  envelopeWidth = 1;

  procCoordinates[0] = procCoordinates_[0];
  procCoordinates[1] = procCoordinates_[1];
  procCoordinates[2] = procCoordinates_[2];

  c = 1.0;
  cs = 1.0/sqrt(3.0);
  csPow2 = 1.0/3.0;
  nx = nx_/decomposition[0] + envelopeWidth*2;
  ny = ny_/decomposition[1] + envelopeWidth*2;
  q = q_;

  f = vector<double>(nx*ny*q,0.0);
  f0 = vector<double>(nx*ny*q,0.0);
  fcoll = vector<double>(nx*ny*q,0.0);

  B = vector<double>(nx*ny,0.0);
  rho = vector<double>(nx*ny,0.0);
  x = vector<double>(nx*ny,0.0);
  y = vector<double>(nx*ny,0.0);
  
  u = vector<double>(nx*ny*2,0.0);
  us = vector<double>(nx*ny*2,0.0);

  Fhydx = vector<double>(nx*ny,0.0);
  Fhydy = vector<double>(nx*ny,0.0);

  pData.resize(nx*ny);

  if(q == 9){
    e = { 0.0 ,  0.0, 
          1.0 ,  0.0,
          0.0 ,  1.0,
          -1.0,  0.0,
          0.0 , -1.0,
          1.0 ,  1.0,
          -1.0,  1.0,
          -1.0, -1.0,
          1.0 ,  -1.0
        };

    w = { 4.0/9.0,
          1.0/9.0,
          1.0/9.0,
          1.0/9.0,
          1.0/9.0,
          1.0/36.0,
          1.0/36.0,
          1.0/36.0,
          1.0/36.0,
        };
    
    g = { 0.0,
          1.0/3.0,
          1.0/3.0,
          1.0/3.0,
          1.0/3.0,
          1.0/12.0,
          1.0/12.0,
          1.0/12.0,
          1.0/12.0,
        };
    }

};

Lattice2D::~Lattice2D(){};

void Lattice2D::initialise_channel_geometry(double wallHeightHalf_, double eps_, double nychannel_, double dx_, double dy_){
  dx = dx_;
  dy = dy_;
  
  for(int i = 0; i < nx; ++i){
    for(int j = 0; j < ny; ++j){

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

};

void Lattice2D::initialise_domain(double dx_, double dy_){
  dx = dx_;
  dy = dy_;
  
  for(int i = 0; i < nx; ++i){
    for(int j = 0; j < ny; ++j){

      int ind_2D = i*ny + j;

      x[ind_2D] = i*dx - dx*envelopeWidth + procCoordinates[0]*(nx-2*envelopeWidth)*dx;
      y[ind_2D] = j*dy - dy*envelopeWidth + procCoordinates[1]*(ny-2*envelopeWidth)*dy;


      Fhydx[ind_2D] = 0.0;
      Fhydy[ind_2D] = 0.0;
 /*
      pData[ind_2D].particleID[0] = 0;
      pData[ind_2D].particleID[1] = 0;
      pData[ind_2D].solidFraction[0] = 0.0;
      pData[ind_2D].solidFraction[1] = 0.0;
      pData[ind_2D].particleVelocity[0] = 0.0;
      pData[ind_2D].particleVelocity[1] = 0.0;
      pData[ind_2D].particleVelocity[1] = 0.0;
      pData[ind_2D].particleVelocity[3] = 0.0;
      pData[ind_2D].particleVelocity[4] = 0.0;
      pData[ind_2D].particleVelocity[5] = 0.0;
      pData[ind_2D].hydrodynamicForce[0] = 0.0;
      pData[ind_2D].hydrodynamicForce[1] = 0.0;
      pData[ind_2D].hydrodynamicForce[2] = 0.0;
      pData[ind_2D].hydrodynamicForce[3] = 0.0;
      pData[ind_2D].hydrodynamicForce[4] = 0.0;
      pData[ind_2D].hydrodynamicForce[5] = 0.0;
  */
    }
  }

};

void Lattice2D::set_f(int i_, int j_, int iq_, double value_){
  f[i_*ny*q + j_*q + iq_] = value_;
};

void Lattice2D::set_f(int ind_iq_, double value_){
  f[ind_iq_] = value_;
};

void Lattice2D::cp_fcoll_f(){
  f = fcoll;
};

double Lattice2D::get_f(int i_, int j_, int iq_){
  return f[i_*ny*q + j_*q + iq_];
};

double Lattice2D::get_f(int ind_iq_){
  return f[ind_iq_];
};

 void Lattice2D::set_f0(int i_, int j_, int iq_, double value_){
  f0[i_*ny*q + j_*q + iq_] = value_;
};

double Lattice2D::get_f0(int i_, int j_, int iq_){
  return f0[i_*ny*q + j_*q + iq_];
};

double Lattice2D::get_f0(int ind_iq_){
  return f0[ind_iq_];
};

 void Lattice2D::set_fcoll(int i_, int j_, int iq_, double value_){
  fcoll[i_*ny*q + j_*q + iq_] = value_;
};

double Lattice2D::get_fcoll(int i_, int j_, int iq_){
  return fcoll[i_*ny*q + j_*q + iq_];
};

double Lattice2D::get_fcoll(int ind_iq_){
  return fcoll[ind_iq_];
};

 vector<double> Lattice2D::get_B(){
  return B;
};

 vector<double> Lattice2D::get_rho(){
  return rho;
};

 vector<double> Lattice2D::get_x(){
  return x;
};

 vector<double> Lattice2D::get_y(){
  return y;
};

 vector<double> Lattice2D::get_u(){
  return u;
};

int Lattice2D::get_nx(){
  return nx;
};

int Lattice2D::get_ny(){
  return ny;
};

void Lattice2D::set_B(int index, double B_){
  B[index] = B_;
};

double Lattice2D::get_B(int index){
  return B[index];
};

double Lattice2D::get_rho(int index){
  return rho[index];
};

double Lattice2D::get_Fhydx(int index){
  return Fhydx[index];
};

double Lattice2D::get_Fhydy(int index){
  return Fhydy[index];
};

void Lattice2D::set_Fhydx(int index, double Fhydx_)
{
  Fhydx[index] = Fhydx_;
}

void Lattice2D::set_Fhydy(int index, double Fhydy_)
{
  Fhydy[index] = Fhydy_;
}

void Lattice2D::add_Fhydx(int index, double Fhydx_)
{
  Fhydx[index] += Fhydx_;
}

void Lattice2D::add_Fhydy(int index, double Fhydy_)
{
  Fhydy[index] += Fhydy_;
}


// Extend to two or morge particles
void Lattice2D::setParticleOnLattice(int index, int pID, double uP[2], double eps)
{
  pData[index].particleID[0] = pID;
  pData[index].solidFraction[0] = eps;
  pData[index].particleVelocity[0] = uP[0];
  pData[index].particleVelocity[1] = uP[1];
  pData[index].hydrodynamicForce[0] = 0.0;
  pData[index].hydrodynamicForce[1] = 0.0;
  B[index] = pData[index].solidFraction[0];
}

// Todo modify so that it gets the solid fraction of a specific particle ID covering a node
double Lattice2D::getSolidFractionOnLattice(int index, int pID)
{
  return pData[index].solidFraction[0];
}

vector<double> Lattice2D::getSolidVelocityOnLattice(int index, int pID)
{
  return pData[index].particleVelocity;
}
