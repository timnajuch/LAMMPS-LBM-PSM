/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#include "LBM_PSM_lattice.h"


LBMPSMLattice::LBMPSMLattice(int nx_, int ny_, int nz_, int decomposition[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_){
  envelopeWidth = 1;
  dimension = dimension_;

  procCoordinates[0] = procCoordinates_[0];
  procCoordinates[1] = procCoordinates_[1];
  procCoordinates[2] = procCoordinates_[2];

  nxLocal.resize(decomposition[0]*decomposition[1]*decomposition[2]);
  nyLocal.resize(decomposition[0]*decomposition[1]*decomposition[2]);
  nzLocal.resize(decomposition[0]*decomposition[1]*decomposition[2]);

  cs = 1.0/sqrt(3.0);
  csPow2 = 1.0/3.0;
  currentStep = 0;
  nextStep = 1;
  
  if ( nx_ % decomposition[0] == 0){
    nx = nx_/decomposition[0] + envelopeWidth*2;
  }else{
    if (procCoordinates[0] != decomposition[0]-1){
      nx = nx_/decomposition[0] + envelopeWidth*2;
    }else{
      nx = ((int)(nx_/decomposition[0])) + nx_ - ((int)(nx_/decomposition[0]))*decomposition[0] + envelopeWidth*2;
    }
  }

  if ( ny_ % decomposition[1] == 0){
    ny = ny_/decomposition[1] + envelopeWidth*2;
  }else{
    if (procCoordinates[1] != decomposition[1]-1){
      ny = ny_/decomposition[1] + envelopeWidth*2;
    }else{
      ny = ((int)(ny_/decomposition[1])) + ny_ - ((int)(ny_/decomposition[1]))*decomposition[1] + envelopeWidth*2;
    }
  }

  nz = 1;
  q = 9;
  if (dimension == 3)
  {
    if ( nz_ % decomposition[2] == 0){
      nz = nz_/decomposition[2] + envelopeWidth*2;
    }else{
      if (procCoordinates[2] != decomposition[2]-1){
        nz = nz_/decomposition[2] + envelopeWidth*2;
      }else{
        nz = ((int)(nz_/decomposition[2])) + nz_ - ((int)(nz_/decomposition[2]))*decomposition[2] + envelopeWidth*2;
      }
    }

    q = 19;
  }

  MPI_Allgather(&nx, 1, MPI_INT, &(nxLocal[0]), 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&ny, 1, MPI_INT, &(nyLocal[0]), 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&nz, 1, MPI_INT, &(nzLocal[0]), 1, MPI_INT, MPI_COMM_WORLD);

  nxLocalGrid.resize(decomposition[0]);
  nyLocalGrid.resize(decomposition[1]);
  nzLocalGrid.resize(decomposition[2]);

  for(int iproc = 0; iproc < decomposition[0]; iproc++){
    int jproc = 0;
    int kproc = 0;
    int procIndex = iproc*decomposition[1]*decomposition[2] + jproc*decomposition[2] + kproc;
    nxLocalGrid[iproc] = nxLocal[procIndex];
  }

  for(int jproc = 0; jproc < decomposition[1]; jproc++){
    int iproc = 0;
    int kproc = 0;
    int procIndex = iproc*decomposition[1]*decomposition[2] + jproc*decomposition[2] + kproc;
    nyLocalGrid[jproc] = nyLocal[procIndex];
  }

  for(int kproc = 0; kproc < decomposition[2]; kproc++){
    int iproc = 0;
    int jproc = 0;
    int procIndex = iproc*decomposition[1]*decomposition[2] + jproc*decomposition[2] + kproc;
    nzLocalGrid[kproc] = nzLocal[procIndex];
  }

  f = vector<double>(nx*ny*nz*q*2,0.0); // Store populations at two different timesteps (even and odd to avoid using a temporary variable and copying)
  f0 = vector<double>(nx*ny*nz*q,0.0);

  origin_global.push_back(origin_[0]);
  origin_global.push_back(origin_[1]);
  origin_global.push_back(origin_[2]);

  B = vector<double>(nx*ny*nz,0.0);
  rho = vector<double>(nx*ny*nz,0.0);
  x = vector<double>(nx*ny*nz,0.0);
  y = vector<double>(nx*ny*nz,0.0);
  z = vector<double>(nx*ny*nz,0.0);

  u = vector<double>(nx*ny*nz*3,0.0);
  us = vector<double>(nx*ny*nz*3,0.0);

  pData.resize(nx*ny*nz);

  if(q == 9){ // D2Q9 lattice
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
    }

  if(q == 19){ // D3Q19 lattice
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

}


LBMPSMLattice::~LBMPSMLattice(){}


void LBMPSMLattice::initialise_domain(double dx_, double dy_, double dz_){
  dx = dx_;
  dy = dy_;
  dz = dz_;
  
  for(int i = 0; i < nx; ++i){
    for(int j = 0; j < ny; ++j){
      for(int k = 0; k < nz; ++k){

        int nxTotal = 0;
        for(int iproc = 0; iproc < procCoordinates[0]; iproc++){ nxTotal += get_nxLocal(iproc); }

        int nyTotal = 0;
        for(int jproc = 0; jproc < procCoordinates[1]; jproc++){ nyTotal += get_nyLocal(jproc); }

        int nzTotal = 0;
        for(int kproc = 0; kproc < procCoordinates[2]; kproc++){ nzTotal += get_nzLocal(kproc); }

        int index = index_1D(i, j, k);

        x[index] = origin_global[0] + i*dx - dx*envelopeWidth + (nxTotal - procCoordinates[0]*2*envelopeWidth)*dx;
        y[index] = origin_global[1] + j*dy - dy*envelopeWidth + (nyTotal - procCoordinates[1]*2*envelopeWidth)*dy;
        z[index] = origin_global[2] + k*dz - dz*envelopeWidth + (nzTotal - procCoordinates[2]*2*envelopeWidth)*dz;

        pData[index].particleID[0] = 0;
        pData[index].solidFraction[0] = 0.0;
        pData[index].particleVelocity[0] = 0.0;
        pData[index].particleVelocity[1] = 0.0;
        pData[index].particleVelocity[2] = 0.0;
        pData[index].hydrodynamicForce[0] = 0.0;
        pData[index].hydrodynamicForce[1] = 0.0;
        pData[index].hydrodynamicForce[2] = 0.0;
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
  }

}


int LBMPSMLattice::index_1D(int i, int j, int k){ return i*ny*nz + j*nz + k; }

int LBMPSMLattice::index_2D(int i, int j, int k, int direction){ return (i*ny*nz + j*nz + k)*3 + direction; }

int LBMPSMLattice::index_fi(int i, int j, int k, int iq, int step){ return nx*ny*nz*q*step + (i*ny*nz + j*nz + k)*q + iq; }

int LBMPSMLattice::get_currentStep(){ return currentStep; }

void LBMPSMLattice::set_currentStep(int currentStep_){ currentStep = currentStep_; nextStep = 1 - currentStep_; }

void LBMPSMLattice::set_f(int i_, int j_, int k_, int iq_, int step_, double value_){ f[index_fi(i_, j_, k_, iq_, step_)] = value_; }

void LBMPSMLattice::set_f(int ind_iq_, double value_){ f[ind_iq_] = value_; }

double LBMPSMLattice::get_f(int i_, int j_, int k_, int iq_, int step_){ return f[index_fi(i_, j_, k_, iq_, step_)]; }

double LBMPSMLattice::get_f(int ind_iq_){ return f[ind_iq_]; }

 void LBMPSMLattice::set_f0(int i_, int j_, int k_, int iq_, double value_){ f0[index_fi(i_, j_, k_, iq_, 0)] = value_; }

double LBMPSMLattice::get_f0(int i_, int j_, int k_, int iq_){ return f0[index_fi(i_, j_, k_, iq_, 0)]; }

double LBMPSMLattice::get_f0(int ind_iq_){ return f0[ind_iq_]; }

vector<double> LBMPSMLattice::get_B(){ return B; }

vector<double> LBMPSMLattice::get_rho(){ return rho; }

vector<double> LBMPSMLattice::get_x(){ return x; }

vector<double> LBMPSMLattice::get_y(){ return y; }

vector<double> LBMPSMLattice::get_z(){ return z; }

vector<double> LBMPSMLattice::get_u(){ return u; }

vector<double>& LBMPSMLattice::get_B_reference(){ return B; }

vector<double>& LBMPSMLattice::get_rho_reference(){ return rho; }

vector<double>& LBMPSMLattice::get_x_reference(){ return x; }

vector<double>& LBMPSMLattice::get_y_reference(){ return y; }

vector<double>& LBMPSMLattice::get_z_reference(){ return z; }

vector<double>& LBMPSMLattice::get_u_reference(){ return u; }

int LBMPSMLattice::get_nx(){ return nx; }

int LBMPSMLattice::get_ny(){ return ny; }

int LBMPSMLattice::get_nz(){ return nz; }

int LBMPSMLattice::get_envelopeWidth(){ return envelopeWidth; }

int LBMPSMLattice::get_q(){ return q; }

int LBMPSMLattice::get_nxLocal(int iProcIndex){ return nxLocalGrid[iProcIndex]; }

int LBMPSMLattice::get_nyLocal(int jProcIndex){ return nyLocalGrid[jProcIndex]; }

int LBMPSMLattice::get_nzLocal(int kProcIndex){ return nzLocalGrid[kProcIndex]; }

void LBMPSMLattice::set_B(int index, double B_){ B[index] = B_; }

double LBMPSMLattice::get_B(int index){ return B[index]; }

double LBMPSMLattice::get_rho(int index){ return rho[index]; }

double LBMPSMLattice::get_u(int index){ return u[index]; }

double LBMPSMLattice::get_u_at_node(int index_node_1D, int direction){ return u[index_node_1D*3+direction]; }

vector<double>& LBMPSMLattice::getVector_f(){ return f; }

void LBMPSMLattice::setVector_f(vector<double>& fcopy){ f = fcopy; };


ParticleDataOnLattice LBMPSMLattice::getParticleDataOnLatticeNode(int index){ return pData[index]; }

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
  else if(pData[index].particleID[1] == pID){
    pData[index].particleID[1] = pID;
    pData[index].solidFraction[1] = eps;
    pData[index].particleVelocity[3] = uP[0];
    pData[index].particleVelocity[4] = uP[1];
    pData[index].particleVelocity[5] = uP[2];
    pData[index].hydrodynamicForce[3] = 0.0;
    pData[index].hydrodynamicForce[4] = 0.0;
    pData[index].hydrodynamicForce[5] = 0.0;
  }
  else if(pID != pData[index].particleID[0] && pID != pData[index].particleID[1] && pData[index].particleID[0] == 0){
    pData[index].particleID[0] = pID;
    pData[index].solidFraction[0] = eps;
    pData[index].particleVelocity[0] = uP[0];
    pData[index].particleVelocity[1] = uP[1];
    pData[index].particleVelocity[2] = uP[2];
    pData[index].hydrodynamicForce[0] = 0.0;
    pData[index].hydrodynamicForce[1] = 0.0;
    pData[index].hydrodynamicForce[2] = 0.0;
  }
  else if(pID != pData[index].particleID[0] && pID != pData[index].particleID[1] && pData[index].particleID[1] == 0){
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
  else if(pData[index].particleID[1] == pID){
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


double LBMPSMLattice::getSolidFractionOnLattice(int index, int pID){ return pData[index].solidFraction[pID]; }

vector<double> LBMPSMLattice::getSolidVelocityOnLattice(int index){ return pData[index].particleVelocity; }

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


void LBMPSMLattice::add_Fhyd(int index, LAMMPS_NS::tagint pID, double Fhyd, int dir)
{
  if(pData[index].particleID[0] == pID){
    pData[index].hydrodynamicForce[3*0+dir] += Fhyd;
  }
  else if(pData[index].particleID[1] == pID){
    pData[index].hydrodynamicForce[3*1+dir] += Fhyd;
  }
  else if(pID != pData[index].particleID[0] && pID != pData[index].particleID[1] && pData[index].particleID[0] == 0){
    pData[index].hydrodynamicForce[3*0+dir] += Fhyd;
  }
  else if(pID != pData[index].particleID[0] && pID != pData[index].particleID[1] && pData[index].particleID[1] == 0){
    pData[index].hydrodynamicForce[3*1+dir] += Fhyd;
  }
}


vector<int> LBMPSMLattice::get_procCoordinates()
{
  vector<int> returnProcCoordinates {procCoordinates[0], procCoordinates[1], procCoordinates[2] };
  return returnProcCoordinates;
}