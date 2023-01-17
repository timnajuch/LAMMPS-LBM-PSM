/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#include "fix_LBM_PSM_writeVTK.h"

WriteVTK::WriteVTK(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {

  if (narg < 5) error->all(FLERR,"Illegal fix lbm-psm-vtk command");

  for(int ifix=0; ifix<modify->nfix; ifix++)
    if(strcmp(modify->fix[ifix]->style,"lbm-psm")==0)
      fixLBMPSM = (fix_LBM_PSM *)modify->fix[ifix];

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lbm-psm-vtk command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix lbm-psm-vtk command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix lbm-psm-vtk command");
  }
};


WriteVTK::~WriteVTK() {};


int WriteVTK::setmask()
{
  int mask =0;
  mask |= PRE_FORCE;
  return mask;
}


void WriteVTK::init()
{
  decomposition[0] = comm->procgrid[0];
  decomposition[1] = comm->procgrid[1];
  decomposition[2] = comm->procgrid[2];

  nx = fixLBMPSM->dynamics->get_nx();
  ny = fixLBMPSM->dynamics->get_ny();
  nz = fixLBMPSM->dynamics->get_nz();

  vector<double> x_vtk = fixLBMPSM->dynamics->get_x();
  vector<double> y_vtk = fixLBMPSM->dynamics->get_y();
  vector<double> z_vtk = fixLBMPSM->dynamics->get_z();
  vector<double> rho_vtk = fixLBMPSM->dynamics->get_rho();
  vector<double> u_vtk = fixLBMPSM->dynamics->get_u();

  vector<double> B_vtk;
  B_vtk.resize(nx*ny*nz);
  for (int i = 0; i < nx*ny*nz; ++i){
    B_vtk[i] = fixLBMPSM->dynamics->getParticleDataOnLatticeNode(i).solidFraction[0] + fixLBMPSM->dynamics->getParticleDataOnLatticeNode(i).solidFraction[1];
  }
  double u_infty = 1.0;
  double ly = 1.0;

  write_vtk("init.vtk", x_vtk, 1.0/ly, y_vtk, 1.0/ly, z_vtk, 1.0/ly, B_vtk, 1.0, rho_vtk, 1000.0, u_vtk, 1.0/u_infty);
}


void WriteVTK::pre_force(int)
{
  if (update->ntimestep % nevery) return;
  vector<double> x_vtk = fixLBMPSM->dynamics->get_x();
  vector<double> y_vtk = fixLBMPSM->dynamics->get_y();
  vector<double> z_vtk = fixLBMPSM->dynamics->get_z();
  vector<double> rho_vtk = fixLBMPSM->dynamics->get_rho();
  vector<double> u_vtk = fixLBMPSM->dynamics->get_u();
  vector<double> B_vtk;
  B_vtk.resize(nx*ny*nz);
  for (int i = 0; i < nx*ny*nz; ++i){
    B_vtk[i] = fixLBMPSM->dynamics->getParticleDataOnLatticeNode(i).solidFraction[0] + fixLBMPSM->dynamics->getParticleDataOnLatticeNode(i).solidFraction[1];
  }
  double u_infty = fixLBMPSM->unitConversion->get_u_lb();
  double Uc = fixLBMPSM->unitConversion->get_Uc();
  double ly = 1.0;

  ostringstream timeStringTmp;
  timeStringTmp << "flowField_" << setw(10) << setfill('0') << update->ntimestep << ".vtk";
  string timeString(timeStringTmp.str());

  write_vtk(timeString, x_vtk, 1.0/ly, y_vtk, 1.0/ly, z_vtk, 1.0/ly, B_vtk, 1.0, rho_vtk, 1000.0, u_vtk, Uc/u_infty);
}


void WriteVTK::write_vtk(string name_, vector<double> &x_, double x0_, vector<double> &y_, double y0_, vector<double> &z_, double z0_, vector<double> &B_, double B0_, vector<double> &rho_, double rho0_, vector<double> &u_, double u0_){
  int envelopeWidth = 1;
  int nzLoopStart = 0;
  int nzLoopEnd = 1;
  if(domain->dimension == 3){
    nzLoopStart = envelopeWidth;
    nzLoopEnd = nz-envelopeWidth;
  }

  vector<double> x0;
  vector<double> y0;
  vector<double> z0;
  vector<double> B0;
  vector<double> rho0;
  vector<double> u0;
  if(comm->me == 0){
    x0.resize(nx*decomposition[0]*ny*decomposition[1]*nz*decomposition[2]);
    y0.resize(nx*decomposition[0]*ny*decomposition[1]*nz*decomposition[2]);
    z0.resize(nx*decomposition[0]*ny*decomposition[1]*nz*decomposition[2]);
    B0.resize(nx*decomposition[0]*ny*decomposition[1]*nz*decomposition[2]);
    rho0.resize(nx*decomposition[0]*ny*decomposition[1]*nz*decomposition[2]);
    u0.resize(nx*decomposition[0]*ny*decomposition[1]*nz*decomposition[2]*3);
  }

  MPI_Gather(&x_[0], nx*ny*nz, MPI_DOUBLE, &x0[0], nx*ny*nz, MPI_DOUBLE, 0, world);
  MPI_Gather(&y_[0], nx*ny*nz, MPI_DOUBLE, &y0[0], nx*ny*nz, MPI_DOUBLE, 0, world);
  MPI_Gather(&z_[0], nx*ny*nz, MPI_DOUBLE, &z0[0], nx*ny*nz, MPI_DOUBLE, 0, world);
  MPI_Gather(&B_[0], nx*ny*nz, MPI_DOUBLE, &B0[0], nx*ny*nz, MPI_DOUBLE, 0, world);
  MPI_Gather(&rho_[0], nx*ny*nz, MPI_DOUBLE, &rho0[0], nx*ny*nz, MPI_DOUBLE, 0, world);
  MPI_Gather(&u_[0], nx*ny*nz*3, MPI_DOUBLE, &u0[0], nx*ny*nz*3, MPI_DOUBLE, 0, world);

  if(comm->me == 0){

  vector<double> point_id(nx*decomposition[0]*ny*decomposition[1]*nz*decomposition[2], 0.0);

  vector<double> x_vtk = x0;
  scale_vector(x_vtk, x0_);

  vector<double> y_vtk = y0;
  scale_vector(y_vtk, y0_);

  vector<double> z_vtk = z0;
  scale_vector(z_vtk, z0_);

  vector<double> B_vtk = B0;

  vector<double> rho_vtk = rho0;
  scale_vector(rho_vtk, rho0_);

  vector<double> u_vtk = u0;
  scale_vector(u_vtk, u0_);

  ofstream ovel;
  ovel.open(name_);
  ovel << "# vtk DataFile Version 3.1 \n";
  ovel << "Fluid flow field from LBM simulation\n";
  ovel << "ASCII\n";
  ovel << "DATASET RECTILINEAR_GRID\n";
  if(domain->dimension == 2){
    ovel << "DIMENSIONS " << (nx-2*envelopeWidth)*decomposition[0] << " " << (ny-2*envelopeWidth)*decomposition[1] << " " << 1 << "\n";
  }else{
    ovel << "DIMENSIONS " << (nx-2*envelopeWidth)*decomposition[0] << " " << (ny-2*envelopeWidth)*decomposition[1] << " " << (nz-2*envelopeWidth)*decomposition[2] << "\n";
  }

  ovel << "X_COORDINATES " << (nx-2*envelopeWidth)*decomposition[0] << " float" << endl;
  for(int iproc = 0; iproc < decomposition[0]; iproc++){
    for(int i = envelopeWidth; i < nx-envelopeWidth; i++){
      int j = 0; int jproc = 0;
      int k = 0; int kproc = 0;
      int index = i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*nz*decomposition[2] + kproc*nx*ny*nz;
      ovel << x_vtk[index] << "\n";
    }
  }

  ovel << "Y_COORDINATES " << (ny-2*envelopeWidth)*decomposition[1] << " float" << endl;
  for(int jproc = 0; jproc < decomposition[1]; jproc++){
    for(int j = envelopeWidth; j < ny-envelopeWidth; j++){
      int i = 0; int iproc = 0;
      int k = 0; int kproc = 0;
      int index = i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*nz*decomposition[2] + kproc*nx*ny*nz;
      ovel << y_vtk[index] << "\n";
    }
  }

  if(domain->dimension == 2){
    ovel << "Z_COORDINATES " << 1 << " float" << endl;
    ovel << "0.0\n";
  }else{
    ovel << "Z_COORDINATES " << (nz-2*envelopeWidth)*decomposition[2] << " float" << endl;
    for(int kproc = 0; kproc < decomposition[2]; kproc++){
      for(int k = envelopeWidth; k < nz-envelopeWidth; k++){
        int i = 0; int iproc = 0;
        int j = 0; int jproc = 0;
        int index = i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*nz*decomposition[2] + kproc*nx*ny*nz;
        ovel << z_vtk[index] << "\n";
      }
    }
  }

  // POINT_DATA
  if(domain->dimension == 2){
    ovel << "\nPOINT_DATA " << (nx-envelopeWidth*2)*decomposition[0]*(ny-envelopeWidth*2)*decomposition[1] << "\n";
  }else{
    ovel << "\nPOINT_DATA " << (nx-envelopeWidth*2)*decomposition[0]*(ny-envelopeWidth*2)*decomposition[1]*(nz-envelopeWidth*2)*decomposition[2] << "\n";
  }
  ovel << "SCALARS SolidFraction FLOAT" << "\n";
  ovel << "LOOKUP_TABLE default" << "\n";
  for(int kproc = 0; kproc<decomposition[2]; kproc++){
    for(int k=nzLoopStart; k<nzLoopEnd; k++){
      for(int jproc = 0; jproc<decomposition[1]; jproc++){
        for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
          for(int iproc = 0; iproc<decomposition[0]; iproc++){
            for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
              int index = i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*nz*decomposition[2] + kproc*nx*ny*nz;
              ovel << B_vtk[index] <<"\n";
            }
          }
        }
      }
    }
  }

  ovel << "\nSCALARS Density FLOAT" << "\n";
  ovel << "LOOKUP_TABLE default" << "\n";
  for(int kproc = 0; kproc<decomposition[2]; kproc++){
    for(int k=nzLoopStart; k<nzLoopEnd; k++){
      for(int jproc = 0; jproc<decomposition[1]; jproc++){
        for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
          for(int iproc = 0; iproc<decomposition[0]; iproc++){
            for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
              int index = i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*nz*decomposition[2] + kproc*nx*ny*nz;
              ovel << rho_vtk[index] << "\n";
            }
          }
        }
      }
    }
  }

  ovel << "\nSCALARS Velocity-x FLOAT" << "\n";
  ovel << "LOOKUP_TABLE default" << "\n";
  for(int kproc = 0; kproc<decomposition[2]; kproc++){
    for(int k=nzLoopStart; k<nzLoopEnd; k++){
      for(int jproc = 0; jproc<decomposition[1]; jproc++){
        for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
          for(int iproc = 0; iproc<decomposition[0]; iproc++){
            for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
              int index = (i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*nz*decomposition[2] + kproc*nx*ny*nz)*3;
              ovel << u_vtk[index] <<"\n";
            }
          }
        }
      }
    }
  }

  ovel << "\nSCALARS Velocity-y FLOAT" << "\n";
  ovel << "LOOKUP_TABLE default" << "\n";
  for(int kproc = 0; kproc<decomposition[2]; kproc++){
    for(int k=nzLoopStart; k<nzLoopEnd; k++){
      for(int jproc = 0; jproc<decomposition[1]; jproc++){
        for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
          for(int iproc = 0; iproc<decomposition[0]; iproc++){
            for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
              int index = (i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*nz*decomposition[2] + kproc*nx*ny*nz)*3+1;
              ovel << u_vtk[index] <<"\n";
            }
          }
        }
      }
    }
  }

  if(domain->dimension == 3){
    ovel << "\nSCALARS Velocity-z FLOAT" << "\n";
    ovel << "LOOKUP_TABLE default" << "\n";
    for(int kproc = 0; kproc<decomposition[2]; kproc++){
      for(int k=nzLoopStart; k<nzLoopEnd; k++){
        for(int jproc = 0; jproc<decomposition[1]; jproc++){
          for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
            for(int iproc = 0; iproc<decomposition[0]; iproc++){
              for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
                int index = (i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*nz*decomposition[2] + kproc*nx*ny*nz)*3+2;
                ovel << u_vtk[index] <<"\n";
              }
            }
          }
        }
      }
    }
  }

  ovel.close();
  }
};


// TODO extend to 3D if needed
void WriteVTK::write_profile(string name_, int ix_, vector<double> &y_, double y0_, vector<double> &B_, double B0_, vector<double> &rho_, double rho0_, vector<double> &u_, double u0_){
  ofstream ovel;
  ovel.open(name_);
  ovel << "y ux_(Dimensionless) uy_(Dimensionless) rho B" << "\n";
  for(int j=0; j<ny*decomposition[1]; j++){
    ovel << y_[ix_*ny*decomposition[1]+j]*y0_ << " " << u_[ix_*ny*decomposition[1]*2 + j*2]*u0_ << " " << u_[ix_*ny*decomposition[1]*2 + j*2+1]*u0_ << " " << rho_[ix_*ny*decomposition[1]+j]*rho0_ << 
      " " << B_[ix_*ny*decomposition[1]+j]*B0_ << "\n";
  }
  ovel.close();
};


void WriteVTK::scale_vector(vector<double> &vec_, double scaling_){
  std::transform(vec_.begin(), vec_.end(), vec_.begin(),
               std::bind(std::multiplies<double>(), std::placeholders::_1, scaling_));
};
