/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2019
------------------------------------------------------*/

#include "fix_writeVTK.h"

WriteVTK::WriteVTK(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {
  for(int ifix=0; ifix<modify->nfix; ifix++)
    if(strcmp(modify->fix[ifix]->style,"lbm-psm")==0)
      fixPSMLBM = (fix_PSM_LBM *)modify->fix[ifix];
};

WriteVTK::~WriteVTK() {};

int WriteVTK::setmask()
{
  int mask =0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= END_OF_STEP;
  return mask;
}


void WriteVTK::init()
{
  decomposition[0] = comm->procgrid[0];
  decomposition[1] = comm->procgrid[1];
  decomposition[2] = comm->procgrid[2];

  nx = fixPSMLBM->dynamics->get_nx();
  ny = fixPSMLBM->dynamics->get_ny();
  ny = fixPSMLBM->dynamics->get_nz();

  vector<double> x_vtk = fixPSMLBM->dynamics->get_x();
  vector<double> y_vtk = fixPSMLBM->dynamics->get_y();
  vector<double> z_vtk = fixPSMLBM->dynamics->get_z();
  vector<double> rho_vtk = fixPSMLBM->dynamics->get_rho();
  vector<double> u_vtk = fixPSMLBM->dynamics->get_u();
  vector<double> B_vtk = fixPSMLBM->dynamics->get_B();
  //vector<double> x_vtk = fixPSMLBM->get_x();
  //vector<double> y_vtk = fixPSMLBM->get_y();
  //vector<double> rho_vtk = fixPSMLBM->get_rho();
  //vector<double> u_vtk = fixPSMLBM->get_u();
  //vector<double> B_vtk = fixPSMLBM->get_B();
  double u_infty = 1.0;
  double ly = 1.0;
  write_vtk("test.vtk", x_vtk, 1.0/ly, y_vtk, 1.0/ly, z_vtk, 1.0/ly, B_vtk, 1.0, rho_vtk, 1000.0, u_vtk, 1.0/u_infty);
}


void WriteVTK::pre_force(int)
{
  vector<double> x_vtk = fixPSMLBM->dynamics->get_x();
  vector<double> y_vtk = fixPSMLBM->dynamics->get_y();
  vector<double> z_vtk = fixPSMLBM->dynamics->get_z();
  vector<double> rho_vtk = fixPSMLBM->dynamics->get_rho();
  vector<double> u_vtk = fixPSMLBM->dynamics->get_u();
  vector<double> B_vtk = fixPSMLBM->dynamics->get_B();
  double u_infty = fixPSMLBM->unitConversion->get_u_lb();
  double ly = 1.0;
  write_vtk("run.vtk", x_vtk, 1.0/ly, y_vtk, 1.0/ly, z_vtk, 1.0/ly, B_vtk, 1.0, rho_vtk, 1000.0, u_vtk, 1.0/u_infty);
}


//void WriteVTK::write_vtk(string name_, vector<double> &x_, vector<double> &y_, vector<double> &B_, vector<double> &rho_, vector<double> &u_, MPICommunication &mpicomm){
void WriteVTK::write_vtk(string name_, vector<double> &x_, vector<double> &y_, vector<double> &z_, vector<double> &B_, vector<double> &rho_, vector<double> &u_){
// TODO only first POINT block extended to 3D. check if this method is used/useful anyway
vector<double> x0;
vector<double> y0;
vector<double> z0;
vector<double> B0;
vector<double> rho0;
vector<double> u0;

vector<double> point_id(nx*decomposition[0]*ny*decomposition[1]*nz*decomposition[2], 0.0);

ofstream ovel;
ovel.open(name_);
ovel << "# vtk DataFile Version 2.0 \n";
ovel << "Fluid flow field from LBM simulation\n";
ovel << "ASCII\n";
ovel << "DATASET UNSTRUCTURED_GRID\n\n";
// POINTS
int pnt_id = -1;
ovel << "POINTS " << nx*decomposition[0]*ny*decomposition[1]*nz*decomposition[2] << " FLOAT\n";
for(int i=0; i<nx*decomposition[0]; i++){
  for(int j=0; j<ny*decomposition[1]; j++){
    for(int k=0; k<nz*decomposition[2]; k++){
      ovel << x_[i*ny*decomposition[1]*nz*decomposition[2]+j*nz*decomposition[2]+k] << " " << y_[i*ny*decomposition[1]*nz*decomposition[2]+j*nz*decomposition[2]+k] << " " << z_[i*ny*decomposition[1]*nz*decomposition[2]+j*nz*decomposition[2]+k]<< "\n";
      pnt_id += 1;
      point_id[i*ny*decomposition[1]*nz*decomposition[2]+j*nz*decomposition[2]+k] = pnt_id;
    }
  }
}
// CELLS
// Calculate no of cells
int no_cells = (nx*decomposition[0]-1)*(ny*decomposition[1]-1);
int no_cell_blocks = no_cells*5;
ovel << "\nCELLS " << no_cells << " " << no_cell_blocks;
for(int j=0; j<(ny*decomposition[1]-1); j++){
  for(int i=0; i<(nx*decomposition[0]-1); i++){
    ovel << "\n4 " << point_id[i*ny*decomposition[1]+j] << " " << point_id[(i+1)*ny*decomposition[1]+j] << " " << point_id[(i+1)*ny*decomposition[1]+(j+1)] << " " << point_id[i*ny*decomposition[1]+(j+1)];
  }
}
// CELL_TYPES
ovel << "\n\nCELL_TYPES " << no_cells << "\n";
for(int i = 0; i<no_cells; i++)
  ovel << "9 ";
// POINT_DATA
// Density
ovel << "\n\nPOINT_DATA " << nx*decomposition[0]*ny*decomposition[1] << "\n";
ovel << "SCALARS Density FLOAT" << "\n";
ovel << "LOOKUP_TABLE default" << "\n";
for(int i=0; i<nx*decomposition[0]; i++){
  for(int j=0; j<ny*decomposition[1]; j++){
    ovel << rho_[i*ny*decomposition[1]+j] << "\n";
  }
}

ovel << "\nSCALARS Velocity-x FLOAT" << "\n";
ovel << "LOOKUP_TABLE default" << "\n";
for(int i=0; i<nx*decomposition[0]; i++){
  for(int j=0; j<ny*decomposition[1]; j++){
    ovel << u_[i*ny*decomposition[1]*2 + j*2] <<"\n";
  }
}

ovel << "\nSCALARS Velocity-y FLOAT" << "\n";
ovel << "LOOKUP_TABLE default" << "\n";
for(int i=0; i<nx*decomposition[0]; i++){
  for(int j=0; j<ny*decomposition[1]; j++){
    ovel << u_[i*ny*decomposition[1]*2 + j*2+1] <<"\n";
  }
}

ovel << "\nSCALARS B-Field FLOAT" << "\n";
ovel << "LOOKUP_TABLE default" << "\n";
for(int i=0; i<nx*decomposition[0]; i++){
    for(int j=0; j<ny*decomposition[1]; j++){
        ovel << B_[i*ny*decomposition[1]+j] <<"\n";
    }
}

ovel.close();
};

//void WriteVTK::write_vtk(string name_, vector<double> &x_, double x0_, vector<double> &y_, double y0_, vector<double> &B_, double B0_, vector<double> &rho_, double rho0_, vector<double> &u_, double u0_, MPICommunication &mpicomm){
void WriteVTK::write_vtk(string name_, vector<double> &x_, double x0_, vector<double> &y_, double y0_, vector<double> &z_, double z0_, vector<double> &B_, double B0_, vector<double> &rho_, double rho0_, vector<double> &u_, double u0_){

  int envelopeWidth = 1;

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
//  MPI_Gather(&x_[0], nx*ny, MPI_DOUBLE, &x0[0], nx*ny, MPI_DOUBLE, 0, mpicomm.cartComm3D);
//  MPI_Gather(&y_[0], nx*ny, MPI_DOUBLE, &y0[0], nx*ny, MPI_DOUBLE, 0, mpicomm.cartComm3D);
//  MPI_Gather(&B_[0], nx*ny, MPI_DOUBLE, &B0[0], nx*ny, MPI_DOUBLE, 0, mpicomm.cartComm3D);
//  MPI_Gather(&rho_[0], nx*ny, MPI_DOUBLE, &rho0[0], nx*ny, MPI_DOUBLE, 0, mpicomm.cartComm3D);
//  MPI_Gather(&u_[0], nx*ny*2, MPI_DOUBLE, &u0[0], nx*ny*2, MPI_DOUBLE, 0, mpicomm.cartComm3D);

  MPI_Gather(&x_[0], nx*ny*nz, MPI_DOUBLE, &x0[0], nx*ny*nz, MPI_DOUBLE, 0, world);
  MPI_Gather(&y_[0], nx*ny*nz, MPI_DOUBLE, &y0[0], nx*ny*nz, MPI_DOUBLE, 0, world);
  MPI_Gather(&B_[0], nx*ny*nz, MPI_DOUBLE, &B0[0], nx*ny*nz, MPI_DOUBLE, 0, world);
  MPI_Gather(&rho_[0], nx*ny*nz, MPI_DOUBLE, &rho0[0], nx*ny*nz, MPI_DOUBLE, 0, world);
  MPI_Gather(&u_[0], nx*ny*nz*3, MPI_DOUBLE, &u0[0], nx*ny*nz*3, MPI_DOUBLE, 0, world);

  if(comm->me == 0){

  //vector<double> point_id((nx-envelopeWidth*2)*decomposition[0]*(ny-envelopeWidth*2)*decomposition[1], 0.0);
  vector<double> point_id(nx*decomposition[0]*ny*decomposition[1]*nz*decomposition[2], 0.0);

  vector<double> x_vtk = x0;
  //x_vtk = scale_vector(x_vtk, x0_);
  scale_vector(x_vtk, x0_);

  vector<double> y_vtk = y0;
  //y_vtk = scale_vector(y_vtk, y0_);
  scale_vector(y_vtk, y0_);

  vector<double> z_vtk = z0;
  //z_vtk = scale_vector(z_vtk, z0_);
  scale_vector(z_vtk, z0_);

  vector<double> B_vtk = B0;
  //B_vtk = scale_vector(B_vtk, x0_);

  vector<double> rho_vtk = rho0;
  //rho_vtk = scale_vector(rho_vtk, rho0_);
  scale_vector(rho_vtk, rho0_);

  vector<double> u_vtk = u0;
  //u_vtk = scale_vector(u_vtk, u0_);
  scale_vector(u_vtk, u0_);
  ofstream ovel;
  ovel.open(name_);
  ovel << "# vtk DataFile Version 3.1 \n";
  ovel << "Fluid flow field from LBM simulation\n";
  ovel << "ASCII\n";
  ovel << "DATASET RECTILINEAR_GRID\n";
  ovel << "DIMENSIONS " << (nx-2*envelopeWidth)*decomposition[0] << " " << (ny-2*envelopeWidth)*decomposition[1] << " " << (nz-2*envelopeWidth)*decomposition[2] << "\n";
  //ovel << "DATASET UNSTRUCTURED_GRID\n\n";

  ovel << "X_COORDINATES " << (nx-2*envelopeWidth)*decomposition[0] << " float" << endl;
  for(int iproc = 0; iproc < decomposition[0]; iproc++){
    for(int i = envelopeWidth; i < nx-envelopeWidth; i++){
      int j = 0; int jproc = 0;
      int k = 0; int kproc = 0;
      //int index = i*ny + j + iproc*nx*ny*decomposition[1] + jproc*nx*ny;
      int index = i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*nz*decomposition[2] + kproc*nx*ny*nz;
      ovel << x_vtk[index] << "\n";
    }
  }

  ovel << "Y_COORDINATES " << (ny-2*envelopeWidth)*decomposition[1] << " float" << endl;
  for(int jproc = 0; jproc < decomposition[1]; jproc++){
    for(int j = envelopeWidth; j < ny-envelopeWidth; j++){
      int i = 0; int iproc = 0;
      int k = 0; int kproc = 0;
      //int index = i*ny + j + iproc*nx*ny*decomposition[1] + jproc*nx*ny;
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
        //int index = i*ny + j + iproc*nx*ny*decomposition[1] + jproc*nx*ny;
        int index = i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*nz*decomposition[2] + kproc*nx*ny*nz;
        ovel << z_vtk[index] << "\n";
      }
    }
  }

  // POINT_DATA
  // Density
  if(domain->dimension == 2){
    ovel << "\nPOINT_DATA " << (nx-envelopeWidth*2)*decomposition[0]*(ny-envelopeWidth*2)*decomposition[1] << "\n";
  }else{
    ovel << "\nPOINT_DATA " << (nx-envelopeWidth*2)*decomposition[0]*(ny-envelopeWidth*2)*decomposition[1]*(nz-envelopeWidth*2)*decomposition[2] << "\n";
  }
  ovel << "SCALARS B-Field FLOAT" << "\n";
  ovel << "LOOKUP_TABLE default" << "\n";
  for(int kproc = 0; kproc<decomposition[2]; kproc++){
    for(int k=envelopeWidth; k<(nz-envelopeWidth); k++){
      for(int jproc = 0; jproc<decomposition[1]; jproc++){
        for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
          for(int iproc = 0; iproc<decomposition[0]; iproc++){
            for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
//              int index = i*ny + j + iproc*nx*ny*decomposition[1] + jproc*nx*ny;
              int index = i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*decomposition[1] + kproc*nx*ny*nz;
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
    for(int k=envelopeWidth; k<(nz-envelopeWidth); k++){
      for(int jproc = 0; jproc<decomposition[1]; jproc++){
        for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
          for(int iproc = 0; iproc<decomposition[0]; iproc++){
            for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
              //int index = i*ny + j + iproc*nx*ny*decomposition[1] + jproc*nx*ny;
              int index = i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*decomposition[1] + kproc*nx*ny*nz;
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
    for(int k=envelopeWidth; k<(nz-envelopeWidth); k++){
      for(int jproc = 0; jproc<decomposition[1]; jproc++){
        for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
          for(int iproc = 0; iproc<decomposition[0]; iproc++){
            for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
              //int index = i*ny*2 + j*2 + iproc*nx*ny*decomposition[1]*2 + jproc*nx*ny*2;
              int index = (i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*decomposition[1] + kproc*nx*ny*nz)*3;
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
    for(int k=envelopeWidth; k<(nz-envelopeWidth); k++){
      for(int jproc = 0; jproc<decomposition[1]; jproc++){
        for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
          for(int iproc = 0; iproc<decomposition[0]; iproc++){
            for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
              //int index = i*ny*2 + j*2 + iproc*nx*ny*decomposition[1]*2 + jproc*nx*ny*2 + 1;
              int index = (i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*decomposition[1] + kproc*nx*ny*nz)*3+1;
              ovel << u_vtk[index] <<"\n";
            }
          }
        }
      }
    }
  }

  ovel << "\nSCALARS Velocity-z FLOAT" << "\n";
  ovel << "LOOKUP_TABLE default" << "\n";
  for(int kproc = 0; kproc<decomposition[2]; kproc++){
    for(int k=envelopeWidth; k<(nz-envelopeWidth); k++){
      for(int jproc = 0; jproc<decomposition[1]; jproc++){
        for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
          for(int iproc = 0; iproc<decomposition[0]; iproc++){
            for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
              //int index = i*ny*2 + j*2 + iproc*nx*ny*decomposition[1]*2 + jproc*nx*ny*2 + 1;
              int index = (i*ny*nz + j*nz + k + iproc*nx*ny*decomposition[1]*nz*decomposition[2] + jproc*nx*ny*decomposition[1] + kproc*nx*ny*nz)*3+2;
              ovel << u_vtk[index] <<"\n";
            }
          }
        }
      }
    }
  }

/*ovel << "Z_COORDINATES " << (nz-2*envelopeWidth)*decomposition[2] << " float" << endl;  
  for(int kproc = 0; jproc < decomposition[2]; kproc++){
    for(int k = envelopeWidth; k < nz-envelopeWidth; k++){
     // int i = 0; int iproc = 0;
      int index = i*ny + j + iproc*nx*ny*decomposition[1] + jproc*nx*ny;
      ovel << z_vtk[index] << "\n";
    }
  }
*/
  /*
  // POINTS
  int pnt_id = -1;
  ovel << "POINTS " << (nx-envelopeWidth*2)*decomposition[0]*(ny-envelopeWidth*2)*decomposition[1] << " FLOAT\n";
  for(int iproc = 0; iproc<decomposition[0]; iproc++){
    for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
      for(int jproc = 0; jproc<decomposition[1]; jproc++){
        for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
          int index = i*ny + j + iproc*nx*ny*decomposition[1] + jproc*nx*ny;
          ovel << x_vtk[index] << " " << y_vtk[index] << " 0" << "\n";
          pnt_id += 1;
          point_id[index] = pnt_id;
        }
      }
    }
  }
  // CELLS
  // Calculate no of cells
  int no_cells = ((nx-envelopeWidth*2)*decomposition[0]-1)*((ny-envelopeWidth*2)*decomposition[1]-1);
  int no_cell_blocks = no_cells*5;
  ovel << "\nCELLS " << no_cells << " " << no_cell_blocks;
  for(int jproc = 0; jproc<decomposition[1]; jproc++){
    for(int j=envelopeWidth; j<(ny-1-envelopeWidth); j++){
      for(int iproc = 0; iproc<decomposition[0]; iproc++){
        for(int i=envelopeWidth; i<(nx-1-envelopeWidth); i++){
          int index1 = i*ny + j + iproc*nx*ny*decomposition[1] + jproc*nx*ny;
          int index2 = (i+1)*ny + j + iproc*nx*ny*decomposition[1] + jproc*nx*ny;
          ovel << "\n4 " << point_id[index1] << " " << point_id[index2] << " " << point_id[index2+1] << " " << point_id[index1+1];
          //ovel << "\n4 " << point_id[i*ny*decomposition[1]+j] << " " << point_id[(i+1)*ny*decomposition[1]+j] << " " << point_id[(i+1)*ny*decomposition[1]+(j+1)] << " " << point_id[i*ny*decomposition[1]+(j+1)];
        }
      }
    }
  }
  // CELL_TYPES
  ovel << "\n\nCELL_TYPES " << no_cells << "\n";
  for(int i = 0; i<no_cells; i++)
    ovel << "9 ";*/
/*  // POINT_DATA
  // Density
  ovel << "\n\nPOINT_DATA " << (nx-envelopeWidth*2)*decomposition[0]*(ny-envelopeWidth*2)*decomposition[1] << "\n";
  ovel << "SCALARS B-Field FLOAT" << "\n";
  ovel << "LOOKUP_TABLE default" << "\n";
  for(int iproc = 0; iproc<decomposition[0]; iproc++){
    for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
      for(int jproc = 0; jproc<decomposition[1]; jproc++){
        for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
          int index = i*ny + j + iproc*nx*ny*decomposition[1] + jproc*nx*ny;
          ovel << B_vtk[index] <<"\n";
        }
      }
    }
  }

  ovel << "\nSCALARS Density FLOAT" << "\n";
  ovel << "LOOKUP_TABLE default" << "\n";
  for(int iproc = 0; iproc<decomposition[0]; iproc++){
    for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
      for(int jproc = 0; jproc<decomposition[1]; jproc++){
        for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
          int index = i*ny + j + iproc*nx*ny*decomposition[1] + jproc*nx*ny;
          ovel << rho_vtk[index] << "\n";
        }
      }
    }
  }

  ovel << "\nSCALARS Velocity-x FLOAT" << "\n";
  ovel << "LOOKUP_TABLE default" << "\n";
  for(int iproc = 0; iproc<decomposition[0]; iproc++){
    for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
      for(int jproc = 0; jproc<decomposition[1]; jproc++){
        for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
          int index = i*ny*2 + j*2 + iproc*nx*ny*decomposition[1]*2 + jproc*nx*ny*2;
          ovel << u_vtk[index] <<"\n";
        }
      }
    }
  }

  ovel << "\nSCALARS Velocity-y FLOAT" << "\n";
  ovel << "LOOKUP_TABLE default" << "\n";
  for(int iproc = 0; iproc<decomposition[0]; iproc++){
    for(int i=envelopeWidth; i<(nx-envelopeWidth); i++){
      for(int jproc = 0; jproc<decomposition[1]; jproc++){
        for(int j=envelopeWidth; j<(ny-envelopeWidth); j++){
          int index = i*ny*2 + j*2 + iproc*nx*ny*decomposition[1]*2 + jproc*nx*ny*2 + 1;
          ovel << u_vtk[index] <<"\n";
        }
      }
    }
  }
*/

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

//vector<double> WriteVTK::scale_vector(vector<double> &vec_, doublescaling_){
void WriteVTK::scale_vector(vector<double> &vec_, double scaling_){
  std::transform(vec_.begin(), vec_.end(), vec_.begin(),
               std::bind(std::multiplies<double>(), std::placeholders::_1, scaling_));
  //return vec_;
};
