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
}


void WriteVTK::pre_force(int)
{
  if (update->ntimestep % nevery) return;
  double u_infty = fixLBMPSM->unitConversion->get_u_lb();
  double Uc = fixLBMPSM->unitConversion->get_Uc();

  ostringstream timeStringTmp;
  timeStringTmp << "flowField_" << setw(10) << setfill('0') << update->ntimestep << ".vtk";
  string timeString(timeStringTmp.str());

  write_vtk(timeString,
            fixLBMPSM->dynamics->get_x_reference(), 1.0,
            fixLBMPSM->dynamics->get_y_reference(), 1.0,
            fixLBMPSM->dynamics->get_z_reference(), 1.0,
            fixLBMPSM->dynamics->get_B_reference(), 1.0,
            fixLBMPSM->dynamics->get_rho_reference(), fixLBMPSM->get_rho(),
            fixLBMPSM->dynamics->get_u_reference(), Uc/u_infty);

}


void WriteVTK::write_vtk(string name_, vector<double> &x_, double x0_, vector<double> &y_, double y0_, vector<double> &z_, double z0_, vector<double> &B_, double B0_, vector<double> &rho_, double rho0_, vector<double> &u_, double u0_){
  int envelopeWidth = 1;
  int nzLoopStart = 0;
  int nzLoopEnd = 1;
  if(domain->dimension == 3){
    nzLoopStart = envelopeWidth;
    nzLoopEnd = nz-envelopeWidth;
  }

  vector<double> B_clip;
  B_clip.resize(nx*ny*nz);
  for (int i = 0; i < nx*ny*nz; ++i){
    B_clip[i] = fixLBMPSM->dynamics->getParticleDataOnLatticeNode(i).solidFraction[0] + fixLBMPSM->dynamics->getParticleDataOnLatticeNode(i).solidFraction[1];
    if (B_clip[i] > 1.0)
      { B_clip[i] = 1.0; }
  }

  // Variables for MPI communication
  int *MPIGathervCounts; // recvcounts definition / parameter in for MPI_Gatherv(...)
  int *MPIGathervDispls; // displs definition for / parameter in MPI_Gatherv(...)
  int *MPIGathervCounts3D; // recvcounts definition / parameter in for MPI_Gatherv(...)
  int *MPIGathervDispls3D; // displs definition for / parameter in MPI_Gatherv(...)
 
  MPIGathervCounts = new int[comm->nprocs];
  MPIGathervDispls = new int[comm->nprocs];
  MPIGathervCounts3D = new int[comm->nprocs];
  MPIGathervDispls3D = new int[comm->nprocs];

  int nelements = nx*ny*nz;
  int nelements3D = nx*ny*nz*3;
  // Each process tells the root how many elements it holds
  MPI_Gather(&nelements, 1, MPI_INT, MPIGathervCounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&nelements3D, 1, MPI_INT, MPIGathervCounts3D, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Displacements in the receive buffer for MPI_GATHERV
  for (int i = 0; i < comm->nprocs; i++){
    MPIGathervDispls[i] = (i > 0) ? (MPIGathervDispls[i-1] + MPIGathervCounts[i-1]) : 0;
    MPIGathervDispls3D[i] = (i > 0) ? (MPIGathervDispls3D[i-1] + MPIGathervCounts3D[i-1]) : 0;
  }

  vector<double> x0;
  vector<double> y0;
  vector<double> z0;
  vector<double> B0;
  vector<double> rho0;
  vector<double> u0;
  if (comm->me == 0){
    x0.resize(MPIGathervDispls[comm->nprocs-1]+MPIGathervCounts[comm->nprocs-1]);
    y0.resize(MPIGathervDispls[comm->nprocs-1]+MPIGathervCounts[comm->nprocs-1]);
    z0.resize(MPIGathervDispls[comm->nprocs-1]+MPIGathervCounts[comm->nprocs-1]);
    B0.resize(MPIGathervDispls[comm->nprocs-1]+MPIGathervCounts[comm->nprocs-1]);
    rho0.resize(MPIGathervDispls[comm->nprocs-1]+MPIGathervCounts[comm->nprocs-1]);
    u0.resize(MPIGathervDispls3D[comm->nprocs-1]+MPIGathervCounts3D[comm->nprocs-1]);
  }

  MPI_Gatherv(&x_[0], nelements, MPI_DOUBLE, &x0[0], MPIGathervCounts, MPIGathervDispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(&y_[0], nelements, MPI_DOUBLE, &y0[0], MPIGathervCounts, MPIGathervDispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(&z_[0], nelements, MPI_DOUBLE, &z0[0], MPIGathervCounts, MPIGathervDispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(&B_clip[0], nelements, MPI_DOUBLE, &B0[0], MPIGathervCounts, MPIGathervDispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(&rho_[0], nelements, MPI_DOUBLE, &rho0[0], MPIGathervCounts, MPIGathervDispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(&u_[0], nelements3D, MPI_DOUBLE, &u0[0], MPIGathervCounts3D, MPIGathervDispls3D, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int nxTotal = 0;
  for(int iproc = 0; iproc < decomposition[0]; iproc++){
      nxTotal += fixLBMPSM->dynamics->get_nxLocal(iproc); }

  int nyTotal = 0;
    for(int jproc = 0; jproc < decomposition[1]; jproc++){
      nyTotal += fixLBMPSM->dynamics->get_nyLocal(jproc); }

  int nzTotal = 0;
    for(int kproc = 0; kproc < decomposition[2]; kproc++){
      nzTotal += fixLBMPSM->dynamics->get_nzLocal(kproc); }

  if(comm->me == 0){
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
      ovel << "DIMENSIONS " << nxTotal-2*envelopeWidth*decomposition[0] << " " << nyTotal-2*envelopeWidth*decomposition[1] << " " << 1 << "\n";
    }else{
      ovel << "DIMENSIONS " << nxTotal-2*envelopeWidth*decomposition[0] << " " << nyTotal-2*envelopeWidth*decomposition[1] << " " << nzTotal-2*envelopeWidth*decomposition[2] << "\n";
    }

    ovel << "X_COORDINATES " << nxTotal-2*envelopeWidth*decomposition[0] << " float" << endl;
    for(int iproc = 0; iproc < decomposition[0]; iproc++){
      for(int i = envelopeWidth; i < fixLBMPSM->dynamics->get_nxLocal(iproc)-envelopeWidth; i++){
        int j = 0; int jproc = 0;
        int k = 0; int kproc = 0;
        int procIndex = iproc*decomposition[1]*decomposition[2] + jproc*decomposition[2] + kproc;
        int index = i*fixLBMPSM->dynamics->get_nyLocal(jproc)*fixLBMPSM->dynamics->get_nzLocal(kproc) + j*fixLBMPSM->dynamics->get_nzLocal(kproc) + k + MPIGathervDispls[procIndex];
        ovel << x_vtk[index] << "\n";
      }
    }

    ovel << "Y_COORDINATES " << nyTotal-2*envelopeWidth*decomposition[1] << " float" << endl;
    for(int jproc = 0; jproc < decomposition[1]; jproc++){
      for(int j = envelopeWidth; j < fixLBMPSM->dynamics->get_nyLocal(jproc)-envelopeWidth; j++){
        int i = 0; int iproc = 0;
        int k = 0; int kproc = 0;
        int procIndex = iproc*decomposition[1]*decomposition[2] + jproc*decomposition[2] + kproc;
        int index = i*fixLBMPSM->dynamics->get_nyLocal(jproc)*fixLBMPSM->dynamics->get_nzLocal(kproc) + j*fixLBMPSM->dynamics->get_nzLocal(kproc) + k + MPIGathervDispls[procIndex];
        ovel << y_vtk[index] << "\n";
      }
    }

    if(domain->dimension == 2){
      ovel << "Z_COORDINATES " << 1 << " float" << endl;
      ovel << "0.0\n";
    }else{
      ovel << "Z_COORDINATES " << nzTotal-2*envelopeWidth*decomposition[2] << " float" << endl;
      for(int kproc = 0; kproc < decomposition[2]; kproc++){
        for(int k = envelopeWidth; k < fixLBMPSM->dynamics->get_nzLocal(kproc)-envelopeWidth; k++){
          int i = 0; int iproc = 0;
          int j = 0; int jproc = 0;
          int procIndex = iproc*decomposition[1]*decomposition[2] + jproc*decomposition[2] + kproc;
          int index = i*fixLBMPSM->dynamics->get_nyLocal(jproc)*fixLBMPSM->dynamics->get_nzLocal(kproc) + j*fixLBMPSM->dynamics->get_nzLocal(kproc) + k + MPIGathervDispls[procIndex];
          ovel << z_vtk[index] << "\n";
        }
      }
    }

    // POINT_DATA
    if(domain->dimension == 2){
      ovel << "\nPOINT_DATA " << (nxTotal-envelopeWidth*2*decomposition[0])*(nyTotal-envelopeWidth*2*decomposition[1]) << "\n";
    }else{
      ovel << "\nPOINT_DATA " << (nxTotal-envelopeWidth*2*decomposition[0])*(nyTotal-envelopeWidth*2*decomposition[1])*(nzTotal-envelopeWidth*2*decomposition[2]) << "\n";
    }
    ovel << "SCALARS SolidFraction FLOAT" << "\n";
    ovel << "LOOKUP_TABLE default" << "\n";
    for(int kproc = 0; kproc<decomposition[2]; kproc++){
      if(domain->dimension == 3){
        nzLoopStart = envelopeWidth;
        nzLoopEnd = fixLBMPSM->dynamics->get_nzLocal(kproc)-envelopeWidth;
      }
      for(int k=nzLoopStart; k<nzLoopEnd; k++){
        for(int jproc = 0; jproc<decomposition[1]; jproc++){
          for(int j=envelopeWidth; j<(fixLBMPSM->dynamics->get_nyLocal(jproc)-envelopeWidth); j++){
            for(int iproc = 0; iproc<decomposition[0]; iproc++){
              for(int i=envelopeWidth; i<(fixLBMPSM->dynamics->get_nxLocal(iproc)-envelopeWidth); i++){
                int procIndex = iproc*decomposition[1]*decomposition[2] + jproc*decomposition[2] + kproc;
                int index = i*fixLBMPSM->dynamics->get_nyLocal(jproc)*fixLBMPSM->dynamics->get_nzLocal(kproc) + j*fixLBMPSM->dynamics->get_nzLocal(kproc) + k + MPIGathervDispls[procIndex];
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
      if(domain->dimension == 3){
        nzLoopStart = envelopeWidth;
        nzLoopEnd = fixLBMPSM->dynamics->get_nzLocal(kproc)-envelopeWidth;
      }
      for(int k=nzLoopStart; k<nzLoopEnd; k++){
        for(int jproc = 0; jproc<decomposition[1]; jproc++){
          for(int j=envelopeWidth; j<(fixLBMPSM->dynamics->get_nyLocal(jproc)-envelopeWidth); j++){
            for(int iproc = 0; iproc<decomposition[0]; iproc++){
              for(int i=envelopeWidth; i<(fixLBMPSM->dynamics->get_nxLocal(iproc)-envelopeWidth); i++){
                int procIndex = iproc*decomposition[1]*decomposition[2] + jproc*decomposition[2] + kproc;
                int index = i*fixLBMPSM->dynamics->get_nyLocal(jproc)*fixLBMPSM->dynamics->get_nzLocal(kproc) + j*fixLBMPSM->dynamics->get_nzLocal(kproc) + k + MPIGathervDispls[procIndex];
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
      if(domain->dimension == 3){
        nzLoopStart = envelopeWidth;
        nzLoopEnd = fixLBMPSM->dynamics->get_nzLocal(kproc)-envelopeWidth;
      }
      for(int k=nzLoopStart; k<nzLoopEnd; k++){
        for(int jproc = 0; jproc<decomposition[1]; jproc++){
          for(int j=envelopeWidth; j<(fixLBMPSM->dynamics->get_nyLocal(jproc)-envelopeWidth); j++){
            for(int iproc = 0; iproc<decomposition[0]; iproc++){
              for(int i=envelopeWidth; i<(fixLBMPSM->dynamics->get_nxLocal(iproc)-envelopeWidth); i++){
                int procIndex = iproc*decomposition[1]*decomposition[2] + jproc*decomposition[2] + kproc;
                int index = (i*fixLBMPSM->dynamics->get_nyLocal(jproc)*fixLBMPSM->dynamics->get_nzLocal(kproc) + j*fixLBMPSM->dynamics->get_nzLocal(kproc) + k)*3 + MPIGathervDispls3D[procIndex];
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
      if(domain->dimension == 3){
        nzLoopStart = envelopeWidth;
        nzLoopEnd = fixLBMPSM->dynamics->get_nzLocal(kproc)-envelopeWidth;
      }
      for(int k=nzLoopStart; k<nzLoopEnd; k++){
        for(int jproc = 0; jproc<decomposition[1]; jproc++){
          for(int j=envelopeWidth; j<(fixLBMPSM->dynamics->get_nyLocal(jproc)-envelopeWidth); j++){
            for(int iproc = 0; iproc<decomposition[0]; iproc++){
              for(int i=envelopeWidth; i<(fixLBMPSM->dynamics->get_nxLocal(iproc)-envelopeWidth); i++){
                int procIndex = iproc*decomposition[1]*decomposition[2] + jproc*decomposition[2] + kproc;
                int index = (i*fixLBMPSM->dynamics->get_nyLocal(jproc)*fixLBMPSM->dynamics->get_nzLocal(kproc) + j*fixLBMPSM->dynamics->get_nzLocal(kproc) + k)*3 + 1 + MPIGathervDispls3D[procIndex];
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
      if(domain->dimension == 3){
        nzLoopStart = envelopeWidth;
        nzLoopEnd = fixLBMPSM->dynamics->get_nzLocal(kproc)-envelopeWidth;
      }
      for(int k=nzLoopStart; k<nzLoopEnd; k++){
        for(int jproc = 0; jproc<decomposition[1]; jproc++){
          for(int j=envelopeWidth; j<(fixLBMPSM->dynamics->get_nyLocal(jproc)-envelopeWidth); j++){
            for(int iproc = 0; iproc<decomposition[0]; iproc++){
              for(int i=envelopeWidth; i<(fixLBMPSM->dynamics->get_nxLocal(iproc)-envelopeWidth); i++){
                  int procIndex = iproc*decomposition[1]*decomposition[2] + jproc*decomposition[2] + kproc;
                  int index = (i*fixLBMPSM->dynamics->get_nyLocal(jproc)*fixLBMPSM->dynamics->get_nzLocal(kproc) + j*fixLBMPSM->dynamics->get_nzLocal(kproc) + k)*3 + 2 + MPIGathervDispls3D[procIndex];
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



void WriteVTK::scale_vector(vector<double> &vec_, double scaling_){
  std::transform(vec_.begin(), vec_.end(), vec_.begin(),
               std::bind(std::multiplies<double>(), std::placeholders::_1, scaling_));
};
