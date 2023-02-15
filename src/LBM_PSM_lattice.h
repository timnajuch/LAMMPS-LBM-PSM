/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch, 2022
------------------------------------------------------*/

#ifndef LBM_PSM_LATTICE_2D_H
#define LBM_PSM_LATTICE_2D_H

#include <iostream>
#include <math.h>
#include <vector>

#include "lmptype.h"

#include "LBM_PSM_particleDataOnLattice.h"
#include "LBM_PSM_MPICOMM.h"

using namespace std;

class LBMPSMLattice{

  protected:
    int nx, ny, nz;                   // Global number of lattice nodes
  
    vector<int> nxLocal;              // Local (on processor) number of lattice nodes in x-direction. Processor id is continuous over all dimensions (i.e. procIndex = iproc*decomposition[1]*decomposition[2] + jproc*decomposition[2] + kproc; )
    vector<int> nyLocal;              // Local (on processor) number of lattice nodes in y-direction. Processor id is continuous over all dimensions
    vector<int> nzLocal;              // Local (on processor) number of lattice nodes in z-direction. Processor id is continuous over all dimensions
    vector<int> nxLocalGrid;          // Local (on processor) number of lattice nodes in x-direction. Processor id defines processor in grid location
    vector<int> nyLocalGrid;          // Local (on processor) number of lattice nodes in y-direction. Processor id defines processor in grid location
    vector<int> nzLocalGrid;          // Local (on processor) number of lattice nodes in z-direction. Processor id defines processor in grid location

    int envelopeWidth;                // Envelope / halo width defining the number of lattice nodes which are communicated between processors via MPI
    int dimension;                    // 2D or 3D system

    int q;                            // Defines the number of velocities in LB (e.g. D2Q9 (q = 9) or D3Q19 (q = 19))
    vector<double> w;                 // Weighting factor in LB equilibrium function
    vector<double> e;                 // LB velocity sets

    double cs, csPow2, csPow4;        // Speed of sound in LBM and to power of 2/4
    double invCsPow2, invCsPow4;      // Inverse of speed of sound to power of 2/4

    double dx, dy, dz;                // Lattice discretisation in 3D
  
    int currentStep, nextStep;        // Populations of two different time steps are stored in same std::vector at different places in memory. Need this variables to keep track

    vector<double> f;                 // Populations 
    vector<double> f0;                // Equilibrium populations

    vector<double> origin_global;     // Global origin of lattice 
    vector<double> x, y, z;           // Coordinates of lattice nodes 
    vector<double> rho;               // Density 
    vector<double> B;                 // Weighting factor which is based on the solid fraction
    
    vector<double> u;                 // Fluid velocity
    vector<double> us;                // Solid velocity

    vector<ParticleDataOnLattice> pData;  // Particle data stored on lattice node

    void set_f(int i_, int j_, int k_, int iq_, int step_, double value_);
    void set_f(int ind_iq_, double value_);
    double get_f(int i_, int j_, int k_, int iq_, int step_);
    double get_f(int ind_iq_);

    void set_f0(int i_, int j_, int k_, int iq_, double value_);
    double get_f0(int i_, int j_, int k_, int iq_);
    double get_f0(int ind_iq_);

    int decomposition[3];
    int procCoordinates[3];

  public:
    LBMPSMLattice(int nx_, int ny_, int nz_, int decomposition[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_);
    ~LBMPSMLattice();

    void initialise_domain(double dx_, double dy_, double dz_);

    int index_1D(int i, int j, int k);
    int index_2D(int i, int j, int k, int direction);
    int index_fi(int i, int j, int k, int iq, int step);

    int get_currentStep();
    void set_currentStep(int currentStep);

    vector<double> get_B();
    vector<double> get_rho();
    vector<double> get_x();
    vector<double> get_y();
    vector<double> get_z();
    vector<double> get_u();

    vector<double>& get_B_reference();
    vector<double>& get_rho_reference();
    vector<double>& get_x_reference();
    vector<double>& get_y_reference();
    vector<double>& get_z_reference();
    vector<double>& get_u_reference();

    int get_nx();
    int get_ny();
    int get_nz();
    int get_envelopeWidth();
    int get_q();

    int get_nxLocal(int iProcIndex);
    int get_nyLocal(int jProcIndex);
    int get_nzLocal(int kProcIndex);

    void set_B(int index, double B_);
    double get_B(int index);

    double get_rho(int index);
    double get_u(int index);
    double get_u_at_node(int index_node_1D, int direction);

    vector<double>& getVector_f();
    void setVector_f(vector<double>& fcopy); // Function used to copy the population values which are read from a restart file

    ParticleDataOnLattice getParticleDataOnLatticeNode(int index);
    void setParticleOnLattice(int index, LAMMPS_NS::tagint pID, double uP[3], double eps);
    void setToZero(int index, LAMMPS_NS::tagint pID);
    double getSolidFractionOnLattice(int index, int pID);
    vector<double> getSolidVelocityOnLattice(int index);
    vector<double> getSolidVelocityOnLattice(int index, int pID);
    void add_Fhyd(int index, LAMMPS_NS::tagint pID, double Fhyd, int dir);

    vector<int> get_procCoordinates();

    friend class ZouHeBC;
};

#endif