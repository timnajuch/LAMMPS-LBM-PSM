/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch
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
    double* w;                        // Weighting factor in LB equilibrium function
    int* ex;                          // LB velocity set in x-direction
    int* ey;                          // LB velocity set in y-direction
    int* ez;                          // LB velocity set in z-direction
    

    double cs, csPow2, csPow4;        // Speed of sound in LBM and to power of 2/4
    double invCsPow2, invCsPow4;      // Inverse of speed of sound to power of 2/4

    double dx, dy, dz;                // Lattice discretisation in 3D
  
    vector<double> f_curr;            // Populations at current time step
    vector<double> f_next;            // Populations at next time step (ping-pong buffer)
    vector<double> f0;                // Equilibrium populations

    vector<double> origin_global;     // Global origin of lattice 
    vector<double> x, y, z;           // Coordinates of lattice nodes 
    vector<double> rho;               // Density 
    vector<double> B;                 // Weighting factor which is based on the solid fraction
    
    vector<double> u;                 // Fluid velocity
    vector<double> us;                // Solid velocity

    vector<ParticleDataOnLattice> pData;  // Particle data stored on lattice node

    void setLattice2D();
    void setLattice3D();

    void set_f(int ind_iq_, double value_);
    double get_f(int ind_iq_);

    void set_f0(int i_, int j_, int k_, int iq_, double value_);
    void set_f0(int ind_iq_, double value_);
    double get_f0(int i_, int j_, int k_, int iq_);
    double get_f0(int ind_iq_);

    int decomposition[3];
    int procCoordinates[3];

  public:
    LBMPSMLattice(int nx_, int ny_, int nz_, int decomposition[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_);
    ~LBMPSMLattice();

    void initialise_domain(double dx_, double dy_, double dz_);

    inline int index_1D(int i, int j, int k){ return i*ny*nz + j*nz + k; }
    inline int index_2D(int i, int j, int k, int direction){ return (i*ny*nz + j*nz + k)*3 + direction; }
    inline int index_fi(int i, int j, int k, int iq){ return (i*ny*nz + j*nz + k)*q + iq; }

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

    vector<double>& getVector_f_curr();
    void setVector_f_curr(vector<double>& fcopy); // Function used to copy the population values which are read from a restart file

    ParticleDataOnLattice getParticleDataOnLatticeNode(int index);
    inline ParticleDataOnLattice& getReferenceParticleDataOnLatticeNode(int index);
    inline const ParticleDataOnLattice& getReferenceParticleDataOnLatticeNode(int index) const;
    void setParticleOnLattice(int index, LAMMPS_NS::tagint pID, double uP[3], double eps);
    void setToZero(int index, LAMMPS_NS::tagint pID);
    double getSolidFractionOnLattice(int index, int pID);
    //vector<double> getSolidVelocityOnLattice(int index);
    //vector<double> getSolidVelocityOnLattice(int index, int pID);
    void add_Fhyd(int index, LAMMPS_NS::tagint pID, double Fhyd, int dir);

    vector<int> get_procCoordinates();

    friend class ZouHeBC;
};

/*
// Perhaps use in future. But need to refactor LBMPSMLattice class as template class (incl. children) and check other classes (e.g. boundaries) calling lattice class
struct D2Q9 {
    static constexpr int q = 9;

    static constexpr int ex[q] = {0,1,0,-1,0,1,-1,-1,1};
    static constexpr int ey[q] = {0,0,1,0,-1,1,1,-1,-1};
    static constexpr int ez[q] = {0,0,0,0,0,0,0,0,0};

    static constexpr double w[q] = {
        4.0/9.0,
        1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,
        1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0
    };
};


struct D3Q19 {
    static constexpr int q = 19;

    static constexpr int ex[q] = {
         0, 1,-1, 0, 0, 0, 0,
         1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0
    };

    static constexpr int ey[q] = {
         0, 0, 0, 1,-1, 0, 0,
         1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1
    };

    static constexpr int ez[q] = {
         0, 0, 0, 0, 0, 1,-1,
         0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1
    };

    static constexpr double w[q] = {
         1.0/3.0,
         1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0,
         1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
         1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
         1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
    };
};
*/

struct LatticeD2Q9 {
    int q = 9; // non-constant, runtime visible if needed
    int ex[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    int ey[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    int ez[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}; // can be ignored for 2D
    double w[9] = {
        4.0/9.0,
        1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,
        1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0
    };
};

struct LatticeD3Q19 {
    int q = 19;
    int ex[19] = {0,1,-1,0,0,0,0,1,-1,1,-1,0,0,1,-1,1,-1,0,0};
    int ey[19] = {0,0,0,1,-1,0,0,1,-1,0,0,1,-1,-1,1,0,0,1,-1};
    int ez[19] = {0,0,0,0,0,1,-1,0,0,1,-1,1,-1,0,0,-1,1,-1,1};
    double w[19] = {
        1.0/3.0,
        1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,
        1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,
        1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0
    };
};


inline ParticleDataOnLattice& LBMPSMLattice::getReferenceParticleDataOnLatticeNode(int index){ return pData[index]; }

inline const ParticleDataOnLattice& LBMPSMLattice::getReferenceParticleDataOnLatticeNode(int index) const { return pData[index]; }


#endif
