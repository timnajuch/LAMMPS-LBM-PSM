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
    int nx, ny, nz;
    int q;

    vector<int> nxLocal;
    vector<int> nyLocal;
    vector<int> nzLocal;
    vector<int> nxLocalGrid;
    vector<int> nyLocalGrid;
    vector<int> nzLocalGrid;

    int envelopeWidth;
    int dimension;

    vector<double> w;
    vector<double> g;
    vector<double> e;

    double c, cs , csPow2;
    
    double dx, dy, dz;

  
    vector<double> f;
    vector<double> f0;
    vector<double> fcoll;


    vector<double> procOrigin;
    vector<double> procLength;
    vector<double> x, y, z;
    vector<double> rho;
    vector<double> B;
    
    vector<double> u;
    vector<double> us;

    vector<double> Fhydx;
    vector<double> Fhydy;
    vector<double> Fhydz;

    vector<ParticleDataOnLattice> pData;

    void set_f(int i_, int j_, int k_, int iq_, double value_);
    void set_f(int ind_iq_, double value_);
    void cp_fcoll_f();
    double get_f(int i_, int j_, int k_, int iq_);
    double get_f(int ind_iq_);

    void set_f0(int i_, int j_, int k_, int iq_, double value_);
    double get_f0(int i_, int j_, int k_, int iq_);
    double get_f0(int ind_iq_);

    void set_fcoll(int i_, int j_, int k_, int iq_, double value_);
    double get_fcoll(int i_, int j_, int k_, int iq_);
    double get_fcoll(int ind_iq_);

    int decomposition[3];
    int procCoordinates[3];

  public:
    LBMPSMLattice(int nx_, int ny_, int nz_, int q_, int decomposition[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_, double dx);
    ~LBMPSMLattice();

    void initialise_channel_geometry(double wallHeightHalf_, double eps_, double nychannel_, double dx_, double dy_, double dz_);
    void initialise_domain(double dx_, double dy_, double dz_);

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

    vector<double> getProcOrigin(){ return procOrigin; }
    vector<double> getProcLength(){ return procLength; }
    int* getProcCoordinates(){ return procCoordinates; }

    int get_nx();
    int get_ny();
    int get_nz();
    int get_envelopeWidth(){ return envelopeWidth; }

    int get_nxLocal(int iProcIndex);
    int get_nyLocal(int jProcIndex);
    int get_nzLocal(int kProcIndex);

    vector<double>& getVector_f(){ return f; }

    void set_B(int index, double B_);
    double get_B(int index);
    
    double get_rho(int index);
    double get_u(int index);

    void setParticleOnLattice(int index, LAMMPS_NS::tagint pID, double uP[3], double eps);
    void setToZero(int index, LAMMPS_NS::tagint pID);
    double getSolidFractionOnLattice(int index, int pID);
    vector<double> getSolidVelocityOnLattice(int index, int pID);
    vector<double> getSolidVelocityOnLattice(int index);

    ParticleDataOnLattice getParticleDataOnLatticeNode(int index);

    void set_Fhydx(int index, int pID, double Fhydx_);
    void set_Fhydy(int index, int pID, double Fhydy_);
    void set_Fhydz(int index, int pID, double Fhydz_);
    void add_Fhydx(int index, int pID, double Fhydx_);
    void add_Fhydy(int index, int pID, double Fhydy_);
    void add_Fhydz(int index, int pID, double Fhydz_);
    double get_Fhydx(int index);
    double get_Fhydy(int index);
    double get_Fhydz(int index);
 
    vector<int> get_procCoordinates();
   
    friend class ZouHeBC;
};

#endif
