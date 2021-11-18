/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#ifndef LATTICE_2D_H
#define LATTICE_2D_H

#include <iostream>
#include <math.h>
#include <vector>

#include "lmptype.h"

#include "particleDataOnLattice.h"

using namespace std;

class Lattice2D{

  protected:
    int nx, ny, nz;
    int q;

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
    Lattice2D(int nx_, int ny_, int nz_, int q_, int decomposition[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_, double dx);
    ~Lattice2D();

    void initialise_channel_geometry(double wallHeightHalf_, double eps_, double nychannel_, double dx_, double dy_, double dz_);
    void initialise_domain(double dx_, double dy_, double dz_);

    vector<double> get_B();
    vector<double> get_rho();
    vector<double> get_x();
    vector<double> get_y();
    vector<double> get_z();
    vector<double> get_u();

    vector<double> getProcOrigin(){ return procOrigin; }
    vector<double> getProcLength(){ return procLength; }
    int* getProcCoordinates(){ return procCoordinates; }

    int get_nx();
    int get_ny();
    int get_nz();
    int get_envelopeWidth(){ return envelopeWidth; }

    vector<double>& getVector_f(){ return f; }

    void set_B(int index, double B_);
    double get_B(int index);
    
    double get_rho(int index);

    //void setParticleOnLattice(int index, LAMMPS_NS::tagint pID, double uP[2], double eps);
    void setParticleOnLattice(int index, LAMMPS_NS::tagint pID, double uP[3], double eps);
    double getSolidFractionOnLattice(int index, int pID);
    vector<double> getSolidVelocityOnLattice(int index, int pID);

    ParticleDataOnLattice getParticleDataOnLatticeNode(int index);

    void set_Fhydx(int index, double Fhydx_);
    void set_Fhydy(int index, double Fhydy_);
    void set_Fhydz(int index, double Fhydz_);
    void add_Fhydx(int index, double Fhydx_);
    void add_Fhydy(int index, double Fhydy_);
    void add_Fhydz(int index, double Fhydz_);
    double get_Fhydx(int index);
    double get_Fhydy(int index);
    double get_Fhydz(int index);
    
    friend class ZouHeBC2D;
};

#endif
