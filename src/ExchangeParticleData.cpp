/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2019
------------------------------------------------------*/

#include "ExchangeParticleData.h"

//ExchangeParticleData::ExchangeParticleData(double dp_, std::vector<double> xp_, std::vector<double> us_) : dp(dp_), xp(xp_), us(us_) {};
ExchangeParticleData::ExchangeParticleData() {};

ExchangeParticleData::~ExchangeParticleData() {};


void ExchangeParticleData::setParticlesOnLattice(Lattice2D *lattice2D_){
  for(int i = 0; i < lattice2D_->get_nx(); ++i){
    for(int j = 0; j < lattice2D_->get_ny(); ++j){

      int ind_phys_1D = i * lattice2D_->get_ny() + j;
      
      lattice2D_->set_B(ind_phys_1D, calcSolidFraction(i, j, xp[0], xp[1], dp/2.0));

    }
  }
};


void ExchangeParticleData::setParticlesOnLattice(Lattice2D *lattice2D_, Unit_Conversion *unitConversion, int numberParticles, LAMMPS_NS::tagint *tag, double **xPart, double **uPart, double *rp, vector<double> boxLength, vector<double> origin)
{
  for(int iPart = 0; iPart < numberParticles; ++iPart){
      double x_lb = fmin(fmax(unitConversion->get_pos_lb(xPart[iPart][0]-(lattice2D_->getProcOrigin()[0]-(double)lattice2D_->get_envelopeWidth()*unitConversion->get_dx())), 0.0), (double)lattice2D_->get_nx());
      double y_lb = fmin(fmax(unitConversion->get_pos_lb(xPart[iPart][1]-(lattice2D_->getProcOrigin()[1]-(double)lattice2D_->get_envelopeWidth()*unitConversion->get_dx())), 0.0), (double)lattice2D_->get_ny());
      double r_lb = unitConversion->get_radius_lb(*rp);

      int nodeZone[2][2] = {{(int)(x_lb-r_lb)-3, (int)(x_lb+r_lb)+3}, {(int)(y_lb-r_lb)-3, (int)(y_lb+r_lb)+3}};
      for (int i = 0; i < 2; i++){
          nodeZone[0][i] = fmin(fmax(nodeZone[0][i],0),lattice2D_->get_nx());
          nodeZone[1][i] = fmin(fmax(nodeZone[1][i],0),lattice2D_->get_ny());
      }
      
      for(int i = nodeZone[0][0]; i < nodeZone[0][1]; ++i){
        for(int j = nodeZone[1][0]; j < nodeZone[1][1]; ++j){

          int ind_phys_1D = i * lattice2D_->get_ny() + j;

          lattice2D_->setParticleOnLattice(ind_phys_1D, tag[iPart], uPart[iPart], calcSolidFraction(i, j, x_lb, y_lb, r_lb));

        }
      }
  }
};



double ExchangeParticleData::calcSolidFraction(int i, int j, double xP_LB, double yP_LB, double rP_LB){

    int slicesPerDim = 10;
    double sliceWidth = 1.0/((double)slicesPerDim);
    double fraction = 1.0/((double)(slicesPerDim*slicesPerDim));

    double sqrt2half = sqrt(2.0)/2.0;
    double dx = (double)i - xP_LB;
    double dy = (double)j - yP_LB;

    double dist = dx*dx + dy*dy;

    double rP = rP_LB + sqrt2half;
    if(dist > rP*rP) return 0.0;

    double rM = rP_LB - sqrt2half;
    if(dist < rM*rM) return  1.0;

    double rSq = rP_LB*rP_LB;
    double dx_sq[slicesPerDim], dy_sq[slicesPerDim];

    for(int i = 0; i < slicesPerDim; ++i){
      double delta = -0.5 + ((double)i+0.5)*sliceWidth;
      double dx_i = dx+delta; dx_sq[i] = dx_i*dx_i;
      double dy_i = dy+delta; dy_sq[i] = dy_i*dy_i;
    }   

    int n = 0;
    for(int i = 0; i < slicesPerDim; ++i){
      for(int j = 0; j < slicesPerDim; ++j){
        if(dx_sq[i] + dy_sq[j] < rSq) ++n;
      }   
    }   

    return fraction*((double)n);
};


void ExchangeParticleData::calculateHydrodynamicInteractions(Lattice2D *lattice2D_, Unit_Conversion *unitConversion, int numberParticles, LAMMPS_NS::tagint *tag, double **xPart, double *rp, int hydroForceID, int hydroTorqueID, int stressletID, LAMMPS_NS::Atom *atom, vector<double> origin)
{
  for(int iPart = 0; iPart < numberParticles; ++iPart){
      double x_lb = fmin(fmax(unitConversion->get_pos_lb(xPart[iPart][0]-(lattice2D_->getProcOrigin()[0]-(double)lattice2D_->get_envelopeWidth()*unitConversion->get_dx())), 0.0), (double)lattice2D_->get_nx());
      double y_lb = fmin(fmax(unitConversion->get_pos_lb(xPart[iPart][1]-(lattice2D_->getProcOrigin()[1]-(double)lattice2D_->get_envelopeWidth()*unitConversion->get_dx())), 0.0), (double)lattice2D_->get_ny());
      double r_lb = unitConversion->get_radius_lb(*rp);

      int nodeZone[2][2] = {{(int)(x_lb-r_lb)-3, (int)(x_lb+r_lb)+3}, {(int)(y_lb-r_lb)-3, (int)(y_lb+r_lb)+3}};
      for (int i = 0; i < 2; i++){
          nodeZone[0][i] = fmin(fmax(nodeZone[0][i],0),lattice2D_->get_nx());
          nodeZone[1][i] = fmin(fmax(nodeZone[1][i],0),lattice2D_->get_nx());
      }
      
      for(int i = nodeZone[0][0]; i < nodeZone[0][1]; ++i){
        for(int j = nodeZone[1][0]; j < nodeZone[1][1]; ++j){

          int ind_phys_1D = i * lattice2D_->get_ny() + j;
          
          // Todo extend to two or more particles
          LAMMPS_NS::tagint pID = lattice2D_->getParticleDataOnLatticeNode(ind_phys_1D).particleID[0];
          vector<double> Fhyd = lattice2D_->getParticleDataOnLatticeNode(ind_phys_1D).hydrodynamicForce;

          atom->dvector[hydroForceID][0] += Fhyd[0]*unitConversion->get_forceFactor();
          atom->dvector[hydroForceID][1] += Fhyd[1]*unitConversion->get_forceFactor();
  //        atom->dvector[hydroForceID][2] += Fhyd[2];
    
  // TODO torque and stresslet      
  //        atom->dvector[hydroTorqueID][0] += Fhyd[0]*dx etc;
  //        atom->dvector[hydroTorqueID][1] += Fhyd[1];
  //        atom->dvector[hydroTorqueID][2] += Fhyd[2];

        }
      }
  }
}
