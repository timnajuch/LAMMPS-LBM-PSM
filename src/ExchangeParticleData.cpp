/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2019
------------------------------------------------------*/

#include "ExchangeParticleData.h"
#include "domain.h"

ExchangeParticleData::ExchangeParticleData(int dimension_){
  dimension = dimension_;
};

ExchangeParticleData::~ExchangeParticleData() {};


void ExchangeParticleData::setParticlesOnLattice(Lattice2D *lattice2D_, Unit_Conversion *unitConversion, int numberParticles, LAMMPS_NS::tagint *tag, double **xPart, double **uPart, double **omega, double *rp, vector<double> boxLength, vector<double> origin)
{
  for(int iPart = 0; iPart < numberParticles; ++iPart){
      double x_lb = fmin(fmax(unitConversion->get_pos_lb(xPart[iPart][0]-(lattice2D_->getProcOrigin()[0]-(double)lattice2D_->get_envelopeWidth()*unitConversion->get_dx())), 0.0), (double)lattice2D_->get_nx()-1.0);
      double y_lb = fmin(fmax(unitConversion->get_pos_lb(xPart[iPart][1]-(lattice2D_->getProcOrigin()[1]-(double)lattice2D_->get_envelopeWidth()*unitConversion->get_dx())), 0.0), (double)lattice2D_->get_ny()-1.0);
      double z_lb = fmin(fmax(unitConversion->get_pos_lb(xPart[iPart][2]-(lattice2D_->getProcOrigin()[2]-(double)lattice2D_->get_envelopeWidth()*unitConversion->get_dx())), 0.0), (double)lattice2D_->get_nz()-1.0);
      double r_lb = unitConversion->get_radius_lb(rp[iPart]);

      int nodeZone[3][2] = {{(int)(x_lb-r_lb)-3, (int)(x_lb+r_lb)+3}, 
                            {(int)(y_lb-r_lb)-3, (int)(y_lb+r_lb)+3},
                            {(int)(z_lb-r_lb)-3, (int)(z_lb+r_lb)+3}};
      for (int i = 0; i < 2; i++){
          nodeZone[0][i] = fmin(fmax(nodeZone[0][i],0),lattice2D_->get_nx());
          nodeZone[1][i] = fmin(fmax(nodeZone[1][i],0),lattice2D_->get_ny());
          nodeZone[2][i] = fmin(fmax(nodeZone[2][i],0),lattice2D_->get_nz());
      }
      
      for(int i = nodeZone[0][0]; i < nodeZone[0][1]; ++i){
        for(int j = nodeZone[1][0]; j < nodeZone[1][1]; ++j){
          for(int k = nodeZone[2][0]; k < nodeZone[2][1]; ++k){

            int ind_phys_1D = i * lattice2D_->get_ny() * lattice2D_->get_nz() + j * lattice2D_->get_nz() + k;

            double sf = calcSolidFraction(i, j, k, x_lb, y_lb, z_lb, r_lb);

            int id_old = 0;
            double sf_old = 0.0;
            // Todo: Extend to two or more particles
            if (tag[iPart] == lattice2D_->getParticleDataOnLatticeNode(ind_phys_1D).particleID[0]){
              id_old = lattice2D_->getParticleDataOnLatticeNode(ind_phys_1D).particleID[0];
              sf_old = lattice2D_->getParticleDataOnLatticeNode(ind_phys_1D).solidFraction[0];
            }

            int const decFlag = (sf > 0.001) + 2*(sf_old > 0.001);
            double zeroVel[3] = {0.0, 0.0, 0.0};

            double dx = (double)i - x_lb;
            double dy = (double)j - y_lb;
//            double dz = 0.0;
            double dz = (double)k - z_lb;  
            double uNode[3];
            uNode[0] = unitConversion->get_vel_lb(uPart[iPart][0]);
            uNode[1] = unitConversion->get_vel_lb(uPart[iPart][1]);
            uNode[2] = unitConversion->get_vel_lb(uPart[iPart][3]);
            uNode[0] += unitConversion->get_freq_lb(omega[iPart][1])*dz - unitConversion->get_freq_lb(omega[iPart][2])*dy;
            uNode[1] += -unitConversion->get_freq_lb(omega[iPart][0])*dz + unitConversion->get_freq_lb(omega[iPart][2])*dx;
            uNode[2] += unitConversion->get_freq_lb(omega[iPart][0])*dy - unitConversion->get_freq_lb(omega[iPart][1])*dx;
            //uNode[2] += omega[0]*dy - omega[1]*dx;

            switch(decFlag){
              case 0: // sf == 0 && sf_old == 0
                lattice2D_->setParticleOnLattice(ind_phys_1D, 0, zeroVel, 0.0);
                break;
              case 1: // sf > 0 && sf_old == 0
                lattice2D_->setParticleOnLattice(ind_phys_1D, tag[iPart], uNode, calcSolidFraction(i, j, k, x_lb, y_lb, z_lb, r_lb));
                break;
              case 2: // sf == 0 && sf_old > 0
                if( id_old == tag[iPart] ) // then particle has left this cell
                  lattice2D_->setParticleOnLattice(ind_phys_1D, 0, zeroVel, 0.0);
                break;
              case 3: // sf > 0 && sf_old > 0
                if( sf > sf_old || id_old == tag[iPart] )
                  lattice2D_->setParticleOnLattice(ind_phys_1D, tag[iPart], uNode, calcSolidFraction(i, j, k, x_lb, y_lb, z_lb, r_lb));
                break;
            }

          }
        }
      }
  }
};



double ExchangeParticleData::calcSolidFraction(int i, int j, int k, double xP_LB, double yP_LB, double zP_LB, double rP_LB){

    int slicesPerDim = 10;
    double sliceWidth = 1.0/((double)slicesPerDim);
    double fraction = 1.0/((double)(slicesPerDim*slicesPerDim));
    //if(domain->dimension == 3)
    if (dimension == 3)
      fraction = 1.0/((double)(slicesPerDim*slicesPerDim*slicesPerDim));

    double sqrt2half = sqrt(2.0)/2.0;
    //if(domain->dimension == 3){
    if (dimension == 3){
      sqrt2half = sqrt(3.0)/2.0;
    }
    double dx = (double)i - xP_LB;
    double dy = (double)j - yP_LB;
    double dz = (double)k - zP_LB;

    double dist = dx*dx + dy*dy + dz*dz;

    double rP = rP_LB + sqrt2half;
    if(dist > rP*rP) return 0.0;

    double rM = rP_LB - sqrt2half;
    if(dist < rM*rM) return  1.0;

    double rSq = rP_LB*rP_LB;
    double dx_sq[slicesPerDim], dy_sq[slicesPerDim], dz_sq[slicesPerDim];

    for(int i = 0; i < slicesPerDim; ++i){
      double delta = -0.5 + ((double)i+0.5)*sliceWidth;
      double dx_i = dx+delta; dx_sq[i] = dx_i*dx_i;
      double dy_i = dy+delta; dy_sq[i] = dy_i*dy_i;
      double dz_i = dz+delta; dz_sq[i] = dz_i*dz_i;
    }   

    int n = 0;
    //if(domain->dimension == 2){
    if (dimension == 2){
      for(int i = 0; i < slicesPerDim; ++i){
        for(int j = 0; j < slicesPerDim; ++j){
          if(dx_sq[i] + dy_sq[j] < rSq) ++n;
        }   
      }   
    }else{
      for(int i = 0; i < slicesPerDim; ++i){
        for(int j = 0; j < slicesPerDim; ++j){
          for(int k = 0; k < slicesPerDim; ++k){
            if(dx_sq[i] + dy_sq[j] + dz_sq[k] < rSq) ++n;
          }
        }
      }
    }

    return fraction*((double)n);
};


void ExchangeParticleData::calculateHydrodynamicInteractions(Lattice2D *lattice2D_, Unit_Conversion *unitConversion, double *xPart, double rp, vector<double> &fHydro, vector<double> &tHydro, vector<double> &stresslet)
{
    int envelopeWidth = lattice2D_->get_envelopeWidth();

    double x_lb = fmin(fmax(unitConversion->get_pos_lb(xPart[0]-(lattice2D_->getProcOrigin()[0]-(double)envelopeWidth*unitConversion->get_dx())), 0.0), (double)lattice2D_->get_nx()-1.0);
    double y_lb = fmin(fmax(unitConversion->get_pos_lb(xPart[1]-(lattice2D_->getProcOrigin()[1]-(double)envelopeWidth*unitConversion->get_dx())), 0.0), (double)lattice2D_->get_ny()-1.0);
    double z_lb = fmin(fmax(unitConversion->get_pos_lb(xPart[2]-(lattice2D_->getProcOrigin()[2]-(double)envelopeWidth*unitConversion->get_dx())), 0.0), (double)lattice2D_->get_nz()-1.0);
    double r_lb = unitConversion->get_radius_lb(rp);

    int nodeZone[3][2] = {{(int)(x_lb-r_lb)-3, (int)(x_lb+r_lb)+3},
                          {(int)(y_lb-r_lb)-3, (int)(y_lb+r_lb)+3},
                          {(int)(z_lb-r_lb)-3, (int)(z_lb+r_lb)+3}};
    for (int i = 0; i < 2; i++){
        nodeZone[0][i] = fmin(fmax(nodeZone[0][i],envelopeWidth),lattice2D_->get_nx()-envelopeWidth);
        nodeZone[1][i] = fmin(fmax(nodeZone[1][i],envelopeWidth),lattice2D_->get_ny()-envelopeWidth);
        nodeZone[2][i] = fmin(fmax(nodeZone[2][i],envelopeWidth),lattice2D_->get_nz()-envelopeWidth);
    }
    
    for(int i = nodeZone[0][0]; i < nodeZone[0][1]; ++i){
      for(int j = nodeZone[1][0]; j < nodeZone[1][1]; ++j){
        for(int k = nodeZone[2][0]; k < nodeZone[2][1]; ++k){

          int ind_phys_1D = i * lattice2D_->get_ny() * lattice2D_->get_nz() + j * lattice2D_->get_nz() + k;
          
          // Todo extend to two or more particles
          LAMMPS_NS::tagint pID = lattice2D_->getParticleDataOnLatticeNode(ind_phys_1D).particleID[0];
          vector<double> Fhyd = lattice2D_->getParticleDataOnLatticeNode(ind_phys_1D).hydrodynamicForce;

          fHydro[0] += Fhyd[0]*unitConversion->get_forceFactor();
          fHydro[1] += Fhyd[1]*unitConversion->get_forceFactor();
          fHydro[2] += Fhyd[2]*unitConversion->get_forceFactor();

          double dx = (double)i - x_lb;
          double dy = (double)j - y_lb;
          double dz = (double)k - z_lb;  

          tHydro[0] += (dy*Fhyd[2] - dz*Fhyd[1])*unitConversion->get_torqueFactor();
          tHydro[1] += (dz*Fhyd[0] - dx*Fhyd[2])*unitConversion->get_torqueFactor();
          tHydro[2] += (dx*Fhyd[1] - dy*Fhyd[0])*unitConversion->get_torqueFactor();

          stresslet[0] -= 0.5*(dx*Fhyd[0] + dx*Fhyd[0])*unitConversion->get_torqueFactor();
          stresslet[1] -= 0.5*(dy*Fhyd[1] + dy*Fhyd[1])*unitConversion->get_torqueFactor();
          stresslet[2] -= 0.5*(dz*Fhyd[2] + dz*Fhyd[2])*unitConversion->get_torqueFactor();
          stresslet[3] -= 0.5*(dy*Fhyd[0] + dx*Fhyd[1])*unitConversion->get_torqueFactor();
          stresslet[4] -= 0.5*(dz*Fhyd[0] + dx*Fhyd[2])*unitConversion->get_torqueFactor();
          stresslet[5] -= 0.5*(dz*Fhyd[1] + dy*Fhyd[2])*unitConversion->get_torqueFactor();

        }
      }
    }
};
