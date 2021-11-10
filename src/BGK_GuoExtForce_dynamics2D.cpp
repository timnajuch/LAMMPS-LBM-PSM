/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#include "BGK_GuoExtForce_dynamics2D.h"

BGK_GuoExtForce_Dynamics2D::BGK_GuoExtForce_Dynamics2D(double tau_, int nx_, int ny_, int nz_, int q_, vector<double> F_lbm_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_) :
  Dynamics2D(nx_, ny_, nz_, q_, decomposition_, procCoordinates_, origin_, boxLength_, dimension_), tau(tau_), F_lbm(F_lbm_) {};


BGK_GuoExtForce_Dynamics2D::~BGK_GuoExtForce_Dynamics2D(){};


void BGK_GuoExtForce_Dynamics2D::compute_macro_values(){
  for(int i = 0; i < Lattice2D::nx; ++i){
    for(int j = 0; j < Lattice2D::ny; ++j){
      for(int k = 0; k < Lattice2D::nz; ++k){

        int ind_phys_1D = i * Lattice2D::ny * Lattice2D::nz + j * Lattice2D::nz + k;
        int ind_phys_2D = i * Lattice2D::ny * Lattice2D::nz * 3 + j*Lattice2D::nz*3 + k*3;

        double rho_tmp = 0.0;
        double jx = 0.0;
        double jy = 0.0;
        double jz = 0.0;
      
        for(int iq = 0; iq < Lattice2D::q; ++iq){
          int ind_iq = i * Lattice2D::ny * Lattice2D::nz * Lattice2D::q + j * Lattice2D::nz*Lattice2D::q + k * Lattice2D::q + iq; 
          rho_tmp += Lattice2D::get_f(ind_iq);
          jx += Lattice2D::get_f(ind_iq) * Lattice2D::e[3*iq];
          jy += Lattice2D::get_f(ind_iq) * Lattice2D::e[3*iq+1];
          jz += Lattice2D::get_f(ind_iq) * Lattice2D::e[3*iq+2];
        }
        Lattice2D::rho[ind_phys_1D] = rho_tmp;

  //      jx += F_lbm[0]/2.0;
  //      jy += F_lbm[1]/2.0;
        Lattice2D::u[ind_phys_2D] = jx/rho_tmp;
        Lattice2D::u[ind_phys_2D+1] = jy/rho_tmp;
        Lattice2D::u[ind_phys_2D+2] = jz/rho_tmp;

        //p[i][j] = rho_tmp*csPow2;
      }
    }
  }
};


void BGK_GuoExtForce_Dynamics2D::compute_macro_values(int i_, int j_, int k_){
  int ind_phys_1D = i_ * Lattice2D::ny + j_ * Lattice2D::nz + k_;
  int ind_phys_2D = i_ * Lattice2D::ny * Lattice2D::nz * 3 + j_*Lattice2D::nz*3 + k_*3;

  double rho_tmp = 0.0;
  double jx = 0.0;
  double jy = 0.0;
  double jz = 0.0;

  for(int iq = 0; iq < Lattice2D::q; ++iq){
    int ind_iq = i_ * Lattice2D::ny * Lattice2D::nz * Lattice2D::q + j_ * Lattice2D::nz * Lattice2D::q + k_*Lattice2D::q + iq; 
    rho_tmp += Lattice2D::get_f(ind_iq);
    jx += Lattice2D::get_f(ind_iq) * Lattice2D::e[3*iq];
    jy += Lattice2D::get_f(ind_iq) * Lattice2D::e[3*iq+1];
    jz += Lattice2D::get_f(ind_iq) * Lattice2D::e[3*iq+2];
  }
  Lattice2D::rho[ind_phys_1D] = rho_tmp;

//      jx += F_lbm[0]/2.0;
//      jy += F_lbm[1]/2.0;
  Lattice2D::u[ind_phys_2D] = jx/rho_tmp;
  Lattice2D::u[ind_phys_2D+1] = jy/rho_tmp;
  Lattice2D::u[ind_phys_2D+2] = jz/rho_tmp;

  //p[i][j] = rho_tmp*csPow2;
};


void BGK_GuoExtForce_Dynamics2D::collision(){
  for(int i = 0; i < Lattice2D::nx; ++i){
    for(int j = 0; j < Lattice2D::ny; ++j){
      for(int k = 0; k < Lattice2D::nz; ++k){
        for(int iq = 0; iq < Lattice2D::q; ++iq){

          int ind_iq = i * Lattice2D::ny * Lattice2D::nz * Lattice2D::q + j * Lattice2D::nz * Lattice2D::q + k*Lattice2D::q + iq;
          int ind_phys_1D = i * Lattice2D::ny * Lattice2D::nz + j * Lattice2D::nz + k;
          int ind_phys_2D = i * Lattice2D::ny * Lattice2D::nz * 3 + j*Lattice2D::nz*3 + k*3;
   
          Lattice2D::set_f0(i, j, k, iq, Dynamics2D::feq(iq, ind_phys_1D, ind_phys_2D, Lattice2D::rho, Lattice2D::u) );
          vector<double> uSolid{Lattice2D::getSolidVelocityOnLattice(ind_phys_1D, 0)[0], Lattice2D::getSolidVelocityOnLattice(ind_phys_1D, 0)[1], Lattice2D::getSolidVelocityOnLattice(ind_phys_1D, 0)[2]};
          double f0_solid = Dynamics2D::feq(iq, ind_phys_1D, ind_phys_2D, Lattice2D::get_rho(ind_phys_1D), uSolid);

          // TODO iq0 etc
          int iqO = 0;    
          if ((iq == 1) || (iq == 5) || (iq == 2) || (iq == 6)) 
            iqO = iq + 2;
          if ((iq == 3) || (iq == 7) || (iq == 4) || (iq == 8)) 
            iqO = iq - 2;

          int ind_iqO = i * Lattice2D::ny * Lattice2D::nz * Lattice2D::q + j * Lattice2D::nz * Lattice2D::q + k * Lattice2D::q + iqO;
          
          double f0iqO = Dynamics2D::feq(iqO, ind_phys_1D, ind_phys_2D, Lattice2D::rho, Lattice2D::u);
          double f0iqO_solid = Dynamics2D::feq(iqO, ind_phys_1D, ind_phys_2D, Lattice2D::rho, uSolid);


          double BGKcoll = ( Lattice2D::get_f0(ind_iq) - Lattice2D::get_f(ind_iq) ) / tau;

          double solid_coll = f0_solid - Lattice2D::get_f(ind_iq) + ( 1.0 - 1.0/tau) * (Lattice2D::get_f(ind_iq) - Lattice2D::get_f0(ind_iq) );
  //        double solid_coll = Lattice2D::get_f(ind_iqO) - Lattice2D::get_f(ind_iq) + f0_solid - f0iqO;

          double B = Lattice2D::getSolidFractionOnLattice(ind_phys_1D, 0);
          
          Lattice2D::set_fcoll(i, j, k, iq, Lattice2D::get_f(ind_iq)  + ( 1.0 - B ) * BGKcoll  + B * solid_coll 
                  + ( 1.0 - B ) * (Lattice2D::g[iq] / Lattice2D::c) * ( Lattice2D::e[3*iq] * F_lbm[0] + Lattice2D::e[3*iq+1] * F_lbm[1] + Lattice2D::e[3*iq+2] * F_lbm[2] ) ); //TODO forcing term not correct. not g. just use w divided by cs**2

          if(iq == 0)
          {
            Lattice2D::set_Fhydx(ind_phys_1D, B * solid_coll * Lattice2D::e[3*iq]);
            Lattice2D::set_Fhydy(ind_phys_1D, B * solid_coll * Lattice2D::e[3*iq+1]);
            Lattice2D::set_Fhydz(ind_phys_1D, B * solid_coll * Lattice2D::e[3*iq+2]);
          }else{
            Lattice2D::add_Fhydx(ind_phys_1D, B * solid_coll * Lattice2D::e[3*iq]);
            Lattice2D::add_Fhydy(ind_phys_1D, B * solid_coll * Lattice2D::e[3*iq+1]);
            Lattice2D::add_Fhydz(ind_phys_1D, B * solid_coll * Lattice2D::e[3*iq+2]);
          }

  /*      
          Lattice2D::set_fcoll(i, j, iq, Lattice2D::get_f(ind_iq) + ( 1.0 - B ) * BGKcoll + B * solid_coll 
                   + (1.0-1.0/(2.0*tau)) * Lattice2D::w[iq] * Lattice2D::rho[ind_phys_1D] * (
                     ( (Lattice2D::e[2*iq] - Lattice2D::u[ind_phys_2D]) / (Lattice2D::cs * Lattice2D::cs) 
                          + (Lattice2D::e[2*iq] * Lattice2D::u[ind_phys_2D] + Lattice2D::e[2*iq+1] * Lattice2D::u[ind_phys_2D+1]) / (pow(Lattice2D::cs,4.0)) 
                                * Lattice2D::e[2*iq] ) * F_lbm[0]
                   + (  (Lattice2D::e[2*iq+1] - Lattice2D::u[ind_phys_2D+1]) / (Lattice2D::cs * Lattice2D::cs) 
                          + (Lattice2D::e[2*iq] * Lattice2D::u[ind_phys_2D] + Lattice2D::e[2*iq+1] * Lattice2D::u[ind_phys_2D+1]) / (pow(Lattice2D::cs,4.0)) 
                                * Lattice2D::e[2*iq+1] ) * F_lbm[1]) ); 
  */
        }
      }
    }
  }

};


void BGK_GuoExtForce_Dynamics2D::collision(int i_, int j_, int k_, int iq_){
        
  int ind_iq = i_ * Lattice2D::ny * Lattice2D::nz * Lattice2D::q + j_ * Lattice2D::nz * Lattice2D::q + k_*Lattice2D::q + iq_;
  int ind_phys_1D = i_ * Lattice2D::ny *Lattice2D::nz + j_*Lattice2D::nz + k_;
  int ind_phys_2D = i_ * Lattice2D::ny *Lattice2D::nz* 3 + j_*3*Lattice2D::nz + k_*3;

  Lattice2D::set_f0(i_, j_, k_, iq_, Dynamics2D::feq(iq_, ind_phys_1D, ind_phys_2D, Lattice2D::rho, Lattice2D::u) );
  vector<double> uSolid{Lattice2D::getSolidVelocityOnLattice(ind_phys_1D, 0)[0], Lattice2D::getSolidVelocityOnLattice(ind_phys_1D, 0)[1], Lattice2D::getSolidVelocityOnLattice(ind_phys_1D, 0)[2]};
  double f0_solid = Dynamics2D::feq(iq_, ind_phys_1D, ind_phys_2D, Lattice2D::get_rho(ind_phys_1D), uSolid);


  double BGKcoll = ( Lattice2D::get_f0(ind_iq) - Lattice2D::get_f(ind_iq) ) / tau;

  // TODO: iq_m calc etc. Shoyld be *at leadt for 3D, +1 or -1 for odd and even respectively. Rewrite perhaps also 2D for consistency reasons
  int iq_m = 0;
  if(iq_ > 0)
  {
    if(iq_ < 5)
    {
      iq_m = iq_ < 3 ? iq_ + 2 : iq_ - 2;
    }
    else
    {
      iq_m = iq_ < 7 ? iq_ + 2 : iq_ - 2;
    }
  }
  int ind_iq_m = i_ * Lattice2D::ny * Lattice2D::nz * Lattice2D::q + j_ * Lattice2D::nz * Lattice2D::q + k_*Lattice2D::q + iq_m;
  Lattice2D::set_f0(i_, j_, k_, iq_m, Dynamics2D::feq(iq_m, ind_phys_1D, ind_phys_2D, Lattice2D::rho, Lattice2D::u) );

  double solid_coll = f0_solid - Lattice2D::get_f(ind_iq) + ( 1.0 - 1.0/tau) * (Lattice2D::get_f(ind_iq) - Lattice2D::get_f0(ind_iq) );
  //double solid_coll = f0_solid - Lattice2D::get_f(ind_iq) + Lattice2D::get_f(ind_iq_m) - Lattice2D::get_f0(ind_iq_m);

  double B = Lattice2D::getSolidFractionOnLattice(ind_phys_1D, 0);
  //double B = Lattice2D::get_B(ind_phys_1D);
        
  Lattice2D::set_f(i_, j_, k_, iq_, Lattice2D::get_f(ind_iq)  + ( 1.0 - B ) * BGKcoll  + B * solid_coll 
          + ( 1.0 - B ) * (Lattice2D::g[iq_] / Lattice2D::c) * ( Lattice2D::e[3*iq_] * F_lbm[0] + Lattice2D::e[3*iq_+1] * F_lbm[1] + Lattice2D::e[3*iq_+2] * F_lbm[2]) ); //TODO forcing term. See above

  // todo extend to two or more particles
  if(iq_ == 0)
  {
    Lattice2D::set_Fhydx(ind_phys_1D, -B * solid_coll * Lattice2D::e[3*iq_]);
    Lattice2D::set_Fhydy(ind_phys_1D, -B * solid_coll * Lattice2D::e[3*iq_+1]);
    Lattice2D::set_Fhydz(ind_phys_1D, -B * solid_coll * Lattice2D::e[3*iq_+2]);
  }else{
    Lattice2D::add_Fhydx(ind_phys_1D, -B * solid_coll * Lattice2D::e[3*iq_]);
    Lattice2D::add_Fhydy(ind_phys_1D, -B * solid_coll * Lattice2D::e[3*iq_+1]);
    Lattice2D::add_Fhydz(ind_phys_1D, -B * solid_coll * Lattice2D::e[3*iq_+2]);
  }

};


void  BGK_GuoExtForce_Dynamics2D::initialise_dynamics(double rho_, double ux_, double uy_, double uz_){

  for(int i = 0; i < Lattice2D::nx; ++i){
    for(int j = 0; j < Lattice2D::ny; ++j){
      for(int k = 0; k < Lattice2D::nz; ++k){
        int ind_phys_1D = i * Lattice2D::ny * Lattice2D::nz + j*Lattice2D::nz + k;
        int ind_phys_2D = i * Lattice2D::ny *Lattice2D::nz * 3 + j*Lattice2D::nz*3 + k*3;
        Lattice2D::rho[ind_phys_1D] = rho_;
        Lattice2D::u[ind_phys_2D] = ux_;
        Lattice2D::u[ind_phys_2D+1] = uy_;
        Lattice2D::u[ind_phys_2D+3] = uz_;

        for(int iq = 0; iq < Lattice2D::q; ++iq){
            int ind_iq = i * Lattice2D::ny * Lattice2D::nz * Lattice2D::q + j * Lattice2D::nz *Lattice2D::q + k*Lattice2D::q + iq; 

            Lattice2D::set_f0(i, j, k, iq, Dynamics2D::feq(iq, ind_phys_1D, ind_phys_2D, Lattice2D::rho, Lattice2D::u) );
            Lattice2D::set_f(i, j, k, iq, Dynamics2D::feq(iq, ind_phys_1D, ind_phys_2D, Lattice2D::rho, Lattice2D::u) );
            Lattice2D::set_fcoll(i, j, k, iq, Dynamics2D::feq(iq, ind_phys_1D, ind_phys_2D, Lattice2D::rho, Lattice2D::u) );

        }
      }
    }
  }

};

void BGK_GuoExtForce_Dynamics2D::macroCollideStream(){
// TODO extend to 3D corenr, edges etc
  for(int i = 0; i < Lattice2D::nx; ++i){
    for(int j = 0; j < Lattice2D::ny; ++j){
      for(int k = 0; k < Lattice2D::nz; ++k){

        int corner = 0;

        //if(domain->dimension == 2){
        if(Lattice2D::dimension == 2){
          if(i == 0                && j == 0)                 corner = 1;
          if(i == 0                && j == Lattice2D::ny-1)   corner = 2;
          if(i == Lattice2D::nx-1  && j == 0)                 corner = 1;
          if(i == Lattice2D::nx-1  && j == Lattice2D::ny-1)   corner = 2;
        }
/*else{
          // Edges x-y
          if(i == 0                && j == 0                && k > 0  && k < Lattice2D::nz-1)
            corner = 1;
          if(i == 0                && j == Lattice2D::ny-1  && k > 0  && k < Lattice2D::nz-1)
            corner = 2;
          if(i == Lattice2D::nx-1  && j == 0                && k > 0  && k < Lattice2D::nz-1)
            corner = 1;
          if(i == Lattice2D::nx-1  && j == Lattice2D::ny-1  && k > 0  && k < Lattice2D::nz-1)
            corner = 2;
          // Edges y-z
          if(k == 0                && j == 0                && i > 0  && i < Lattice2D::nx-1)
            corner = 4;
          if(k == 0                && j == Lattice2D::ny-1  && i > 0  && i < Lattice2D::nx-1)
            corner = 5;
          if(k == Lattice2D::nx-1  && j == 0                && i > 0  && i < Lattice2D::nx-1)
            corner = 6;
          if(k == Lattice2D::nx-1  && j == Lattice2D::ny-1  && i > 0  && i < Lattice2D::nx-1)
            corner = 7;

          if(i == 0                && j == 0                && k == 0)
            corner = 3;
          if(i == 0                && j == 0                && k == Lattice2D::nz-1)
            corner = 4;
          if(i == 0                && j == Lattice2D::ny-1  && k == 0)
            corner = 4;
          if(i == 0                && j == Lattice2D::ny-1  && k == Lattice2D::nz-1)
            corner = 5;
          if(i == Lattice2D::nx-1  && j == 0                && k == 0)
            corner = 6;
          if(i == Lattice2D::nx-1  && j == 0                && k == Lattice2D::nz-1)
            corner = 7;
          if(i == Lattice2D::nx-1  && j == Lattice2D::ny-1  && k == 0)
            corner = 8;
          if(i == Lattice2D::nx-1  && j == Lattice2D::ny-1  && k == Lattice2D::nz-1)
            corner = 9;
        }
*/
        compute_macro_values(i, j, k);
        for(int iq = 0; iq < Lattice2D::q; ++iq){
          
          collision(i, j, k, iq);
          
          //if(domain->dimension == 2){
          if(Lattice2D::dimension == 2){
            if( (i > 0) && (i < Lattice2D::nx-1) && (j > 0) && (j < Lattice2D::ny-1) )
              Dynamics2D::streamBulk(i, j, k, iq);
            if(i == 0)                  Dynamics2D::streamBC_xn(i,j,iq, corner);
            if(i == Lattice2D::nx - 1)  Dynamics2D::streamBC_xp(i,j,iq, corner);
            if(j == 0)                  Dynamics2D::streamBC_yn(i,j,iq, corner);
            if(j == Lattice2D::ny - 1)  Dynamics2D::streamBC_yp(i,j,iq, corner);
          }else{
            if( (i > 0) && (i < Lattice2D::nx-1) && (j > 0) && (j < Lattice2D::ny-1) && (k > 0) && (k < Lattice2D::nz-1) )
              Dynamics2D::streamBulk(i, j, k, iq);
            if(i == 0)                  Dynamics2D::streamBC(i,j,k,iq);
            if(i == Lattice2D::nx - 1)  Dynamics2D::streamBC(i,j,k,iq);
            if(j == 0)                  Dynamics2D::streamBC(i,j,k,iq);
            if(j == Lattice2D::ny - 1)  Dynamics2D::streamBC(i,j,k,iq);
            if(k == 0)                  Dynamics2D::streamBC(i,j,k,iq);
            if(k == Lattice2D::nz - 1)  Dynamics2D::streamBC(i,j,k,iq);
            //TODO combine bulk and BC streaming?
          }
            

//          if(domain->dimension == 3){
//            if(k == 0)                  Dynamics2D::streamBC_zn(i,j,k,iq, corner);
//            if(k == Lattice2D::nz - 1)  Dynamics2D::streamBC_zp(i,j,k,iq, corner);
//          }
        }

      }
    }
  }

  Lattice2D::cp_fcoll_f();

};
