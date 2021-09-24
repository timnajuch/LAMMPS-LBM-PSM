/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#include "dynamics2D.h"

Dynamics2D::Dynamics2D(int nx_, int ny_, int q_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_) :
  Lattice2D(nx_, ny_, q_, decomposition_, procCoordinates_, origin_, boxLength_)
{
  for(int i = 0; i < nx; ++i){
    for(int j = 0; j < ny; ++j){
      for(int iq = 0; iq < q; ++iq){
        int ind_phys_1D = i * Lattice2D::ny + j;
        int ind_phys_2D = i * Lattice2D::ny * 2 + j*2;
        int ind_iq = i * Lattice2D::ny * Lattice2D::q + j * Lattice2D::q + iq;

        Lattice2D::set_f0(i, j, iq, Dynamics2D::feq(iq, ind_phys_1D, ind_phys_2D, Lattice2D::rho, Lattice2D::u) );
        Lattice2D::set_f(i, j, iq, Dynamics2D::feq(iq, ind_phys_1D, ind_phys_2D, Lattice2D::rho, Lattice2D::u) );
        Lattice2D::set_fcoll(i, j, iq, Dynamics2D::feq(iq, ind_phys_1D, ind_phys_2D, Lattice2D::rho, Lattice2D::u) );
      }
    }
  }
};


Dynamics2D::~Dynamics2D(){};
    

double Dynamics2D::feq(int iq_, int ind_phys_1D_, int ind_phys_2D_, vector<double> &rho_, vector<double> &u_){
  return rho_[ind_phys_1D_] * Lattice2D::w[iq_] * 
        (1.0 + ( Lattice2D::e[2*iq_] * u_[ind_phys_2D_] + Lattice2D::e[2*iq_+1] * u_[ind_phys_2D_+1] ) / (pow(Lattice2D::cs, 2.0))
        + pow( Lattice2D::e[2*iq_] * u_[ind_phys_2D_] + Lattice2D::e[iq_*2+1] * u_[ind_phys_2D_+1] , 2.0) / (2.0*pow(Lattice2D::cs, 4.0))
        - 1.0/2.0 * ( u_[ind_phys_2D_] * u_[ind_phys_2D_] + u_[ind_phys_2D_+1] * u_[ind_phys_2D_+1] ) / (pow(Lattice2D::cs,2.0) ) );
}


double Dynamics2D::feq(int iq_, int ind_phys_1D_, int ind_phys_2D_, double rho, vector<double> u){
  return rho * Lattice2D::w[iq_] * 
        (1.0 + ( Lattice2D::e[2*iq_] * u[0] + Lattice2D::e[2*iq_+1] * u[1] ) / (pow(Lattice2D::cs, 2.0))
        + pow( Lattice2D::e[2*iq_] * u[0] + Lattice2D::e[iq_*2+1] * u[1] , 2.0) / (2.0*pow(Lattice2D::cs, 4.0))
        - 1.0/2.0 * ( u[0] * u[0] + u[1] * u[1] ) / (pow(Lattice2D::cs,2.0) ) );
}


void Dynamics2D::streaming_periodic_2D(){
  int in, ip, jn, jp;
  in = 0;
  ip = 0;
  jn = 0;
  jp = 0;
  if(Lattice2D::q == 9){

    for(int i = 0; i < Lattice2D::nx; ++i){
      in = (i>0   )               ? (i-1):(Lattice2D::nx-1);
      ip = (i<Lattice2D::nx-1) ? (i+1):(0   );

      for(int j = 0; j < Lattice2D::ny; ++j){
        jn = (j>0   )               ? (j-1):(Lattice2D::ny-1);
        jp = (j<Lattice2D::ny-1) ? (j+1):(0   );

        Lattice2D::set_f( i ,j ,0, Lattice2D::get_fcoll( i,j,0 ) );
        Lattice2D::set_f( ip,j ,1, Lattice2D::get_fcoll( i,j,1 ) );
        Lattice2D::set_f( i ,jp,2, Lattice2D::get_fcoll( i,j,2 ) );
        Lattice2D::set_f( in,j ,3, Lattice2D::get_fcoll( i,j,3 ) );
        Lattice2D::set_f( i ,jn,4, Lattice2D::get_fcoll( i,j,4 ) );
        Lattice2D::set_f( ip,jp,5, Lattice2D::get_fcoll( i,j,5 ) );
        Lattice2D::set_f( in,jp,6, Lattice2D::get_fcoll( i,j,6 ) );
        Lattice2D::set_f( in,jn,7, Lattice2D::get_fcoll( i,j,7 ) );
        Lattice2D::set_f( ip,jn,8, Lattice2D::get_fcoll( i,j,8 ) );
      }
    }
  }else{
    std::cout << "Only D2Q9 model (explicitly) implemented. Other models require a more general streaming loop!" << endl;
  }
};


void Dynamics2D::streaming(){
  int in, ip, jn, jp;
  in = 0;
  ip = 0;
  jn = 0;
  jp = 0;
  if(Lattice2D::q == 9){

    int nx1 = Lattice2D::nx;
    int ny1 = Lattice2D::ny;

    for(int i = 0; i < nx1; ++i){
      in = (i>0   ) ?(i-1):(0     );
      ip = (i<nx1-1)?(i+1):(nx1-1 );

      for(int j = 0; j < ny1; ++j){
        jn = (j>0    )?(j-1):(ny1-1);
        jp = (j<ny1-1)?(j+1):(0    );

        Lattice2D::set_f( i ,j ,0, Lattice2D::get_fcoll( i,j,0 ) );
        Lattice2D::set_f( ip,j ,1, Lattice2D::get_fcoll( i,j,1 ) );
        Lattice2D::set_f( i ,jp,2, Lattice2D::get_fcoll( i,j,2 ) );
        Lattice2D::set_f( in,j ,3, Lattice2D::get_fcoll( i,j,3 ) );
        Lattice2D::set_f( i ,jn,4, Lattice2D::get_fcoll( i,j,4 ) );
        Lattice2D::set_f( ip,jp,5, Lattice2D::get_fcoll( i,j,5 ) );
        Lattice2D::set_f( in,jp,6, Lattice2D::get_fcoll( i,j,6 ) );
        Lattice2D::set_f( in,jn,7, Lattice2D::get_fcoll( i,j,7 ) );
        Lattice2D::set_f( ip,jn,8, Lattice2D::get_fcoll( i,j,8 ) );
      }
    }
  }else{
    std::cout << "Only D2Q9 model (explicitly) implemented. Other models require more general streaming loop!" << endl;
  }
};


void Dynamics2D::streamBulk(int i_, int j_, int iq_){

Lattice2D::set_fcoll( i_ + Lattice2D::e[2*iq_] , j_ + Lattice2D::e[2*iq_+1] , iq_, Lattice2D::get_f(i_, j_, iq_ ) );

};


void Dynamics2D::streamBC_xn(int i_, int j_, int iq_, int corner_){
  if(corner_ == 0 && Lattice2D::e[2*iq_] > -1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[2*iq_] ,j_ + Lattice2D::e[2*iq_+1] , iq_, Lattice2D::get_f(i_, j_, iq_ ) );

  // Bottom corner
  if(corner_ == 1 && Lattice2D::e[2*iq_] > -1 && Lattice2D::e[2*iq_+1] > -1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[2*iq_] ,j_ + Lattice2D::e[2*iq_+1] , iq_, Lattice2D::get_f(i_, j_, iq_ ) );

  // Top corner
  if(corner_ == 2 && Lattice2D::e[2*iq_] > -1 && Lattice2D::e[2*iq_+1] < 1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[2*iq_] ,j_ + Lattice2D::e[2*iq_+1] , iq_, Lattice2D::get_f(i_, j_, iq_ ) );
};


void Dynamics2D::streamBC_xp(int i_, int j_, int iq_, int corner_){
  if(corner_ == 0 && Lattice2D::e[2*iq_] < 1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[2*iq_] ,j_ + Lattice2D::e[2*iq_+1] , iq_, Lattice2D::get_f(i_, j_, iq_ ) );

  // Bottom corner
  if(corner_ == 1 && Lattice2D::e[2*iq_] < 1 && Lattice2D::e[2*iq_+1] > -1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[2*iq_] ,j_ + Lattice2D::e[2*iq_+1] , iq_, Lattice2D::get_f(i_, j_, iq_ ) );

  // Top corner
  if(corner_ == 2 && Lattice2D::e[2*iq_] < 1 && Lattice2D::e[2*iq_+1] < 1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[2*iq_] ,j_ + Lattice2D::e[2*iq_+1] , iq_, Lattice2D::get_f(i_, j_, iq_ ) );
};


void Dynamics2D::streamBC_yn(int i_, int j_, int iq_, int corner_){
  if(corner_ == 0 && Lattice2D::e[2*iq_+1] > -1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[2*iq_] ,j_ + Lattice2D::e[2*iq_+1] , iq_, Lattice2D::get_f(i_, j_, iq_ ) );
};

void Dynamics2D::streamBC_yp(int i_, int j_, int iq_, int corner_){
  if(corner_ == 0 && Lattice2D::e[2*iq_+1] < 1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[2*iq_] ,j_ + Lattice2D::e[2*iq_+1] , iq_, Lattice2D::get_f(i_, j_, iq_ ) );
};
