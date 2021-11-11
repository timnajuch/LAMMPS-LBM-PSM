/*------------------------------------------------------ 
This file is part of the LBM-PSM project.

See the README file in the top-level LBM-PSM directory.

Tim Najuch, 2021
------------------------------------------------------*/

#include "dynamics2D.h"

Dynamics2D::Dynamics2D(int nx_, int ny_, int nz_, int q_, int decomposition_[3], int procCoordinates_[3], vector<double> origin_, vector<double> boxLength_, int dimension_) :
  Lattice2D(nx_, ny_, nz_, q_, decomposition_, procCoordinates_, origin_, boxLength_, dimension_)
{
  for(int i = 0; i < nx; ++i){
    for(int j = 0; j < ny; ++j){
      for(int k = 0; k < nz; ++k){
        for(int iq = 0; iq < q; ++iq){
          int ind_phys_1D = i * Lattice2D::ny * Lattice2D::nz + j * Lattice2D::nz + k;
          int ind_phys_2D = (i * Lattice2D::ny * Lattice2D::nz + j * Lattice2D::nz + k)*3;
          int ind_iq = i * Lattice2D::ny * Lattice2D::nz * Lattice2D::q + j * Lattice2D::nz * Lattice2D::q + k*Lattice2D::q + iq;

          Lattice2D::set_f0(i, j, k, iq, Dynamics2D::feq(iq, ind_phys_1D, ind_phys_2D, Lattice2D::rho, Lattice2D::u) );
          Lattice2D::set_f(i, j, k, iq, Dynamics2D::feq(iq, ind_phys_1D, ind_phys_2D, Lattice2D::rho, Lattice2D::u) );
          Lattice2D::set_fcoll(i, j, k, iq, Dynamics2D::feq(iq, ind_phys_1D, ind_phys_2D, Lattice2D::rho, Lattice2D::u) );
        }
      }
    }
  }
};


Dynamics2D::~Dynamics2D(){};
    

double Dynamics2D::feq(int iq_, int ind_phys_1D_, int ind_phys_2D_, vector<double> &rho_, vector<double> &u_){
  return rho_[ind_phys_1D_] * Lattice2D::w[iq_] * 
        (1.0 + ( Lattice2D::e[3*iq_] * u_[ind_phys_2D_] + Lattice2D::e[3*iq_+1] * u_[ind_phys_2D_+1] + Lattice2D::e[3*iq_+2] * u_[ind_phys_2D_+2]) / (pow(Lattice2D::cs, 2.0))
        + pow( Lattice2D::e[3*iq_] * u_[ind_phys_2D_] + Lattice2D::e[3*iq_+1] * u_[ind_phys_2D_+1] + Lattice2D::e[3*iq_+2] * u_[ind_phys_2D_+2] , 2.0) / (2.0*pow(Lattice2D::cs, 4.0))
        - 1.0/2.0 * ( u_[ind_phys_2D_] * u_[ind_phys_2D_] + u_[ind_phys_2D_+1] * u_[ind_phys_2D_+1] + u_[ind_phys_2D_+2] * u_[ind_phys_2D_+2]) / (pow(Lattice2D::cs,2.0) ) );
}


double Dynamics2D::feq(int iq_, int ind_phys_1D_, int ind_phys_2D_, double rho, vector<double> u){
  return rho * Lattice2D::w[iq_] * 
        (1.0 + ( Lattice2D::e[3*iq_] * u[0] + Lattice2D::e[3*iq_+1] * u[1] + Lattice2D::e[3*iq_+2] * u[2]) / (pow(Lattice2D::cs, 2.0))
        + pow( Lattice2D::e[3*iq_] * u[0] + Lattice2D::e[3*iq_+1] * u[1] + Lattice2D::e[iq_*3+2] * u[2], 2.0) / (2.0*pow(Lattice2D::cs, 4.0))
        - 1.0/2.0 * ( u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) / (pow(Lattice2D::cs,2.0) ) );
}


void Dynamics2D::streaming_periodic_2D(){
  int in, ip, jn, jp, kn, kp;
  in = 0;
  ip = 0;
  jn = 0;
  jp = 0;
  kn = 0;
  kp = 0;
  if(Lattice2D::q == 9){

    for(int i = 0; i < Lattice2D::nx; ++i){
      in = (i>0   )               ? (i-1):(Lattice2D::nx-1);
      ip = (i<Lattice2D::nx-1) ? (i+1):(0   );

      for(int j = 0; j < Lattice2D::ny; ++j){
        jn = (j>0   )               ? (j-1):(Lattice2D::ny-1);
        jp = (j<Lattice2D::ny-1) ? (j+1):(0   );

        Lattice2D::set_f( i ,j ,0, 0, Lattice2D::get_fcoll( i,j,0,0 ) );
        Lattice2D::set_f( ip,j ,0, 1, Lattice2D::get_fcoll( i,j,0,1 ) );
        Lattice2D::set_f( i ,jp,0, 2, Lattice2D::get_fcoll( i,j,0,2 ) );
        Lattice2D::set_f( in,j ,0, 3, Lattice2D::get_fcoll( i,j,0,3 ) );
        Lattice2D::set_f( i ,jn,0, 4, Lattice2D::get_fcoll( i,j,0,4 ) );
        Lattice2D::set_f( ip,jp,0, 5, Lattice2D::get_fcoll( i,j,0,5 ) );
        Lattice2D::set_f( in,jp,0, 6, Lattice2D::get_fcoll( i,j,0,6 ) );
        Lattice2D::set_f( in,jn,0, 7, Lattice2D::get_fcoll( i,j,0,7 ) );
        Lattice2D::set_f( ip,jn,0, 8, Lattice2D::get_fcoll( i,j,0,8 ) );
      }
    }
  }
  else if(Lattice2D::q == 19){
    for(int i = 0; i < Lattice2D::nx; ++i){
      //in = (i>0   )            ? (i-1):(Lattice2D::nx-1);
      //ip = (i<Lattice2D::nx-1) ? (i+1):(0   );

      for(int j = 0; j < Lattice2D::ny; ++j){
        //jn = (j>0   )            ? (j-1):(Lattice2D::ny-1);
        //jp = (j<Lattice2D::ny-1) ? (j+1):(0   );

        for(int k = 0; j < Lattice2D::nz; ++k){
          //kn = (k>0   )            ? (k-1):(Lattice2D::nz-1);
          //kp = (k<Lattice2D::nz-1) ? (k+1):(0   );

          for(int iq = 0; iq < Lattice2D::q; ++iq){

            int iShift = Lattice2D::e[3*iq];
            if(iShift > 0 && i < Lattice2D::nx-1)
            { 
              iShift += i;
            }else if(iShift > 0 && i > Lattice2D::nx-2)
            {
              iShift = 0;
            }else if(iShift < 0 && i > 0)
            {
              iShift = i - iShift;
            }else if(iShift < 0 && i < 1)
            {
              iShift = Lattice2D::nx-1;
            }

            int jShift = Lattice2D::e[3*iq+1];
            if(jShift > 0 && j < Lattice2D::ny-1)
            { 
              jShift += j;
            }else if(jShift > 0 && j > Lattice2D::ny-2)
            {
              jShift = 0;
            }else if(jShift < 0 && j > 0)
            {
              jShift = j - jShift;
            }else if(jShift < 0 && j < 1)
            {
              jShift = Lattice2D::ny-1;
            }
        
            int kShift = Lattice2D::e[3*iq+2];
            if(kShift > 0 && k < Lattice2D::nz-1)
            { 
              kShift += k;
            }else if(kShift > 0 && k > Lattice2D::nz-2)
            {
              kShift = 0;
            }else if(kShift < 0 && k > 0)
            {
              kShift = j - kShift;
            }else if(kShift < 0 && k < 1)
            {
              kShift = Lattice2D::nz-1;
            }
          
            Lattice2D::set_f( iShift, jShift, kShift, 0, Lattice2D::get_fcoll( i,j,k,iq ) );

/*
            Lattice2D::set_f( i ,j ,0, Lattice2D::get_fcoll( i,j,k,0 ) );
            Lattice2D::set_f( ip,j ,1, Lattice2D::get_fcoll( i,j,k,1 ) );
            Lattice2D::set_f( i ,jp,2, Lattice2D::get_fcoll( i,j,k,2 ) );
            Lattice2D::set_f( in,j ,3, Lattice2D::get_fcoll( i,j,k,3 ) );
            Lattice2D::set_f( i ,jn,4, Lattice2D::get_fcoll( i,j,k,4 ) );
            Lattice2D::set_f( ip,jp,5, Lattice2D::get_fcoll( i,j,k,5 ) );
            Lattice2D::set_f( in,jp,6, Lattice2D::get_fcoll( i,j,k,6 ) );
            Lattice2D::set_f( in,jn,7, Lattice2D::get_fcoll( i,j,k,7 ) );
            Lattice2D::set_f( ip,jn,8, Lattice2D::get_fcoll( i,j,k,8 ) );
*/
          }
        }
      }
    }
  
  }else{
    std::cout << "Only D2Q9 and D3Q19 model (explicitly) implemented. Other models require a more general streaming loop!" << endl;
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

        Lattice2D::set_f( i ,j ,0,0, Lattice2D::get_fcoll( i,j,0,0 ) );
        Lattice2D::set_f( ip,j ,0,1, Lattice2D::get_fcoll( i,j,0,1 ) );
        Lattice2D::set_f( i ,jp,0,2, Lattice2D::get_fcoll( i,j,0,2 ) );
        Lattice2D::set_f( in,j ,0,3, Lattice2D::get_fcoll( i,j,0,3 ) );
        Lattice2D::set_f( i ,jn,0,4, Lattice2D::get_fcoll( i,j,0,4 ) );
        Lattice2D::set_f( ip,jp,0,5, Lattice2D::get_fcoll( i,j,0,5 ) );
        Lattice2D::set_f( in,jp,0,6, Lattice2D::get_fcoll( i,j,0,6 ) );
        Lattice2D::set_f( in,jn,0,7, Lattice2D::get_fcoll( i,j,0,7 ) );
        Lattice2D::set_f( ip,jn,0,8, Lattice2D::get_fcoll( i,j,0,8 ) );
      }
    }
  }
  else if(Lattice2D::q == 19){ // todo not entirely correct
    for(int i = 0; i < Lattice2D::nx; ++i){
      //in = (i>0   )            ? (i-1):(Lattice2D::nx-1);
      //ip = (i<Lattice2D::nx-1) ? (i+1):(0   );

      for(int j = 0; j < Lattice2D::ny; ++j){
        //jn = (j>0   )            ? (j-1):(Lattice2D::ny-1);
        //jp = (j<Lattice2D::ny-1) ? (j+1):(0   );

        for(int k = 0; j < Lattice2D::nz; ++k){
          //kn = (k>0   )            ? (k-1):(Lattice2D::nz-1);
          //kp = (k<Lattice2D::nz-1) ? (k+1):(0   );

          for(int iq = 0; iq < Lattice2D::q; ++iq){

            int iShift = Lattice2D::e[3*iq];
            if(iShift > 0 && i < Lattice2D::nx-1)
            { 
              iShift += i;
            }else if(iShift > 0 && i > Lattice2D::nx-2)
            {
              iShift = 0;
            }else if(iShift < 0 && i > 0)
            {
              iShift = i - iShift;
            }else if(iShift < 0 && i < 1)
            {
              iShift = Lattice2D::nx-1;
            }

            int jShift = Lattice2D::e[3*iq+1];
            if(jShift > 0 && j < Lattice2D::ny-1)
            { 
              jShift += j;
            }else if(jShift > 0 && j > Lattice2D::ny-2)
            {
              jShift = 0;
            }else if(jShift < 0 && j > 0)
            {
              jShift = j - jShift;
            }else if(jShift < 0 && j < 1)
            {
              jShift = Lattice2D::ny-1;
            }
        
            int kShift = Lattice2D::e[3*iq+2];
            if(kShift > 0 && k < Lattice2D::nz-1)
            { 
              kShift += k;
            }else if(kShift > 0 && k > Lattice2D::nz-2)
            {
              kShift = 0;
            }else if(kShift < 0 && k > 0)
            {
              kShift = j - kShift;
            }else if(kShift < 0 && k < 1)
            {
              kShift = Lattice2D::nz-1;
            }
          
            //Lattice2D::set_f( iShift, jShift, kShift, 0, Lattice2D::get_fcoll( i,j,k,iq ) );
            Lattice2D::set_f( iShift, jShift, kShift, iq, Lattice2D::get_fcoll( i,j,k,iq ) );

/*
            Lattice2D::set_f( i ,j ,0, Lattice2D::get_fcoll( i,j,k,0 ) );
            Lattice2D::set_f( ip,j ,1, Lattice2D::get_fcoll( i,j,k,1 ) );
            Lattice2D::set_f( i ,jp,2, Lattice2D::get_fcoll( i,j,k,2 ) );
            Lattice2D::set_f( in,j ,3, Lattice2D::get_fcoll( i,j,k,3 ) );
            Lattice2D::set_f( i ,jn,4, Lattice2D::get_fcoll( i,j,k,4 ) );
            Lattice2D::set_f( ip,jp,5, Lattice2D::get_fcoll( i,j,k,5 ) );
            Lattice2D::set_f( in,jp,6, Lattice2D::get_fcoll( i,j,k,6 ) );
            Lattice2D::set_f( in,jn,7, Lattice2D::get_fcoll( i,j,k,7 ) );
            Lattice2D::set_f( ip,jn,8, Lattice2D::get_fcoll( i,j,k,8 ) );
*/
          }
        }
      }
    }
  
  }else{
    std::cout << "Only D2Q9 model (explicitly) implemented. Other models require more general streaming loop!" << endl;
  }
};


void Dynamics2D::streamBulk(int i_, int j_, int k_, int iq_){

//Lattice2D::set_fcoll( i_ + Lattice2D::e[2*iq_] , j_ + Lattice2D::e[2*iq_+1] , iq_, Lattice2D::get_f(i_, j_, iq_ ) );
//Lattice2D::set_fcoll( i_ + Lattice2D::e[3*iq_] , j_ + Lattice2D::e[3*iq_+1] , iq_, k_ + Lattice2D::e[3*iq_+2], Lattice2D::get_f(i_, j_, k_, iq_ ) );
Lattice2D::set_fcoll( i_ + Lattice2D::e[3*iq_] , j_ + Lattice2D::e[3*iq_+1] , k_ + Lattice2D::e[3*iq_+2], iq_, Lattice2D::get_f(i_, j_, k_, iq_ ) );

};


// TODO
void Dynamics2D::streamBC_xn(int i_, int j_, int iq_, int corner_){
  int k_ = 0;
  if(corner_ == 0 && Lattice2D::e[3*iq_] > -1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[3*iq_], j_ + Lattice2D::e[3*iq_+1], k_ + Lattice2D::e[3*iq_+2], iq_, Lattice2D::get_f(i_, j_, k_, iq_ ) );
    //Lattice2D::set_fcoll( i_ + Lattice2D::e[2*iq_] ,j_ + Lattice2D::e[2*iq_+1] , iq_, Lattice2D::get_f(i_, j_, iq_ ) );

  // Bottom corner/edge
  if(corner_ == 1 && Lattice2D::e[3*iq_] > -1 && Lattice2D::e[3*iq_+1] > -1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[3*iq_], j_ + Lattice2D::e[3*iq_+1], k_ + Lattice2D::e[3*iq_+2], iq_, Lattice2D::get_f(i_, j_, k_, iq_ ) );
//    Lattice2D::set_fcoll( i_ + Lattice2D::e[2*iq_] ,j_ + Lattice2D::e[2*iq_+1] , iq_, Lattice2D::get_f(i_, j_, iq_ ) );

  // Top corner/edge
  if(corner_ == 2 && Lattice2D::e[3*iq_] > -1 && Lattice2D::e[3*iq_+1] < 1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[3*iq_], j_ + Lattice2D::e[3*iq_+1], k_ + Lattice2D::e[3*iq_+2], iq_, Lattice2D::get_f(i_, j_, k_, iq_ ) );
    //Lattice2D::set_fcoll( i_ + Lattice2D::e[2*iq_] ,j_ + Lattice2D::e[2*iq_+1] , iq_, Lattice2D::get_f(i_, j_, iq_ ) );
/*
  if(domain->dimension == 3){
  // Bottom back corner
    if(corner_ == 3 && Lattice2D::e[3*iq_] > -1 && Lattice2D::e[3*iq_+1] > -1)
      Lattice2D::set_fcoll( i_ + Lattice2D::e[3*iq_], j_ + Lattice2D::e[3*iq_+1], k_ + Lattice2D::e[3*iq_+2], iq_, Lattice2D::get_f(i_, j_, k_, iq_ ) );
*/
  //}
};

void Dynamics2D::streamBC(int i_, int j_, int k_, int iq_)
{
  if(   Lattice2D::e[3*iq_]+i_ >= 0 && Lattice2D::e[3*iq_]+i_ <= Lattice2D::nx-1
    &&  Lattice2D::e[3*iq_+1]+j_ >= 0 && Lattice2D::e[3*iq_+1]+j_ <= Lattice2D::ny-1
    &&  Lattice2D::e[3*iq_+2]+k_ >= 0 && Lattice2D::e[3*iq_+2]+k_ <= Lattice2D::nz-1)
    {
      Lattice2D::set_fcoll( i_ + Lattice2D::e[3*iq_], j_ + Lattice2D::e[3*iq_+1], k_ + Lattice2D::e[3*iq_+2], iq_, Lattice2D::get_f(i_, j_, k_, iq_ ) );
    }
}

void Dynamics2D::streamBC_xp(int i_, int j_, int iq_, int corner_){
  if(corner_ == 0 && Lattice2D::e[3*iq_] < 1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[3*iq_] ,j_ + Lattice2D::e[3*iq_+1] ,0, iq_, Lattice2D::get_f(i_, j_, 0,iq_ ) );

  // Bottom corner
  if(corner_ == 1 && Lattice2D::e[3*iq_] < 1 && Lattice2D::e[3*iq_+1] > -1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[3*iq_] ,j_ + Lattice2D::e[3*iq_+1] ,0, iq_, Lattice2D::get_f(i_, j_, 0,iq_ ) );

  // Top corner
  if(corner_ == 2 && Lattice2D::e[3*iq_] < 1 && Lattice2D::e[3*iq_+1] < 1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[3*iq_] ,j_ + Lattice2D::e[3*iq_+1] ,0, iq_, Lattice2D::get_f(i_, j_, 0,iq_ ) );
};


void Dynamics2D::streamBC_yn(int i_, int j_, int iq_, int corner_){
  if(corner_ == 0 && Lattice2D::e[3*iq_+1] > -1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[3*iq_] ,j_ + Lattice2D::e[3*iq_+1] , 0, iq_, Lattice2D::get_f(i_, j_, 0,iq_ ) );
};

void Dynamics2D::streamBC_yp(int i_, int j_, int iq_, int corner_){
  if(corner_ == 0 && Lattice2D::e[3*iq_+1] < 1)
    Lattice2D::set_fcoll( i_ + Lattice2D::e[3*iq_] ,j_ + Lattice2D::e[3*iq_+1] , 0, iq_, Lattice2D::get_f(i_, j_, 0,iq_ ) );
};
