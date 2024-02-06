/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Basis class: pair_lubricate.cpp / pair_lubricate_GRM.cpp
   Author: Tim Najuch (The University of Edinburgh) 
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_GRM_LBDEM.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "fix_wall.h"
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

// same as fix_deform.cpp

enum{NO_REMAP,X_REMAP,V_REMAP};

// same as fix_wall.cpp

enum{NONE=0,EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

PairLubricateGRMLBDEM::PairLubricateGRMLBDEM(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;

  // set comm size needed by this Pair

  comm_forward = 6;

  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

PairLubricateGRMLBDEM::~PairLubricateGRMLBDEM()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(cut_inner);
  }
}

/* ---------------------------------------------------------------------- */

void PairLubricateGRMLBDEM::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,beta0,beta1,betaInv0,betaInv1,radi,radj;
  double vi[3],vj[3],wi[3],wj[3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  double lamda[3];

  double radijmean;
  double X_A_11, X_A_12;
  double Y_A_11, Y_A_12;
  double Y_B_11, Y_B_12, Y_B_21;
  double Y_C_11, Y_C_12;
  double X_A_11_N, X_A_12_N;
  double Y_A_11_N, Y_A_12_N;
  double Y_B_11_N, Y_B_12_N, Y_B_21_N;
  double Y_C_11_N, Y_C_12_N;
  double vin, vin1, vin2, vin3;
  double vjn, vjn1, vjn2, vjn3;
  double vit1, vit2, vit3;
  double vjt1, vjt2, vjt3;
  double viXei1, viXei2, viXei3;
  double vjXei1, vjXei2, vjXei3;
  double win, wjn;
  double win1, win2, win3, wit1, wit2, wit3;
  double wjn1, wjn2, wjn3, wjt1, wjt2, wjt3;
  double ei1Xwi, ei2Xwi, ei3Xwi;
  double ei1Xwj, ei2Xwj, ei3Xwj;

  double vxmu2f = force->vxmu2f;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    radi = radius[i];

    // translational velocity
    vi[0] = -v[i][0];
    vi[1] = -v[i][1];
    vi[2] = -v[i][2];

    // angular velocity
    wi[0] = -omega[i][0];
    wi[1] = -omega[i][1];
    wi[2] = -omega[i][2];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      // pointing from particle i to j
      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;
      rsq = delx*delx + dely*dely + delz*delz;

      jtype = type[j];
      radj = atom->radius[j];

      double onedx = sqrt(cutsq[itype][jtype]);
      double h_sepN = onedx*cutoff_factor;
      
      r = sqrt(rsq);
      h_sep = r - radi-radj;

      if (h_sep < h_sepN) {

        wj[0] = -omega[j][0];
        wj[1] = -omega[j][1];
        wj[2] = -omega[j][2];
       
        vj[0] = -v[j][0];
        vj[1] = -v[j][1];
        vj[2] = -v[j][2];

        // if less than the minimum gap use the minimum gap instead
        if (h_sep < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype];

        // scale h_sep by mean radius radijmean
		    radijmean = 0.5*(radi+radj);
        h_sep = h_sep/radijmean;
        beta0 = radj/radi;
        beta1 = 1.0 + beta0;
        betaInv0 = 1.0/beta0;
        betaInv1 = 1.0 + betaInv0;
	    
        h_sepN = h_sepN/radijmean;

/* -------------------------------------------------------------------------
        Grand resistance matrix notation.
        1) Jeffrey and Onishi (1984)
        "Resistance functions for two spheres in low-Reynolds-number flow"
        2) Kim & Karilla "Microhydrodynamics" (chapter 7 and 11)

                                                        / U_1 - U(x_1) \
        / F_1 \   / A_11 A_12 Bt_11 Bt_12 Gt_11 Gt_12\ | U_2 - U(x_2) |
        | F_2 |   | A_21 A_22 Bt_21 Bt_22 Gt_21 Gt_22| | w_1 - w(x_1) |
        | L_1 | = | B_11 B_12 C_11  C_12  Ht_11 Ht_12| | w_2 - w(x_2) |
        \ L_2 /   \ B_21 B_22 C_21  C_22  Ht_21 Ht_22/ |    E_infty   |
                                                        \    E_infty   /

        where

        A_(ij)^(ab) = X_(ab)^A*e_i*e_j + Y_(ab)^A (d_(ij) - e_i*e_j)
        B_(ij)^(ab) = Bt_(ji)^(ba) = ...
                    Y_(ab)^B*e_(ijk)*e_k
        C_(ij)^(ab) = X_(ab)^C*e_i*e_j + Y_(ab)^C*(d_(ij) - e_i*e_j)
        G_(ijk)^(ab) = Gt_(kij)_^(ba) = ...
                X_(ab)^G*(e_i*e_j - 1/3*d_(ij))e_k ...
              + Y_(ab)^G*(e_i*d_(jk) + e_j*d_(ik) - 2*e_i*e_j*e_k)
        H_(ijk)^(ab) = Ht_(kij)^(ab) = ...
                Y_(ab)^H*(e_(ikl)*e_l*e_j + e_(ijkl)*e_l*e_i)

        e_i: Unit vector along center to center axis
        d_(ij): Dirac symbol
        e_(ijk): Levi-Cita symbol

        // Twisting terms X_C are not implemented because:
        // a) Rieman zeta function would have to be evaluated
        //      during run-time
        // b) Supposedly neglectable contribution

------------------------------------------------------------------------- */
        X_A_11 = 0.0; X_A_12 = 0.0;
        Y_A_11 = 0.0; Y_A_12 = 0.0;
        Y_B_11 = 0.0; Y_B_12 = 0.0; Y_B_21 = 0.0;
        Y_C_11 = 0.0; Y_C_12 = 0.0;
        X_A_11_N = 0.0; X_A_12_N = 0.0;
        Y_A_11_N = 0.0; Y_A_12_N = 0.0;
        Y_B_11_N = 0.0; Y_B_12_N = 0.0; Y_B_21_N = 0.0;
        Y_C_11_N = 0.0; Y_C_12_N = 0.0;

        double b0p2 = beta0*beta0;
        double b0p3 = b0p2*beta0;
        double b0p4 = b0p3*beta0;
        double b1p2 = beta1*beta1;
        double b1p3 = b1p2*beta1;
        double b1p4 = b1p3*beta1;
        
        double bInv0p2 = betaInv0*betaInv0;
        double bInv0p3 = bInv0p2*betaInv0;
        double bInv0p4 = bInv0p3*betaInv0;
        double bInv1p2 = betaInv1*betaInv1;
        double bInv1p3 = bInv1p2*betaInv1;
        double bInv1p4 = bInv1p3*betaInv1;
        
        double log1oh = log(1.0/h_sep);
        double radip2 = radi*radi;
        double radip3 = radip2*radi;
        double radjp2 = radj*radj;
        double radjp3 = radjp2*radj;

        double log1ohN = log(1.0/h_sepN);

        // Squeezing terms X_A
        if(flagGRM == 1){
          X_A_11 = 2.0*b0p2/b1p3/h_sep;
          X_A_11 *= 6.0*MY_PI*radi;

          X_A_12 = -X_A_11;
        
          X_A_11_N = 2.0*b0p2/b1p3/h_sepN;
          X_A_11_N *= 6.0*MY_PI*radi;

          X_A_12_N = -X_A_11_N;
        }
              
        // Shearing terms Y_A
        if(flagGRM == 2){
          Y_A_11 = 4.0*beta0*(2.0 + beta0 + 2.0*b0p2)/(15.0*b1p3)*log1oh;
          Y_A_11 *= 6.0*MY_PI*radi;

          Y_A_12 = -Y_A_11;

          Y_A_11_N = 4.0*beta0*(2.0 + beta0 + 2.0*b0p2)/(15.0*b1p3)*log1ohN;
          Y_A_11_N *= 6.0*MY_PI*radi;

          Y_A_12_N = -Y_A_11_N;
        }

        // Shearing terms Y_B
        if(flagGRM == 3){
          Y_B_11 = -beta0*(4.0 + beta0)/(5.0*b1p2)*log1oh;
          Y_B_11 *= 4.0*MY_PI*radip2;

          Y_B_12 = -Y_B_11;

          Y_B_11_N = -beta0*(4.0 + beta0)/(5.0*b1p2)*log1ohN;
          Y_B_11_N *= 4.0*MY_PI*radip2;

          Y_B_12_N = -Y_B_11_N;
        }

        // Terms Y_B for forces from rotations
        if(flagGRM == 4){
          Y_B_11 = -beta0*(4.0 + beta0)/(5.0*b1p2)*log1oh;
          Y_B_11 *= 4.0*MY_PI*radip2;

          Y_B_21 = betaInv0*(4.0 + betaInv0)/(5.0*bInv1p2)*log1oh;
          Y_B_21 *= -4.0*MY_PI*radjp2;

          Y_B_11_N = -beta0*(4.0 + beta0)/(5.0*b1p2)*log1ohN;
          Y_B_11_N *= 4.0*MY_PI*radip2;

          Y_B_21_N = betaInv0*(4.0 + betaInv0)/(5.0*bInv1p2)*log1ohN;
          Y_B_21_N *= -4.0*MY_PI*radjp2;
        }

        // Rotational Y_C terms
        if(flagGRM == 5){
          Y_C_11 = 2.0*beta0/(5.0*beta1)*log1oh;
          Y_C_11 *= 8*MY_PI*radip3;

          Y_C_12 = b0p2/(10.0*beta1)*log1oh;
          Y_C_12 *= 8*MY_PI*radip3;

          Y_C_11_N = 2.0*beta0/(5.0*beta1)*log1ohN;
          Y_C_11_N *= 8*MY_PI*radip3;

          Y_C_12_N = b0p2/(10.0*beta1)*log1ohN;
          Y_C_12_N *= 8*MY_PI*radip3;
        }

        // Initialise different variables
        vin = 0.0; vin1 = 0.0; vin2 = 0.0; vin3 = 0.0;
        vit1 = 0.0; vit2 = 0.0; vit3 = 0.0;
        viXei1 = 0.0; viXei2 = 0.0; viXei3 = 0.0;
        win = 0.0; win1 = 0.0; win2 = 0.0; win3 = 0.0; 
        ei1Xwi = 0.0; ei2Xwi = 0.0; ei3Xwi = 0.0;

        vjn = 0.0; vjn1 = 0.0; vjn2 = 0.0; vjn3 = 0.0;
        vjt1 = 0.0; vjt2 = 0.0; vjt3 = 0.0;
        vjXei1 = 0.0; vjXei2 = 0.0; vjXei3 = 0.0;
        wjn = 0.0; wjn1 = 0.0; wjn2 = 0.0; wjn3 = 0.0; 
        ei1Xwj = 0.0; ei2Xwj = 0.0; ei3Xwj = 0.0;

        // Unit vector pointing along center-to-center line (from i to j)
        double ei1 = delx/r;
        double ei2 = dely/r;
        double ei3 = delz/r;

        // Velocity calculations
            
        if( (flagGRM == 1) || (flagGRM == 2)  ){
          vin = vi[0]*ei1 + vi[1]*ei2 + vi[2]*ei3;
          vin1 = vin*ei1;
          vin2 = vin*ei2;
          vin3 = vin*ei3;

          vjn = vj[0]*ei1 + vj[1]*ei2 + vj[2]*ei3;
          vjn1 = vjn*ei1;
          vjn2 = vjn*ei2;
          vjn3 = vjn*ei3;
        }
        
        if( flagGRM == 2 ){
          vit1 = vi[0] - vin1;
          vit2 = vi[1] - vin2;
          vit3 = vi[2] - vin3;

          vjt1 = vj[0] - vjn1;
          vjt2 = vj[1] - vjn2;
          vjt3 = vj[2] - vjn3;
        }

        if( flagGRM == 3 ){
          viXei1 = (vi[1]*ei3 - vi[2]*ei2);
          viXei2 = (vi[2]*ei1 - vi[0]*ei3);
          viXei3 = (vi[0]*ei2 - vi[1]*ei1);

          vjXei1 = (vj[1]*ei3 - vj[2]*ei2);
          vjXei2 = (vj[2]*ei1 - vj[0]*ei3);
          vjXei3 = (vj[0]*ei2 - vj[1]*ei1);
        }

        if( (flagGRM == 4) || (flagGRM == 5) ){
          win = wi[0]*ei1 + wi[1]*ei2 + wi[2]*ei3;
          win1 = win*ei1; 
          win2 = win*ei2; 
          win3 = win*ei3; 

          wit1 = wi[0] - win1;
          wit2 = wi[1] - win2;
          wit3 = wi[2] - win3;

          wjn = wj[0]*ei1 + wj[1]*ei2 + wj[2]*ei3;
          wjn1 = wjn*ei1; 
          wjn2 = wjn*ei2; 
          wjn3 = wjn*ei3; 

          wjt1 = wj[0] - wjn1;
          wjt2 = wj[1] - wjn2;
          wjt3 = wj[2] - wjn3;
        }		

        if( flagGRM == 4 ){
          ei1Xwi = ei2*wi[2] - ei3*wi[1];
          ei2Xwi = ei3*wi[0] - ei1*wi[2];
          ei3Xwi = ei1*wi[1] - ei2*wi[0];

          ei1Xwj = ei2*wj[2] - ei3*wj[1];
          ei2Xwj = ei3*wj[0] - ei1*wj[2];
          ei3Xwj = ei1*wj[1] - ei2*wj[0];
        }
            
        // Force calculations

        fx = 0.0; fy = 0.0; fz = 0.0;
        tx = 0.0; ty = 0.0; tz = 0.0;

        // Squeezing force correction
        if(flagGRM == 1){
          fx = mu*((X_A_11*vin1 + X_A_12*vjn1) - (X_A_11_N*vin1 + X_A_12_N*vjn1));
          fy = mu*((X_A_11*vin2 + X_A_12*vjn2) - (X_A_11_N*vin2 + X_A_12_N*vjn2));
          fz = mu*((X_A_11*vin3 + X_A_12*vjn3) - (X_A_11_N*vin3 + X_A_12_N*vjn3));
        }

        // Shearing force correction
        if(flagGRM == 2){
          fx = mu*((Y_A_11*vit1 + Y_A_12*vjt1) - (Y_A_11_N*vit1 + Y_A_12_N*vjt1));
          fy = mu*((Y_A_11*vit2 + Y_A_12*vjt2) - (Y_A_11_N*vit2 + Y_A_12_N*vjt2));
          fz = mu*((Y_A_11*vit3 + Y_A_12*vjt3) - (Y_A_11_N*vit3 + Y_A_12_N*vjt3));
        }

        // Force arising from rotation
        if(flagGRM == 4){
          fx = mu*((Y_B_11*ei1Xwi + Y_B_21*ei1Xwj) - (Y_B_11_N*ei1Xwi + Y_B_21_N*ei1Xwj));
          fy = mu*((Y_B_11*ei2Xwi + Y_B_21*ei2Xwj) - (Y_B_11_N*ei2Xwi + Y_B_21_N*ei2Xwj));
          fz = mu*((Y_B_11*ei3Xwi + Y_B_21*ei3Xwj) - (Y_B_11_N*ei3Xwi + Y_B_21_N*ei3Xwj));
        }

        // scale force for appropriate units 
        fx *= vxmu2f;
        fy *= vxmu2f;
        fz *= vxmu2f;

        // add to total force
        if( flagGRM == 1 || flagGRM == 2 || flagGRM == 4 ){
          f[i][0] += fx;
          f[i][1] += fy;
          f[i][2] += fz;
        }

        // Torque correction due to shearing motion 
        if(flagGRM == 3){
          tx = mu*((Y_B_11*viXei1 + Y_B_12*vjXei1) - (Y_B_11_N*viXei1 + Y_B_12_N*vjXei1));
          ty = mu*((Y_B_11*viXei2 + Y_B_12*vjXei2) - (Y_B_11_N*viXei2 + Y_B_12_N*vjXei2));
          tz = mu*((Y_B_11*viXei3 + Y_B_12*vjXei3) - (Y_B_11_N*viXei3 + Y_B_12_N*vjXei3));
            }

        // Torque due to rotation
        if(flagGRM == 5){
          tx = mu*((Y_C_11*wit1 + Y_C_12*wjt1) - (Y_C_11_N*wit1 + Y_C_12_N*wjt1));
          ty = mu*((Y_C_11*wit2 + Y_C_12*wjt2) - (Y_C_11_N*wit2 + Y_C_12_N*wjt2));
          tz = mu*((Y_C_11*wit3 + Y_C_12*wjt3) - (Y_C_11_N*wit3 + Y_C_12_N*wjt3));
        }

        // scale torque for appropriate units and add to total torque
        if( flagGRM == 3 || flagGRM == 5 ){
          torque[i][0] += vxmu2f*tx;
          torque[i][1] += vxmu2f*ty;
          torque[i][2] += vxmu2f*tz;
        }

        // set j = nlocal so that only I gets tallied

        if (evflag) ev_tally_xyz(i,nlocal,nlocal,0,
                                  0.0,0.0,fx,fy,fz,-delx,-dely,-delz);

      }
    }
  }

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLubricateGRMLBDEM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(cut_inner,n+1,n+1,"pair:cut_inner");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLubricateGRMLBDEM::settings(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Illegal pair_style command");

  mu = utils::numeric(FLERR,arg[0],false,lmp);
  cut_inner_global = utils::numeric(FLERR,arg[1],false,lmp);
  cut_global = utils::numeric(FLERR,arg[2],false,lmp);
  flagGRM  = utils::inumeric(FLERR,arg[3],false,lmp);
  cutoff_factor =  utils::numeric(FLERR,arg[4],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_inner[i][j] = cut_inner_global;
          cut[i][j] = cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLubricateGRMLBDEM::coeff(int narg, char **arg)
{
  if (narg != 2 && narg != 4)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double cut_inner_one = cut_inner_global;
  double cut_one = cut_global;
  if (narg == 4) {
    cut_inner_one = utils::numeric(FLERR,arg[2],false,lmp);
    cut_one = utils::numeric(FLERR,arg[3],false,lmp);
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut_inner[i][j] = cut_inner_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLubricateGRMLBDEM::init_style()
{
  if (force->newton_pair == 1)
    error->all(FLERR,"Pair lubricate/GRM/LBDEM requires newton pair off");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,
               "Pair lubricate/GRM/LBDEM requires ghost atoms store velocity");
  if (!atom->sphere_flag)
    error->all(FLERR,"Pair lubricate/GRM/LBDEM requires atom style sphere");

  // ensure all particles are finite-size
  // for pair hybrid, should limit test to types using the pair style

  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (radius[i] == 0.0)
      error->one(FLERR,"Pair lubricate/GRM/LBDEM requires extended particles");

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // check for fix deform, if exists it must use "remap v"
  // If box will change volume, set appropriate flag so that volume
  // and v.f. corrections are re-calculated at every step.

  // If available volume is different from box volume
  // due to walls, set volume appropriately; if walls will
  // move, set appropriate flag so that volume and v.f. corrections
  // are re-calculated at every step.

  shearing = flagdeform = flagwall = 0;
  for (int i = 0; i < modify->nfix; i++){
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      shearing = flagdeform = 1;
      if (((FixDeform *) modify->fix[i])->remapflag != Domain::V_REMAP)
        error->all(FLERR,"Using pair lubricate with inconsistent "
                   "fix deform remap option");
    }
    if (strstr(modify->fix[i]->style,"wall") != nullptr) {
      if (flagwall)
        error->all(FLERR,
                   "Cannot use multiple fix wall commands with "
                   "pair lubricate/GRM/LBDEM");
      flagwall = 1; // Walls exist
      wallfix = (FixWall *) modify->fix[i];
      if (wallfix->xflag) flagwall = 2; // Moving walls exist
    }

    if (strstr(modify->fix[i]->style,"wall") != NULL){
      flagwall = 1; // Walls exist
      if (((FixWall *) modify->fix[i])->xflag ) {
        flagwall = 2; // Moving walls exist
        wallfix = (FixWall *) modify->fix[i];
      }
    }
  }

  // check for fix deform, if exists it must use "remap v"

  shearing = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      shearing = 1;
      if (((FixDeform *) modify->fix[i])->remapflag != Domain::V_REMAP)
        error->all(FLERR,"Using pair lubricate/GRM/LBDEM with inconsistent "
                   "fix deform remap option");
    }

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLubricateGRMLBDEM::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut_inner[i][j] = mix_distance(cut_inner[i][i],cut_inner[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  cut_inner[j][i] = cut_inner[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLubricateGRMLBDEM::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cut_inner[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLubricateGRMLBDEM::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&cut_inner[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&cut_inner[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLubricateGRMLBDEM::write_restart_settings(FILE *fp)
{
  fwrite(&mu,sizeof(double),1,fp);
  fwrite(&cut_inner_global,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&flagGRM,sizeof(int),1,fp);
  fwrite(&cutoff_factor,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLubricateGRMLBDEM::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&mu,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_inner_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&flagGRM,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&cutoff_factor,sizeof(double),1,fp,nullptr,error);
  }
  MPI_Bcast(&mu,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_inner_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&flagGRM,1,MPI_INT,0,world);
  MPI_Bcast(&cutoff_factor,1,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

int PairLubricateGRMLBDEM::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
  int i,j,m;

  double **v = atom->v;
  double **omega = atom->omega;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = v[j][0];
    buf[m++] = v[j][1];
    buf[m++] = v[j][2];
    buf[m++] = omega[j][0];
    buf[m++] = omega[j][1];
    buf[m++] = omega[j][2];
  }

  return 6;
}

/* ---------------------------------------------------------------------- */

void PairLubricateGRMLBDEM::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  double **v = atom->v;
  double **omega = atom->omega;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    omega[i][0] = buf[m++];
    omega[i][1] = buf[m++];
    omega[i][2] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   check if name is recognized, return integer index for that name
   if name not recognized, return -1
   if type pair setting, return -2 if no type pairs are set
------------------------------------------------------------------------- */

int PairLubricateGRMLBDEM::pre_adapt(char *name, int ilo, int ihi, int jlo, int jhi)
{
  if (strcmp(name,"mu") == 0) return 0;
  return -1;
}

/* ----------------------------------------------------------------------
   adapt parameter indexed by which
   change all pair variables affected by the reset parameter
   if type pair setting, set I-J and J-I coeffs
------------------------------------------------------------------------- */

void PairLubricateGRMLBDEM::adapt(int which, int ilo, int ihi, int jlo, int jhi,
                          double value)
{
  mu = value;
}
