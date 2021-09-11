/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */


#include "fix_PSM_LBM_BC.h"

fix_PSM_LBM_BC::fix_PSM_LBM_BC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  for(int ifix=0; ifix<modify->nfix; ifix++)
    if(strcmp(modify->fix[ifix]->style,"lbm-psm")==0)
      fixPSMLBM = (fix_PSM_LBM *)modify->fix[ifix];
}

fix_PSM_LBM_BC::~fix_PSM_LBM_BC()
{}

int fix_PSM_LBM_BC::setmask()
{
  int mask =0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= END_OF_STEP;
  return mask;
}


void fix_PSM_LBM_BC::init()
{
  zouHe2D = new ZouHeBC2D(fixPSMLBM->dynamics);
}


void fix_PSM_LBM_BC::pre_force(int)
{
  int envelopeWidth = 1;
  double u_infty = fixPSMLBM->unitConversion->get_u_lb();
  double rho_outlet = 1.0;

  // Shear
  zouHe2D->setZouHeVelBC2D_yn( 0+envelopeWidth,                               1+envelopeWidth, fixPSMLBM->dynamics->get_nx()-2-envelopeWidth, -u_infty );
  zouHe2D->setZouHeVelBC2D_yp( fixPSMLBM->dynamics->get_ny()-1-envelopeWidth, 1+envelopeWidth, fixPSMLBM->dynamics->get_nx()-2-envelopeWidth,  u_infty );

  // External flow
//  if (procCoordinates[0] == 0)
//    { boundaries.setZouHeVelBC2D_xn(  envelopeWidth,                      envelopeWidth, dynamics.get_ny()-1-envelopeWidth, u_infty    ); }
//  if (procCoordinates[0] == decomposition[0]-1)
//    { boundaries.setZouHeDensBC2D_xp( dynamics.get_nx()-1-envelopeWidth,  envelopeWidth, dynamics.get_ny()-1-envelopeWidth, rho_outlet ); }
//  zouHe2D->setZouHeVelBC2D_xn(  envelopeWidth,                                  envelopeWidth,    fixPSMLBM->dynamics->get_ny()-1-envelopeWidth,  u_infty    );
//  zouHe2D->setZouHeDensBC2D_xp( fixPSMLBM->dynamics->get_nx()-1-envelopeWidth,  envelopeWidth,    fixPSMLBM->dynamics->get_ny()-1-envelopeWidth,  rho_outlet );
//  zouHe2D->setZouHeVelBC2D_yn(  0+envelopeWidth,                                1+envelopeWidth,  fixPSMLBM->dynamics->get_nx()-2-envelopeWidth,  u_infty    );
//  zouHe2D->setZouHeVelBC2D_yp(  fixPSMLBM->dynamics->get_ny()-1-envelopeWidth,  1+envelopeWidth,  fixPSMLBM->dynamics->get_nx()-2-envelopeWidth,  u_infty    );
}
