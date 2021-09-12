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
  if (narg < 4) error->all(FLERR,"Illegal fix lbm-psm-bc command");
  
   typeBC = 0;

  int iarg = 3;
  if (strcmp(arg[iarg],"shear") == 0) {
    if (iarg+1 > narg) error->all(FLERR,"Illegal fix lbm-psm-bc command");
    typeBC = 1;
  } else if (strcmp(arg[iarg],"xFlow") == 0) {
    if (iarg+1 > narg) error->all(FLERR,"Illegal fix lbm-psm-bc command");
    typeBC = 2;
  } else error->all(FLERR,"Illegal fix lbm-psm-bc command");


/* TODO Use when genral BC are implemented
  for (int i = 0; i < 3; i++)
  {
    lowerBC[i] = 0;
    upperBC[i] = 0;
    
    lowerRhoBC[i] = -1.0;
    lowerVelBC[i] = -1.0;
    upperRhoBC[i] = -1.0;
    upperVelBC[i] = -1.0;
  }

  

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix lbm-psm-bc command");

      if (strcmp(arg[iarg+1],"rho") == 0) {
        lowerRhoBC[0] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        lowerBC[0] = 2;
      }
      else if ((strcmp(arg[iarg+1],"vel") == 0) {
        lowerVelBC[0] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        lowerBC[0] = 1;
      } else error->all(FLERR,"Illegal fix lbm-psm-bc command");

      if (strcmp(arg[iarg+3],"rho") == 0) {
        upperRhoBC_x = utils::numeric(FLERR,arg[iarg+4],false,lmp);
        upperBC[0] = 2;
      }
      else if ((strcmp(arg[iarg+3],"vel") == 0) {
        upperVelBC[0] = utils::numeric(FLERR,arg[iarg+4],false,lmp);
        upperBC[0] = 1;
      } else error->all(FLERR,"Illegal fix lbm-psm-bc command");

    } else if (strcmp(arg[iarg],"y") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix lbm-psm-bc command");

      if (strcmp(arg[iarg+1],"rho") == 0) {
        lowerRhoBC[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        lowerBC[1] = 2;
      }
      else if ((strcmp(arg[iarg+1],"vel") == 0) {
        lowerVelBC[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        lowerBC[1] = 1;
      } else error->all(FLERR,"Illegal fix lbm-psm-bc command");

      if (strcmp(arg[iarg+3],"rho") == 0) {
        upperRhoBC[1] = utils::numeric(FLERR,arg[iarg+4],false,lmp);
        upperBC[1] = 2;
      }
      else if ((strcmp(arg[iarg+3],"vel") == 0) {
        upperVelBC[1] = utils::numeric(FLERR,arg[iarg+4],false,lmp);
        upperBC[1] = 1;
      } else error->all(FLERR,"Illegal fix lbm-psm-bc command");

    } else if (strcmp(arg[iarg],"z") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix lbm-psm-bc command");

      if (strcmp(arg[iarg+1],"rho") == 0) {
        lowerRhoBC[2] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        lowerBC[2] = 2;
      }
      else if ((strcmp(arg[iarg+1],"vel") == 0) {
        lowerVelBC[2] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        lowerBC[0] = 1;
      } else error->all(FLERR,"Illegal fix lbm-psm-bc command");

      if (strcmp(arg[iarg+3],"rho") == 0) {
        upperRhoBC[2] = utils::numeric(FLERR,arg[iarg+4],false,lmp);
        upperBC[2] = 2;
      }
      else if ((strcmp(arg[iarg+3],"vel") == 0) {
        upperVelBC[2] = utils::numeric(FLERR,arg[iarg+4],false,lmp);
        upperBC[2] = 1;
      } else error->all(FLERR,"Illegal fix lbm-psm-bc command");

    } else error->all(FLERR,"Illegal fix lbm-psm-bc command");
  }
*/

  for(int ifix=0; ifix<modify->nfix; ifix++)
    if(strcmp(modify->fix[ifix]->style,"lbm-psm")==0)
      fixPSMLBM = (fix_PSM_LBM *)modify->fix[ifix];
}


fix_PSM_LBM_BC::~fix_PSM_LBM_BC()
{
//  delete zouHe2D;
//  delete fixPSMLBM;
}


int fix_PSM_LBM_BC::setmask()
{
  int mask = 0;
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
  if (update->ntimestep % fixPSMLBM->nevery) return;

  int envelopeWidth = fixPSMLBM->dynamics->get_envelopeWidth();
  double u_infty = fixPSMLBM->unitConversion->get_u_lb();
  double rho_outlet = 1.0;


/* // TODO General formulation of BC
  /// BC in x-direction ///
  if (lowerBC[0] != 0 && comm->myloc[0] == 0)
  {
    if (lowerVelBC[0] >= 0.0)
    {
      zouHe2D->setZouHeVelBC2D_xn( envelopeWidth, envelopeWidth, fixPSMLBM->dynamics->get_ny()-1-envelopeWidth, u_infty );
    }
// Not implemented yet
//    if (lowerRhoBC[0] >= 0.0)
//    {
//      zouHe2D->setZouHeDensBC2D_xn( envelopeWidth, envelopeWidth, fixPSMLBM->dynamics->get_ny()-1-envelopeWidth,  rho_outlet );
//    }
  }

  if (upperBC[0] != 0 && comm->myloc[0] == comm->procgrid[0]-1)
  {
// Not implemented yet
//    if (upperVelBC[0] >= 0.0)
//    {
//      zouHe2D->setZouHeVelBC2D_xn( envelopeWidth, envelopeWidth, fixPSMLBM->dynamics->get_ny()-1-envelopeWidth, u_infty );
//    }
    if (upperRhoBC[0] >= 0.0)
    {
      zouHe2D->setZouHeDensBC2D_xp( fixPSMLBM->dynamics->get_nx()-1-envelopeWidth, envelopeWidth, fixPSMLBM->dynamics->get_ny()-1-envelopeWidth, rho_outlet );
    }
  }


  /// BC in y-direction ///
  if (lowerBC[1] != 0 && comm->myloc[1] == 0)
  {
    if (lowerVelBC[0] >= 0.0)
    {
      zouHe2D->setZouHeVelBC2D_yn( 0+envelopeWidth, 1+envelopeWidth, fixPSMLBM->dynamics->get_nx()-2-envelopeWidth, u_infty );
    }
// Not implemented yet
//    if (lowerRhoBC[0] >= 0.0)
//    {
//      zouHe2D->setZouHeDensBC2D_yn( );
//    }
  }

  if (upperBC[1] != 0 && comm->myloc[1] == comm->procgrid[1]-1)
  {
    if (upperVelBC[0] >= 0.0)
    {
      zouHe2D->setZouHeVelBC2D_yp( fixPSMLBM->dynamics->get_ny()-1-envelopeWidth, 1+envelopeWidth, fixPSMLBM->dynamics->get_nx()-2-envelopeWidth, u_infty );
    }
// Not implemented yet
//    if (upperRhoBC[0] >= 0.0)
//    {
//      zouHe2D->setZouHeDensBC2D_xp( fixPSMLBM->dynamics->get_nx()-1-envelopeWidth, envelopeWidth, fixPSMLBM->dynamics->get_ny()-1-envelopeWidth, rho_outlet );
//    }
  }


  /// BC in z-direction ///
  // None implemented yet, i.e. periodic (choose periodic in z-direction in input script

*/


  // Shear
  if (typeBC == 1)
  {
    if (comm->myloc[1] == 0)
      { zouHe2D->setZouHeVelBC2D_yn( 0+envelopeWidth, 1+envelopeWidth, fixPSMLBM->dynamics->get_nx()-2-envelopeWidth, -u_infty ); }
    if (comm->myloc[1] == comm->procgrid[1]-1)
      { zouHe2D->setZouHeVelBC2D_yp( fixPSMLBM->dynamics->get_ny()-1-envelopeWidth, 1+envelopeWidth, fixPSMLBM->dynamics->get_nx()-2-envelopeWidth, u_infty ); }
  }
  else if (typeBC == 2)
  {
    if (comm->myloc[0] == 0)
      { zouHe2D->setZouHeVelBC2D_xn( envelopeWidth, envelopeWidth, fixPSMLBM->dynamics->get_ny()-1-envelopeWidth, u_infty ); }
    if (comm->myloc[0] == comm->procgrid[0]-1)
      { zouHe2D->setZouHeDensBC2D_xp( fixPSMLBM->dynamics->get_nx()-1-envelopeWidth, envelopeWidth, fixPSMLBM->dynamics->get_ny()-1-envelopeWidth, rho_outlet ); }
    if (comm->myloc[1] == 0)
      { zouHe2D->setZouHeVelBC2D_yn( 0+envelopeWidth, 1+envelopeWidth, fixPSMLBM->dynamics->get_nx()-2-envelopeWidth, u_infty ); }
    if (comm->myloc[1] == comm->procgrid[1]-1)
      { zouHe2D->setZouHeVelBC2D_yp( fixPSMLBM->dynamics->get_ny()-1-envelopeWidth, 1+envelopeWidth, fixPSMLBM->dynamics->get_nx()-2-envelopeWidth, u_infty ); }
  }

  // External flow
//  if (procCoordinates[0] == 0)
//    { boundaries.setZouHeVelBC2D_xn(  envelopeWidth,                      envelopeWidth, dynamics.get_ny()-1-envelopeWidth, u_infty    ); }
//  if (procCoordinates[0] == decomposition[0]-1)
//    { boundaries.setZouHeDensBC2D_xp( dynamics.get_nx()-1-envelopeWidth,  envelopeWidth, dynamics.get_ny()-1-envelopeWidth, rho_outlet ); }
}
