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

#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "fix_activate_platelet.h"
#include "atom.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "math.h"
using namespace LAMMPS_NS;
using namespace FixConst;

//enum{PF_CALLBACK,PF_ARRAY};

/* ---------------------------------------------------------------------- */

FixActivatePlatelet::FixActivatePlatelet(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 8) error->all(FLERR,"Illegal fix ActivatePlatelet command");
  
    // perform initial allocation of atom-based array
  // register with Atom class
  //napply = force->inumeric(FLERR,arg[3]);
  //rate = force->inumeric(FLERR,arg[3]);
  xi0 = utils::numeric(FLERR,arg[3],false,lmp);
  m0 = utils::numeric(FLERR,arg[4],false,lmp);
  tau = utils::numeric(FLERR,arg[5],false,lmp);
  alpha = utils::numeric(FLERR,arg[6],false,lmp);
  n = utils::numeric(FLERR,arg[7],false,lmp);
  
  fexternal = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  
}

/* ---------------------------------------------------------------------- */

FixActivatePlatelet::~FixActivatePlatelet()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);
  memory->destroy(fexternal);
}

/* ---------------------------------------------------------------------- */

int FixActivatePlatelet::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixActivatePlatelet::init()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
    /*double * act=atom->dvector[1];
    double * tst=atom->dvector[2];
    double * pF=atom->dvector[3];
    std::cout<<" nlocal "<<nlocal<<std::endl;
    int tmp=1;
    int index=atom->find_custom("Tst",tmp);
    std::cout<<" index "<<index<<" tst "<<tst[0]<<" "<<tst[1]<< " "<< tst[2]<<std::endl;
    
        std::cout<<"act "<<atom->dvector[0][0]<<" tst "<<atom->dvector[0][1]<<" F "<<atom->dvector[0][2]<<std::endl;
        std::cout<<"act "<<atom->dvector[1][0]<<" tst "<<atom->dvector[1][1]<<" F "<<atom->dvector[1][2]<<std::endl;
        std::cout<<"act "<<atom->dvector[2][0]<<" tst "<<atom->dvector[2][1]<<" F "<<atom->dvector[2][2]<<std::endl;*/
  for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        fexternal[i][0]=0.;
        fexternal[i][1]=0.;
        fexternal[i][2]=0.;
        //std::cout<<"act "<<act[i]<<" tst "<<tst[i]<<" F "<<pF[i]<<std::endl;
        //std::cout<<"act "<<atom->dvector[1][0]<<" tst "<<atom->dvector[1][1]<<" F "<<atom->dvector[1][2]<<std::endl;
        
      }
}

/* ---------------------------------------------------------------------- */

void FixActivatePlatelet::end_of_step()
{/*
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        fexternal[i][0]=0.;
        fexternal[i][1]=0.;
        fexternal[i][2]=0.;
      }*/
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    double * act=atom->dvector[0];
    double * tst=atom->dvector[1];
    double * release=atom->dvector[2];
    double * actFunc=atom->dvector[3];
    double tmp;
    bigint ntimestep = update->ntimestep;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if ((tst[i]<0) && (act[i]>= xi0))
            tst[i]=ntimestep;
        
        if ((ntimestep > tst[i] ) && (tst[i] >0)){
            release[i]=m0/tau*exp(-(ntimestep-tst[i])/tau);
            tmp = pow(xi0/act[i],n);
            actFunc[i]=alpha+(1.-alpha)/(1+tmp);
        }
      }
    //int tmp=1;
    //int index=atom->find_custom("Tst",tmp);
}
/* ---------------------------------------------------------------------- */

void FixActivatePlatelet::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixActivatePlatelet::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixActivatePlatelet::post_force(int vflag)
{
  bigint ntimestep = update->ntimestep;

  // invoke the callback in driver program
  // it will fill fexternal with forces


  // add forces from fexternal to atoms in group
/*
  if (ntimestep % napply == 0) {
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        f[i][0] += fexternal[i][0];
        f[i][1] += fexternal[i][1];
        f[i][2] += fexternal[i][2];
      }
  }*/
}

/* ---------------------------------------------------------------------- */

void FixActivatePlatelet::min_post_force(int vflag)
{
  post_force(vflag);
}



/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixActivatePlatelet::memory_usage()
{
  double bytes = 3*atom->nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixActivatePlatelet::grow_arrays(int nmax)
{
  memory->grow(fexternal,nmax,3,"rls:fexternall");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixActivatePlatelet::copy_arrays(int i, int j, int delflag)
{
  fexternal[j][0] = fexternal[i][0];
  fexternal[j][1] = fexternal[i][1];
  fexternal[j][2] = fexternal[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixActivatePlatelet::pack_exchange(int i, double *buf)
{
  buf[0] = fexternal[i][0];
  buf[1] = fexternal[i][1];
  buf[2] = fexternal[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixActivatePlatelet::unpack_exchange(int nlocal, double *buf)
{
  fexternal[nlocal][0] = buf[0];
  fexternal[nlocal][1] = buf[1];
  fexternal[nlocal][2] = buf[2];
  return 3;
}


