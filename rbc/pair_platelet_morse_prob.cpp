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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "pair_platelet_morse_prob.h"
#include "update.h"
#include "random_mars.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairPlateletMorseProb::PairPlateletMorseProb(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;

  neveryFlag=0;
    
}

/* ---------------------------------------------------------------------- */

PairPlateletMorseProb::~PairPlateletMorseProb()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(d0);
    memory->destroy(alpha);
    memory->destroy(r0);
    memory->destroy(morse1);
    memory->destroy(offset);
    
    memory->destroy(kon);
    memory->destroy(koff);
    memory->destroy(r_rep);
    memory->destroy(nevery);
    memory->destroy(newType1);
    memory->destroy(newType2);
  }
    delete random;
}

/* ---------------------------------------------------------------------- */

void PairPlateletMorseProb::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,dr,dexp,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  // t0 used for activation time
  double *act = atom->dvector[0]; // activation
  double *t0 = atom->dvector[1];
  double *release = atom->dvector[2];
  double *actFunc = atom->dvector[3];
     
  //printf("lammps: get here\n");
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        double rand_num = random->uniform();
        
        if (rand_num < kon[itype][jtype]){// attraction probability
          r = sqrt(rsq);
          dr = r - r0[itype][jtype];
          dexp = exp(-alpha[itype][jtype] * dr);
          //fpair = factor_lj * morse1[itype][jtype] * (dexp*dexp - dexp) / r;
          fpair = factor_lj *sqrt(actFunc[i]*actFunc[j]) *morse1[itype][jtype] * (dexp*dexp - dexp) / r;

          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;
          if (newton_pair || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          }
          /*-------not needed for changing type ----------6/6/2019--jifu-------
          // change atom type
          if(nevery[itype][jtype]){
            if (update->ntimestep % nevery[itype][jtype] == 0){
              if ((type[i] != newType1[itype][jtype]) && (update->ntimestep - t0[i] >= nevery[itype][jtype])) {
                //printf("id %d type %d %d prob %g rand %g charge %g %g nevery %d \n",i, itype, jtype, prob[itype][jtype],rand_num, q[i], q[j], nevery[itype][jtype]);
                type[i] = newType1[itype][jtype];
                t0[i] =static_cast<double>(update->ntimestep);
              }
              if ((type[j] != newType2[itype][jtype]) && (update->ntimestep - t0[j] >= nevery[itype][jtype])) {
                //printf("id %d type %d %d prob %g rand %g charge %g %g nevery %d \n",j, itype, jtype, prob[itype][jtype],rand_num, q[i], q[j], nevery[itype][jtype]);
                type[j]=newType2[itype][jtype];
                t0[j] =static_cast<double>(update->ntimestep);
              }
            }
          }//eof atom type
        */ 

          if (eflag) {
            evdwl = d0[itype][jtype] * (dexp*dexp - 2.0*dexp) -
              offset[itype][jtype];
            evdwl *= factor_lj;
          }

          if (evflag) ev_tally(i,j,nlocal,newton_pair,
                               evdwl,0.0,fpair,delx,dely,delz);
        }else if(rsq < r_rep[itype][jtype]){// pure repulsive
          r = sqrt(rsq);
          dr = r - r0[itype][jtype];
          dexp = exp(-alpha[itype][jtype] * dr);
          fpair = factor_lj *sqrt(actFunc[i]*actFunc[j])* morse1[itype][jtype] * (dexp*dexp - dexp) / r;
          //fpair = factor_lj * morse1[itype][jtype] * (dexp*dexp - dexp) / r;

          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;
          if (newton_pair || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          }
          
          // change atom type
          //if(nevery[itype][jtype]){
          //  if (update->ntimestep % nevery[itype][jtype] == 0){
          //    type[i]=newType1[itype][jtype];
          //    type[j]=newType2[itype][jtype];
          //  }
          //}//eof atom type
       
          if (eflag) {
            evdwl = d0[itype][jtype] * (dexp*dexp - 2.0*dexp) -
              offset[itype][jtype];
            evdwl *= factor_lj;
          }

          if (evflag) ev_tally(i,j,nlocal,newton_pair,
                               evdwl,0.0,fpair,delx,dely,delz);
        }
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairPlateletMorseProb::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(d0,n+1,n+1,"pair:d0");
  memory->create(alpha,n+1,n+1,"pair:alpha");
  memory->create(r0,n+1,n+1,"pair:r0");
  memory->create(morse1,n+1,n+1,"pair:morse1");
  memory->create(offset,n+1,n+1,"pair:offset");
  
  memory->create(r_rep,n+1,n+1,"pair:repulsive");
  memory->create(kon,n+1,n+1,"pair:kon");
  memory->create(koff,n+1,n+1,"pair:koff");
  memory->create(nevery,n+1,n+1,"pair:nevery");
  memory->create(newType1,n+1,n+1,"pair:itype");
  memory->create(newType2,n+1,n+1,"pair:jtype");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPlateletMorseProb::settings(int narg, char **arg)
{
  if (narg > 2) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  seed = force->inumeric(FLERR,arg[1]);
  MPI_Comm_rank(world, &me);
  //seed=12345;
  random = new RanMars(lmp,seed+me);
  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPlateletMorseProb::coeff(int narg, char **arg)
{
  if (narg < 7) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double d0_one = force->numeric(FLERR,arg[2]);
  double alpha_one = force->numeric(FLERR,arg[3]);
  double r0_one = force->numeric(FLERR,arg[4]);

  double cut_one = cut_global;
  cut_one = force->numeric(FLERR,arg[5]);

  double rep_one = force->numeric(FLERR,arg[6]);
  double kon_one = 1.0;//force->numeric(FLERR,arg[7]);
  double koff_one = 0.0;//force->numeric(FLERR,arg[7]);

  int newType1_one=0;
  int newType2_one=0;
  int nevery_one=0;

  int iarg=7;
  while (iarg<narg){
    if (strcmp(arg[iarg],"nevery")==0){
      if (iarg+4>narg) error->all(FLERR,"Illegal pair_style platelet/morse/prob command");
      neveryFlag=1;
      nevery_one = force->inumeric(FLERR,arg[iarg+1]);
      newType1_one = force->inumeric(FLERR, arg[iarg+2]);
      newType2_one = force->inumeric(FLERR, arg[iarg+3]);
      iarg += 4;
    }else if(strcmp(arg[iarg],"kon")==0){
      if (iarg+2>narg) error->all(FLERR,"Illegal pair_style platelet/morse/prob command");
      probFlag=1;
      kon_one = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    }else if(strcmp(arg[iarg],"koff")==0){
      if (iarg+2>narg) error->all(FLERR,"Illegal pair_style platelet/morse/prob command");
      probFlag=1;
      koff_one = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    }

  }
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      d0[i][j] = d0_one;
      alpha[i][j] = alpha_one;
      r0[i][j] = r0_one;
      cut[i][j] = cut_one;
      r_rep[i][j] = rep_one;
      kon[i][j] = kon_one;
      koff[i][j] = koff_one;
      nevery[i][j] = nevery_one;
      newType1[i][j] = newType1_one;
      newType2[i][j] = newType2_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPlateletMorseProb::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  morse1[i][j] = 2.0*d0[i][j]*alpha[i][j];

  if (offset_flag) {
    double alpha_dr = -alpha[i][j] * (cut[i][j] - r0[i][j]);
    offset[i][j] = d0[i][j] * (exp(2.0*alpha_dr) - 2.0*exp(alpha_dr));
  } else offset[i][j] = 0.0;

  d0[j][i] = d0[i][j];
  alpha[j][i] = alpha[i][j];
  r0[j][i] = r0[i][j];
  morse1[j][i] = morse1[i][j];
  offset[j][i] = offset[i][j];
  
  r_rep[j][i] = r_rep[i][j];
  kon[j][i] = kon[i][j];
  koff[j][i] = koff[i][j];
  nevery[j][i] =nevery[i][j];
  newType1[j][i] = newType2[i][j];
  newType2[j][i] = newType1[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairPlateletMorseProb::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&d0[i][j],sizeof(double),1,fp);
        fwrite(&alpha[i][j],sizeof(double),1,fp);
        fwrite(&r0[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
        fwrite(&r_rep[i][j],sizeof(double),1,fp);
        fwrite(&kon[i][j],sizeof(double),1,fp);
        fwrite(&koff[i][j],sizeof(double),1,fp);
        fwrite(&nevery[i][j],sizeof(int),1,fp);
        fwrite(&newType1[i][j],sizeof(int),1,fp);
        fwrite(&newType2[i][j],sizeof(int),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairPlateletMorseProb::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&d0[i][j],sizeof(double),1,fp);
          fread(&alpha[i][j],sizeof(double),1,fp);
          fread(&r0[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
          fread(&r_rep[i][j],sizeof(double),1,fp);
          fread(&kon[i][j],sizeof(double),1,fp);
          fread(&koff[i][j],sizeof(double),1,fp);
          fread(&nevery[i][j],sizeof(int),1,fp);
          fread(&newType1[i][j],sizeof(int),1,fp);
          fread(&newType2[i][j],sizeof(int),1,fp);
        }
        MPI_Bcast(&d0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&r_rep[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&kon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&koff[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&nevery[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&newType1[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&newType2[i][j],1,MPI_INT,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairPlateletMorseProb::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairPlateletMorseProb::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairPlateletMorseProb::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,d0[i][i],alpha[i][i],r0[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairPlateletMorseProb::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g\n",
              i,j,d0[i][j],alpha[i][j],r0[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairPlateletMorseProb::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double r,dr,dexp,phi;

  r = sqrt(rsq);
  dr = r - r0[itype][jtype];
  dexp = exp(-alpha[itype][jtype] * dr);
  fforce = factor_lj * morse1[itype][jtype] * (dexp*dexp - dexp) / r;

  phi = d0[itype][jtype] * (dexp*dexp - 2.0*dexp) - offset[itype][jtype];
  return factor_lj*phi;
}

/* ---------------------------------------------------------------------- */

void *PairPlateletMorseProb::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"d0") == 0) return (void *) d0;
  if (strcmp(str,"r0") == 0) return (void *) r0;
  if (strcmp(str,"alpha") == 0) return (void *) alpha;
  return NULL;
}
