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
#include "stdlib.h"
#include "bond_wlc_pow.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondWLCPOW::BondWLCPOW(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondWLCPOW::~BondWLCPOW()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(r0);
    memory->destroy(a0);
    memory->destroy(kp);
  }
}

/* ---------------------------------------------------------------------- */

void BondWLCPOW::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,lm,a,a2,a3,a4;

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    
    lm = r0[type]/a0[type];
    a = r/lm;//ratio
    a2 = a*a;
    a3 = a2*a;
    a4 = (1-a)*(1-a);

    // force & energy

    if (r > 0.0) fbond = (-k[type]*(0.25/a4 - 0.25 + a) + kp[type]/rsq)/r;
    else fbond = 0.0;
    
    if (eflag) ebond = 0.25*k[type]*lm*(3*a2 - 2*a3)/(1 - a) + kp[type]/r;

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx*fbond;
      f[i1][1] += dely*fbond;
      f[i1][2] += delz*fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx*fbond;
      f[i2][1] -= dely*fbond;
      f[i2][2] -= delz*fbond;
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondWLCPOW::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(r0,n+1,"bond:r0");
  memory->create(a0,n+1,"bond:a0");
  memory->create(kp,n+1,"bond:kp");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondWLCPOW::coeff(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  //force->bounds(FLERR,arg[0],atom->nbondtypes,ilo,ihi);
  utils::bounds(FLERR,arg[0],1,atom->nbondtypes,ilo,ihi,error);

  double k_one = utils::numeric(FLERR,arg[1],false,lmp);
  double r0_one = utils::numeric(FLERR,arg[2],false,lmp);
  double a0_one = utils::numeric(FLERR,arg[3],false,lmp);
  double kp_one = utils::numeric(FLERR,arg[4],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    a0[i] = a0_one;
    kp[i] = kp_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondWLCPOW::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondWLCPOW::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&a0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&kp[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondWLCPOW::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&a0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&kp[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&a0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kp[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondWLCPOW::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,k[i],r0[i],a0[i],kp[i]);
}

/* ---------------------------------------------------------------------- */

double BondWLCPOW::single(int type, double rsq, int i, int j,
                        double &fforce)
{
  double r = sqrt(rsq);
  double lm = r0[type]/a0[type];
  double a = r/lm;//ratio
  double a2 = a*a;
  double a3 = a2*a;
  double a4 = (1-a)*(1-a);

  fforce = 0;
  if (r > 0.0) fforce = (-k[type]*(0.25/a4 - 0.25 + a) + kp[type]/rsq)/r;
  return 0.25*k[type]*lm*(3*a2 - 2*a3)/(1 - a) + kp[type]/r;
}
