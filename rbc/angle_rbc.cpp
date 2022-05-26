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

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "angle_rbc.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "update.h"//debug only
#include <vector>
#include <set>

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleRbc::AngleRbc(LAMMPS *lmp) : Angle(lmp) {
  At = Vt = NULL;
  at = vt = NULL;
  crossFlag = NULL;
  MAXxyz = maxxyz = MINxyz = minxyz =  NULL;
}

/* ---------------------------------------------------------------------- */

AngleRbc::~AngleRbc()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(Cq);
    memory->destroy(q);
    memory->destroy(ka);
    memory->destroy(A0t);
    memory->destroy(kv);
    memory->destroy(V0t);
    memory->destroy(kd);
    memory->destroy(A0);
  }
  if(nmolecule > 0 && !copymode){
    memory->destroy(At);
    memory->destroy(Vt);
    memory->destroy(at);
    memory->destroy(vt);
    memory->destroy(crossFlag);
    memory->destroy(MAXxyz);
    memory->destroy(maxxyz);
    memory->destroy(MINxyz);
    memory->destroy(minxyz);
  }
}

/* ---------------------------------------------------------------------- */
void AngleRbc::init_style()
{
  int i,n;

  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;
  
  n = 0;
  nmolecule = 0;
  for (i = 0; i < nlocal; i++){
    n = MAX(n,molecule[i]); 
  }
  MPI_Allreduce(&n,&nmolecule,1,MPI_INT,MPI_MAX,world);
  //printf("max # of molecules: %d\n", nmolecule);
  
  memory->create(At,nmolecule+1,"angle:At");
  memory->create(Vt,nmolecule+1,"angle:Vt");
  memory->create(at,nmolecule+1,"angle:At");
  memory->create(vt,nmolecule+1,"angle:Vt");
  memory->create(crossFlag,nmolecule+1,3,"angle:crossFlag");
  memory->create(MAXxyz,nmolecule+1,3,"angle:MAXxyz");
  memory->create(maxxyz,nmolecule+1,3,"angle:maxxyz");
  memory->create(MINxyz,nmolecule+1,3,"angle:MINxyz");
  memory->create(minxyz,nmolecule+1,3,"angle:minxyz");
  check_crossing(crossFlag);
  for (int j=1;j<nmolecule+1;j++)
  { if (((crossFlag[j][0]) && (domain->xperiodic)) || ((crossFlag[j][1]) && (domain->yperiodic))|| ((crossFlag[j][2]) && (domain->zperiodic)))
      if (comm->me ==0) error->warning(FLERR,"Cell is straddled across the boundary or the cell size is bigger than half of the box side, check your initial input files\n");
  }
}

/* ---------------------------------------------------------------------- */
void AngleRbc::computeAreaVol(double *At, double *Vt, int **crossFlag){
  for (int i = 1; i < nmolecule+1; i++) {
    At[i] = 0.0;
    Vt[i] = 0.0;
    crossFlag[i][0]=0;
    crossFlag[i][1]=0;
    crossFlag[i][2]=0;
  }
  int i1,i2,i3,n,type,molId;
  double delx1,dely1,delz1,delx2,dely2,delz2,delx3,dely3,delz3;
  double xi[3],cnt[3],xi2;
  double at,vt;
  double nhat,area,a_xy;
  double rsq1,rsq2,r1,r2;
    
  double x1[3],x2[3],x3[3];
  double **x = atom->x;
  double **f = atom->f;
  tagint *molecule = atom->molecule;
  tagint *image = atom->image;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
    
  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;
  check_crossing(crossFlag);
  // each processor has nanglelist, see neighbor.cpp
  // calculate area and volume
  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];
    molId = molecule[i2];

    //if (molId < 0 || molId > nmolecule+1)  error->all(FLERR,"Wrong molecule ID");

    // 1st bond
    /*    /\3
     *   /  \
     *  /1___\2
     * dx_ij=x_i - x_j
     * */
    /*delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];*/
    /*domain->unmap(x[i1],image[i1],x1);
    domain->unmap(x[i2],image[i2],x2);
    domain->unmap(x[i3],image[i3],x3);*/
    positionShift(x[i1],x1,crossFlag[molId]); 
    positionShift(x[i2],x2,crossFlag[molId]); 
    positionShift(x[i3],x3,crossFlag[molId]); 
    delx1 = x1[0] - x2[0];
    dely1 = x1[1] - x2[1];
    delz1 = x1[2] - x2[2];
    
    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond
    /*delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];*/
    delx2 = x3[0] - x2[0];
    dely2 = x3[1] - x2[1];
    delz2 = x3[2] - x2[2];
    
    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // 3rd bond
    /*delx3 = x[i1][0] - x[i3][0];
    dely3 = x[i1][1] - x[i3][1];
    delz3 = x[i1][2] - x[i3][2];*/
    delx3 = x1[0] - x3[0];
    dely3 = x1[1] - x3[1];
    delz3 = x1[2] - x3[2];
    
    // norm xi=dx_2 x dx_1
    xi[0]=delz1*dely2 - dely1*delz2;
    xi[1]=delx1*delz2 - delz1*delx2;
    xi[2]=dely1*delx2 - delx1*dely2;
    
    xi2=xi[0]*xi[0] + xi[1]*xi[1] + xi[2]*xi[2];
    area = 0.5*sqrt(xi2);
    //At[molId] += 0.5*sqrt(xi2);
    At[molId] += area;

    //new code
    /*if (xperiodic || yperiodic){
      a_xy = 0.5*xi[2];//area projection on xy plane
      Vt[molId] += a_xy*(x[i1][2] + x[i2][2] + x[i3][2])/3.0;
    }else if(zperiodic){
      a_xy = 0.5*xi[1];//area projection on xz plane
      Vt[molId] += a_xy*(x[i1][1] + x[i2][1] + x[i3][1])/3.0;
    }*/


    /*cnt[0]=(x[i1][0] + x[i2][0] + x[i3][0])/3.0;
    cnt[1]=(x[i1][1] + x[i2][1] + x[i3][1])/3.0;
    cnt[2]=(x[i1][2] + x[i2][2] + x[i3][2])/3.0;*/
    cnt[0]=(x1[0] + x2[0] + x3[0])/3.0;
    cnt[1]=(x1[1] + x2[1] + x3[1])/3.0;
    cnt[2]=(x1[2] + x2[2] + x3[2])/3.0;
    
    Vt[molId] += 1.0/6.0*(xi[0]*cnt[0] + xi[1]*cnt[1] + xi[2]*cnt[2]);

   /* 
    bigint ntimestep;
    ntimestep = update->ntimestep;
    if (ntimestep > 1700 && n==13){
                printf("angle %d\n",n);
                printf("crossing %d %d %d\n",crossFlag[0],crossFlag[1],crossFlag[2]);
                printf("Vt %lg V0t %lg At %lg A0t %lg\n",Vt[molId],V0t[type],At[molId],A0t[type]);
                printf("dV %lg\n",1.0/6.0*(xi[0]*cnt[0] + xi[1]*cnt[1] + xi[2]*cnt[2]) );
                printf("dA %lg\n",0.5*sqrt(xi2) );
                printf("image %d %d %d\n",image[i1],image[i2],image[i3] );
                
                printf("cnt %lg %lg %lg xi %lg %lg %lg\n",cnt[0],cnt[1],cnt[2],xi[0],xi[1],xi[2] );
                //printf("x1 %lg %lg %lg x2 %lg %lg %lg x3 %lg %lg %lg\n",x[i1][0],x[i1][1],x[i1][2],x[i2][0],x[i2][1],x[i2][2],x[i3][0],x[i3][1],x[i3][2] );
                printf("x1 %lg %lg %lg x2 %lg %lg %lg x3 %lg %lg %lg\n",x1[0],x1[1],x1[2],x2[0],x2[1],x2[2],x3[0],x3[1],x3[2] );
                printf("dx1 %lg %lg %lg dx2 %lg %lg %lg dx3 %lg %lg %lg\n",delx1,dely1,delz1,delx2,dely2,delz2,delx3,dely3,delz3);
            }*/
  }
  for (int i = 1; i < nmolecule+1; i++) {
    MPI_Allreduce(&At[i],&at,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&Vt[i],&vt,1,MPI_DOUBLE,MPI_SUM,world);
    At[i] = at;
    Vt[i] = vt;
    //if (comm->me == 0) 
    //printf("%d cell: Vt %lg  At %lg \n",i,Vt[i],At[i]);
  }
}

/* ------------------No check on crossing boundary-------------------- */
void AngleRbc::computeAreaVol(double *At, double *Vt){
  for (int i = 1; i < nmolecule+1; i++) {
    At[i] = 0.0;
    Vt[i] = 0.0;
    at[i] = 0.0;
    vt[i] = 0.0;
  }
  int i1,i2,i3,n,type,molId;
  double delx1,dely1,delz1,delx2,dely2,delz2,delx3,dely3,delz3;
  double xi[3],cnt[3],xi2;
  //double at,vt;
  double nhat,area,a_xy;
  double rsq1,rsq2,r1,r2;
    
  double x1[3],x2[3],x3[3];
  double **x = atom->x;
  double **f = atom->f;
  tagint *molecule = atom->molecule;
  tagint *image = atom->image;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
    
  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;
  
  if (xperiodic && yperiodic && zperiodic) error->all(FLERR,"At least one nonperiodic boundary condition for volume calculation");
    //printf("%d processor:  %d angles  At %lg \n",comm->me,nanglelist);
  // each processor has nanglelist, see neighbor.cpp
  // calculate area and volume
  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];
    molId = molecule[i2];

    //if (molId < 0 || molId > nmolecule+1)  error->all(FLERR,"Wrong molecule ID");

    // 1st bond
    /*    /\3
     *   /  \
     *  /1___\2
     * dx_ij=x_i - x_j
     * */
    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];
        
    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond
    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];
    
    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // 3rd bond
    delx3 = x[i1][0] - x[i3][0];
    dely3 = x[i1][1] - x[i3][1];
    delz3 = x[i1][2] - x[i3][2];
    
    // norm xi=dx_2 x dx_1
    xi[0]=delz1*dely2 - dely1*delz2;
    xi[1]=delx1*delz2 - delz1*delx2;
    xi[2]=dely1*delx2 - delx1*dely2;
    
    xi2=xi[0]*xi[0] + xi[1]*xi[1] + xi[2]*xi[2];
    area = 0.5*sqrt(xi2);
    //At[molId] += 0.5*sqrt(xi2);
    //At[molId] += area;
    at[molId] += area;

    //new code
    if (!zperiodic){
      a_xy = 0.5*xi[2];//area projection on xy plane
      //Vt[molId] += a_xy*(x[i1][2] + x[i2][2] + x[i3][2])/3.0;
      vt[molId] += a_xy*(x[i1][2] + x[i2][2] + x[i3][2])*0.333333333;
    }else if(!yperiodic){
      a_xy = 0.5*xi[1];//area projection on xz plane
      vt[molId] += a_xy*(x[i1][1] + x[i2][1] + x[i3][1])*0.333333333;
      //cnt[1]=(x[i1][1] + x[i2][1] + x[i3][1])/3.0;
      //printf("%d cell:%d angle, area:  %lg cnt %lg  vt %lg \n",molId,n,a_xy,cnt[1],vt[molId]);// V0t[type], not i   
    }else if(!xperiodic){
      a_xy = 0.5*xi[0];//area projection on xz plane
      vt[molId] += a_xy*(x[i1][0] + x[i2][0] + x[i3][0])*0.333333333;
    }



   /* 
    bigint ntimestep;
    ntimestep = update->ntimestep;
    if (ntimestep > 1700 && n==13){
                printf("angle %d\n",n);
                printf("crossing %d %d %d\n",crossFlag[0],crossFlag[1],crossFlag[2]);
                printf("Vt %lg V0t %lg At %lg A0t %lg\n",Vt[molId],V0t[type],At[molId],A0t[type]);
                printf("dV %lg\n",1.0/6.0*(xi[0]*cnt[0] + xi[1]*cnt[1] + xi[2]*cnt[2]) );
                printf("dA %lg\n",0.5*sqrt(xi2) );
                printf("image %d %d %d\n",image[i1],image[i2],image[i3] );
                
                printf("cnt %lg %lg %lg xi %lg %lg %lg\n",cnt[0],cnt[1],cnt[2],xi[0],xi[1],xi[2] );
                //printf("x1 %lg %lg %lg x2 %lg %lg %lg x3 %lg %lg %lg\n",x[i1][0],x[i1][1],x[i1][2],x[i2][0],x[i2][1],x[i2][2],x[i3][0],x[i3][1],x[i3][2] );
                printf("x1 %lg %lg %lg x2 %lg %lg %lg x3 %lg %lg %lg\n",x1[0],x1[1],x1[2],x2[0],x2[1],x2[2],x3[0],x3[1],x3[2] );
                printf("dx1 %lg %lg %lg dx2 %lg %lg %lg dx3 %lg %lg %lg\n",delx1,dely1,delz1,delx2,dely2,delz2,delx3,dely3,delz3);
            }*/
  }
    
    MPI_Allreduce(at,At,nmolecule+1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(vt,Vt,nmolecule+1,MPI_DOUBLE,MPI_SUM,world);
    /*if (comm->me == 0){ 
    for (int i = 1; i < nmolecule+1; i++) {
            printf("%d cell: Vt %lg  At %lg \n",i,Vt[i],At[i]);
        }
    }*/
}
/* ---------------------------------------------------------------------- */

void AngleRbc::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type,molId;
  double delx1,dely1,delz1,delx2,dely2,delz2,delx3,dely3,delz3;
  double eangle,f1[3],f2[3],f3[3];
  double fa1[3],fa2[3],fa3[3],fv1[3],fv2[3],fv3[3];//debug
  double dA;
  double xi[3],cnt[3],xi2;//
  double a0,beta_a, beta_ad, beta_area, beta_v;

  eangle = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;
  
  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  tagint *molecule = atom->molecule;
  tagint *image = atom->image;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  tagint g1,g2,g3; //global id for atoms 
  int l1,l2,l3; //local id for atoms 
  double x1[3],x2[3],x3[3];
  // each processor has nanglelist, see neighbor.cpp
  a0=0.0;
  //computeAreaVol(At,Vt,crossFlag);
  //for (int i = 0; i < nlocal; i++) domain->unmap(x[i],image[i],x1);
  computeAreaVol(At,Vt);
  check_crossing(crossFlag);
  /*if (comm->me == 0) {
    for (int j=1;j<nmolecule+1;j++)
        printf("%d molecule cross: %d %d %d\n", j, crossFlag[j][0],crossFlag[j][1],crossFlag[j][2]);
  }*/
  bigint ntimestep;
  ntimestep = update->ntimestep;
  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];
    molId = molecule[i2];
    /*g1=tag[i1];
    g2=tag[i2];
    g3=tag[i3];
    l1=atom->map(g1);
    l2=atom->map(g2);
    l3=atom->map(g3);
    l1 = domain->closest_image(i1,l1);
    l2 = domain->closest_image(i1,l2);
    l3 = domain->closest_image(i1,l3);*/
    //printf("proc: %d, %d/%d angles: i1 i2 i3:  %d %d %d local: %d global id %d %d %d\n",comm->me,n,nanglelist,i1,i2,i3,nlocal,g1,g2,g3); 
    /*positionShift(x[i1],x1,crossFlag[molId]); 
    positionShift(x[i2],x2,crossFlag[molId]); 
    positionShift(x[i3],x3,crossFlag[molId]);*/
    //if ((i1<nlocal) && (i2<nlocal) && (i3<nlocal)){
    /*domain->unmap(x[i1],image[i1],x1);
    domain->unmap(x[i2],image[i2],x2);
    domain->unmap(x[i3],image[i3],x3);*/

    //printf("proc: %d, %d/%d angles: i1 i2 i3:  %d %d %d local: %d global id %d %d %d unwrap x1 %lg %f %f x2 %g %g %g x3 %g %g %g wrap x1 %lg %f %f x2 %g %g %g x3 %g %g %g\n",comm->me,n,nanglelist,i1,i2,i3,nlocal,g1,g2,g3,x1[0],x1[1],x1[2],x2[0],x2[1],x2[2],x3[0],x3[1],x3[2],x[i1][0],x[i1][1],x[i1][2],x[i2][0],x[i2][1],x[i2][2],x[i3][0],x[i3][1],x[i3][2]); 
    
    // 1st bond 
    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];
    /*delx1 = x1[0] - x2[0];
    dely1 = x1[1] - x2[1];
    delz1 = x1[2] - x2[2];*/
    // 2nd bond
    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    /*delx2 = x3[0] - x2[0];
    dely2 = x3[1] - x2[1];
    delz2 = x3[2] - x2[2];*/
    // 3rd bond
    delx3 = x[i1][0] - x[i3][0];
    dely3 = x[i1][1] - x[i3][1];
    delz3 = x[i1][2] - x[i3][2];
    /*delx3 = x1[0] - x3[0];
    dely3 = x1[1] - x3[1];
    delz3 = x1[2] - x3[2];*/
   
    // norm xi=dx_2 x dx_1
    xi[0]=delz1*dely2 - dely1*delz2;
    xi[1]=delx1*delz2 - delz1*delx2;
    xi[2]=dely1*delx2 - delx1*dely2;
    
    xi2=xi[0]*xi[0] + xi[1]*xi[1] + xi[2]*xi[2];
    a0 = 0.5*sqrt(xi2);
    
    cnt[0]=0.333333333*(x[i1][0] + x[i2][0] + x[i3][0]);
    cnt[1]=0.333333333*(x[i1][1] + x[i2][1] + x[i3][1]);
    cnt[2]=0.333333333*(x[i1][2] + x[i2][2] + x[i3][2]);

    if (xperiodic){ 
        if ((crossFlag[molId][0]) && (cnt[0]<domain->xprd_half))
            cnt[0] += domain->xprd;    
    }
    if (yperiodic){ 
        if ((crossFlag[molId][1]) && (cnt[1]<domain->yprd_half))
            cnt[1] += domain->yprd;    
    }
    if (zperiodic){ 
        if ((crossFlag[molId][2]) && (cnt[2]<domain->zprd_half))
            cnt[2] += domain->zprd;    
    }
    /*if (comm->me==0){
        //printf("image %ld %ld %ld,i1-3: %d, %d, %d,  xunwrap %lg %lg %lg, x %lg %lg %lg\n",image[i1],image[i2],image[i3],i1,i2,i3,x1[0],x1[1],x1[2],x[i1][0],x[i1][1],x[i1][2]); 
        //printf("%d processor %d angle: dx1_unwrap: %lg %lg %lg, dx2_org: %lg %lg %lg\n",comm->me, n,delx1,dely1,delz1,x[i1][0] - x[i2][0],x[i1][1] - x[i2][1],x[i1][2] - x[i2][2]);
        printf("%d cell, crossFlag: %d %d %d cnt_unwrap: %lg %lg %lg, cnt_org: %lg %lg %lg\n",molId, crossFlag[molId][0],crossFlag[molId][1],crossFlag[molId][2],cnt[0],cnt[1],cnt[2],(x[i1][0] + x[i2][0] + x[i3][0])/3.0,(x[i1][1] + x[i2][1] + x[i3][1])/3.0,(x[i1][2] + x[i2][2] + x[i3][2])/3.0);
    }*/
    //global area
    beta_a = -0.25 * ka[type] * (At[molId] - A0t[type]) / (A0t[type] * a0);
    beta_ad = -0.25 * kd[type] * (a0 - A0[type])/(A0[type]*a0);
    beta_area = beta_a + beta_ad;
    // force & energy
    // xi x dx2
    f1[0] = beta_area*(xi[1]*delz2 - xi[2]*dely2);
    f1[1] = beta_area*(xi[2]*delx2 - xi[0]*delz2);
    f1[2] = beta_area*(xi[0]*dely2 - xi[1]*delx2);
    // xi x dx3
    f2[0] = beta_area*(xi[1]*delz3 - xi[2]*dely3);
    f2[1] = beta_area*(xi[2]*delx3 - xi[0]*delz3);
    f2[2] = beta_area*(xi[0]*dely3 - xi[1]*delx3);
    // xi x -dx1 = - (xi x dx1)
    f3[0] = -beta_area*(xi[1]*delz1 - xi[2]*dely1);
    f3[1] = -beta_area*(xi[2]*delx1 - xi[0]*delz1);
    f3[2] = -beta_area*(xi[0]*dely1 - xi[1]*delx1);
    
    /*fa1[0] = beta_area*(xi[1]*delz2 - xi[2]*dely2);
    fa1[1] = beta_area*(xi[2]*delx2 - xi[0]*delz2);
    fa1[2] = beta_area*(xi[0]*dely2 - xi[1]*delx2);
    fa2[0] = beta_area*(xi[1]*delz3 - xi[2]*dely3);
    fa2[1] = beta_area*(xi[2]*delx3 - xi[0]*delz3);
    fa2[2] = beta_area*(xi[0]*dely3 - xi[1]*delx3);
    fa3[0] = -beta_area*(xi[1]*delz1 - xi[2]*dely1);
    fa3[1] = -beta_area*(xi[2]*delx1 - xi[0]*delz1);
    fa3[2] = -beta_area*(xi[0]*dely1 - xi[1]*delx1);*/
   
    beta_v = -0.166666667 * kv[type] * (Vt[molId] - V0t[type]) / V0t[type];
    
    f1[0] += beta_v*(0.333333333*xi[0] + cnt[1]*delz2 - cnt[2]*dely2);
    f1[1] += beta_v*(0.333333333*xi[1] + cnt[2]*delx2 - cnt[0]*delz2);
    f1[2] += beta_v*(0.333333333*xi[2] + cnt[0]*dely2 - cnt[1]*delx2);

    f2[0] += beta_v*(0.333333333*xi[0] + cnt[1]*delz3 - cnt[2]*dely3);
    f2[1] += beta_v*(0.333333333*xi[1] + cnt[2]*delx3 - cnt[0]*delz3);
    f2[2] += beta_v*(0.333333333*xi[2] + cnt[0]*dely3 - cnt[1]*delx3);
    
    f3[0] += beta_v*(0.333333333*xi[0] - (cnt[1]*delz1 - cnt[2]*dely1));
    f3[1] += beta_v*(0.333333333*xi[1] - (cnt[2]*delx1 - cnt[0]*delz1));
    f3[2] += beta_v*(0.333333333*xi[2] - (cnt[0]*dely1 - cnt[1]*delx1));

    /*fv1[0] = beta_v*(xi[0]/3.0 + cnt[1]*delz2 - cnt[2]*dely2);
    fv1[1] = beta_v*(xi[1]/3.0 + cnt[2]*delx2 - cnt[0]*delz2);
    fv1[2] = beta_v*(xi[2]/3.0 + cnt[0]*dely2 - cnt[1]*delx2);
    fv2[0] = beta_v*(xi[0]/3.0 + cnt[1]*delz3 - cnt[2]*dely3);
    fv2[1] = beta_v*(xi[1]/3.0 + cnt[2]*delx3 - cnt[0]*delz3);
    fv2[2] = beta_v*(xi[2]/3.0 + cnt[0]*dely3 - cnt[1]*delx3);
    fv3[0] = beta_v*(xi[0]/3.0 - (cnt[1]*delz1 - cnt[2]*dely1));
    fv3[1] = beta_v*(xi[1]/3.0 - (cnt[2]*delx1 - cnt[0]*delz1));
    fv3[2] = beta_v*(xi[2]/3.0 - (cnt[0]*dely1 - cnt[1]*delx1));*/

    //debug
    /*if (comm->me == 0){
      //if (ntimestep > 38275){
      printf("%d cell: V0t %lg Vt %lg A0t %lg At %lg \n",molId,V0t[type],Vt[molId],A0t[type],At[molId]);// V0t[type], not i   
      printf("fa1 %lg %lg %lg, fv1 %lg %lg %lg \n",fa1[0], fa1[1],fa1[2],fv1[0],fv1[1],fv1[2]);// V0t[type], not i   
      printf("fa2 %lg %lg %lg, fv2 %lg %lg %lg \n",fa2[0], fa2[1],fa2[2],fv2[0],fv2[1],fv2[2]);// V0t[type], not i   
      printf("fa3 %lg %lg %lg, fv3 %lg %lg %lg \n",fa3[0], fa3[1],fa3[2],fv3[0],fv3[1],fv3[2]);// V0t[type], not i   
      //error->all(FLERR,"Incorrect args for angle coefficients");
      //}
    }*/
    // only local area conservation energy
    dA = a0 - A0[type];
    if (eflag) eangle = 0.5*kd[type]*dA*dA/A0[type];

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }
    // this may not be correct for angle_rbc style
    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
                         
     // } //if(i1<nlocal)
  } // nanglelist
  //for (int i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);
}

/* ---------------------------------------------------------------------- */

void AngleRbc::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(Cq,n+1,"angle:cq");
  memory->create(q,n+1,"angle:q");
  memory->create(ka,n+1,"angle:ka");
  memory->create(A0t,n+1,"angle:a0t");
  memory->create(kv,n+1,"angle:kv");
  memory->create(V0t,n+1,"angle:v0t");
  memory->create(kd,n+1,"angle:kd");
  memory->create(A0,n+1,"angle:a0");

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleRbc::coeff(int narg, char **arg)
{
  if (narg != 9) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  //force->bounds(FLERR,arg[0],atom->nangletypes,ilo,ihi);
  utils::bounds(FLERR,arg[0],1,atom->nangletypes,ilo,ihi,error);

  /*double Cq_one = force->numeric(FLERR,arg[1]);
  double q_one = force->numeric(FLERR,arg[2]);
  double ka_one = force->numeric(FLERR,arg[3]);
  double A0t_one = force->numeric(FLERR,arg[4]);
  double kv_one = force->numeric(FLERR,arg[5]);
  double V0t_one = force->numeric(FLERR,arg[6]);
  double kd_one = force->numeric(FLERR,arg[7]);
  double A0_one = force->numeric(FLERR,arg[8]);*/

  double Cq_one = utils::numeric(FLERR,arg[1],false,lmp);
  double q_one = utils::numeric(FLERR,arg[2],false,lmp);
  double ka_one = utils::numeric(FLERR,arg[3],false,lmp);
  double A0t_one = utils::numeric(FLERR,arg[4],false,lmp);
  double kv_one = utils::numeric(FLERR,arg[5],false,lmp);
  double V0t_one = utils::numeric(FLERR,arg[6],false,lmp);
  double kd_one = utils::numeric(FLERR,arg[7],false,lmp);
  double A0_one = utils::numeric(FLERR,arg[8],false,lmp);


  // convert theta0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    Cq[i] = Cq_one;
    q[i] = q_one;
    ka[i] = ka_one;
    A0t[i] = A0t_one;
    kv[i] = kv_one;
    V0t[i] = V0t_one;
    kd[i] = kd_one;
    A0[i] = A0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleRbc::equilibrium_angle(int i)
{
  return 0.0;// Not needed for Area/volume conservation
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleRbc::write_restart(FILE *fp)
{
  fwrite(&Cq[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&q[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ka[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&A0t[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&kv[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&V0t[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&kd[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&A0[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleRbc::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&Cq[1],sizeof(double),atom->nangletypes,fp);
    fread(&q[1],sizeof(double),atom->nangletypes,fp);
    fread(&ka[1],sizeof(double),atom->nangletypes,fp);
    fread(&A0t[1],sizeof(double),atom->nangletypes,fp);
    fread(&kv[1],sizeof(double),atom->nangletypes,fp);
    fread(&V0t[1],sizeof(double),atom->nangletypes,fp);
    fread(&kv[1],sizeof(double),atom->nangletypes,fp);
    fread(&A0[1],sizeof(double),atom->nangletypes,fp);
  }
  MPI_Bcast(&Cq[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&q[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ka[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&A0t[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kv[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&V0t[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kd[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&A0[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleRbc::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g %g %g\n",i,Cq[i],q[i],
            ka[i],A0t[i],kv[i],V0t[i],kd[i],A0[i]);
}

/* ---------------------------------------------------------------------- */

double AngleRbc::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;
  double a0, xi[3], xi2;
  
  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);
  //area
  xi[0]=delz1*dely2 - dely1*delz2;
  xi[1]=delx1*delz2 - delz1*delx2;
  xi[2]=dely1*delx2 - delx1*dely2;
    
  xi2=xi[0]*xi[0] + xi[1]*xi[1] + xi[2]*xi[2];
  a0 = 0.5*sqrt(xi2);

  double dtheta = a0 - A0[type];
  double tk = 0.5 * kd[type] * dtheta;
  return tk*dtheta/A0[type]; 
}

void AngleRbc::minimum_image_xyz(double &dx, double &dy, double &dz, int &xflag, int &yflag, int &zflag)
{
  int triclinic = domain->triclinic;

  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;
  
  double xprd_half = domain->xprd_half;
  double yprd_half = domain->yprd_half;
  double zprd_half = domain->zprd_half;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double xy = domain->xy;
  double yz = domain->yz;
  double xz = domain->xz;
 
  if (triclinic == 0) {
    if (xperiodic) {
      if (fabs(dx) > xprd_half) {
        xflag = 1;
        if (dx < 0.0) dx += xprd;
        else dx -= xprd;
      }
    }
    if (yperiodic) {
      if (fabs(dy) > yprd_half) {
        yflag = 1;
        if (dy < 0.0) dy += yprd;
        else dy -= yprd;
      }
    }
    if (zperiodic) {
      if (fabs(dz) > zprd_half) {
        zflag = 1;
        if (dz < 0.0) dz += zprd;
        else dz -= zprd;
      }
    }

  } else {
    if (zperiodic) {
      if (fabs(dz) > zprd_half) {
        zflag = 1;
        if (dz < 0.0) {
          dz += zprd;
          dy += yz;
          dx += xz;
        } else {
          dz -= zprd;
          dy -= yz;
          dx -= xz;
        }
      }
    }
    if (yperiodic) {
      if (fabs(dy) > yprd_half) {
        yflag = 1;
        if (dy < 0.0) {
          dy += yprd;
          dx += xy;
        } else {
          dy -= yprd;
          dx -= xy;
        }
      }
    }
    if (xperiodic) {
      if (fabs(dx) > xprd_half) {
        xflag = 1;
        if (dx < 0.0) dx += xprd;
        else dx -= xprd;
      }
    }
  }
}

void AngleRbc::check_crossing(int **crossFlag)
{
  //check if cell is crossing the boundary
  //return cross[ncell][3]
  //double maxX[3],minX[3];
  //double gmax[3],gmin[3];
  int molID;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  tagint *molecule = atom->molecule;
  int imagePBC;
  double maxImag, minImag;
  tagint *image = atom->image;
  for (int j = 1; j < nmolecule+1; j++) {
    /*maxxyz[j][0]=domain->sublo[0];
    maxxyz[j][1]=domain->sublo[1];
    maxxyz[j][2]=domain->sublo[2];
    minxyz[j][0]=domain->subhi[0];
    minxyz[j][1]=domain->subhi[1];
    minxyz[j][2]=domain->subhi[2];*/
    maxxyz[j][0]=domain->boxlo[0];
    maxxyz[j][1]=domain->boxlo[1];
    maxxyz[j][2]=domain->boxlo[2];
    minxyz[j][0]=domain->boxhi[0];
    minxyz[j][1]=domain->boxhi[1];
    minxyz[j][2]=domain->boxhi[2];
    crossFlag[j][0]=0;
    crossFlag[j][1]=0;
    crossFlag[j][2]=0;
  }
  imagePBC=0; 
  for (int i=0;i<nlocal;i++){
        molID = molecule[i];
        //printf("proc: %d, ilocal %d image "TAGINT_FORMAT"\n",comm->me, i, image[i]);
        //if (image[i] != image[0]) imagePBC=1; //works for 1 periodic, but for multiple PBCs, may have problems. Jifu 3/10/2020
        if (maxxyz[molID][0]<x[i][0]) maxxyz[molID][0] = x[i][0];   
        if (maxxyz[molID][1]<x[i][1]) maxxyz[molID][1] = x[i][1];   
        if (maxxyz[molID][2]<x[i][2]) maxxyz[molID][2] = x[i][2];   
        if (minxyz[molID][0]>x[i][0]) minxyz[molID][0] = x[i][0];   
        if (minxyz[molID][1]>x[i][1]) minxyz[molID][1] = x[i][1];   
        if (minxyz[molID][2]>x[i][2]) minxyz[molID][2] = x[i][2];   
    }
    MPI_Allreduce(&maxxyz[0][0],&MAXxyz[0][0],3*(nmolecule+1),MPI_DOUBLE,MPI_MAX,world);
    MPI_Allreduce(&minxyz[0][0],&MINxyz[0][0],3*(nmolecule+1),MPI_DOUBLE,MPI_MIN,world);
    /*if (comm->me == 0) {
    for (int j=1;j<nmolecule+1;j++)
        printf("%d molecule cross: %f %f %f %f %f %f\n", j, MAXxyz[j][0],MAXxyz[j][1],MAXxyz[j][2], MINxyz[j][0],MINxyz[j][1],MINxyz[j][2]);
  }*/
    //printf("image pbc %d\n",imagePBC);
  for (int j = 1; j < nmolecule+1; j++) {
    if (((MAXxyz[j][0]-MINxyz[j][0])>domain->xprd_half) && (domain->xperiodic)) crossFlag[j][0]=1;
    if (((MAXxyz[j][1]-MINxyz[j][1])>domain->yprd_half) && (domain->yperiodic)) crossFlag[j][1]=1;
    if (((MAXxyz[j][2]-MINxyz[j][2])>domain->zprd_half) && (domain->zperiodic)) crossFlag[j][2]=1;
  }
    //std::vector< std::vector<double> > xx;
  //new code
  /*std::set<int> cell;
  std::set<int>::iterator it; 
  for (int i=0;i<nlocal;i++){
    cell.insert(molecule[i]);
  }
  for (it=cell.begin();it!=cell.end();it++){
    xx.clear();
    for (int i=0;i<nlocal;i++){
      molId = molecule[i];
      if (*it == molId){
        std::vector<double> tmp;
        tmp.push_back(x[i][0]);
        tmp.push_back(x[i][1]);
        tmp.push_back(x[i][2]);
        xx.push_back(tmp);
      }
    }
    maxX[0]=domain->sublo[0];
    maxX[1]=domain->sublo[1];
    maxX[2]=domain->sublo[2];
    minX[0]=domain->subhi[0];
    minX[1]=domain->subhi[1];
    minX[2]=domain->subhi[2];
    for (int i=0;i<xx.size();i++){
      if (xx[i][0] > maxX[0]) maxX[0] = xx[i][0];
      if (xx[i][1] > maxX[1]) maxX[1] = xx[i][1];
      if (xx[i][2] > maxX[2]) maxX[2] = xx[i][2];
      if (xx[i][0] < minX[0]) minX[0] = xx[i][0];
      if (xx[i][1] < minX[1]) minX[1] = xx[i][1];
      if (xx[i][2] < minX[2]) minX[2] = xx[i][2];
    }
    MPI_Allreduce(&maxX[0],&gmax[0],3,MPI_DOUBLE,MPI_MAX,world);
    MPI_Allreduce(&minX[0],&gmin[0],3,MPI_DOUBLE,MPI_MIN,world);
    if ((gmax[0]-gmin[0])>domain->xprd_half) crossFlag[*it][0]=1;
    if ((gmax[1]-gmin[1])>domain->yprd_half) crossFlag[*it][1]=1;
    if ((gmax[2]-gmin[2])>domain->zprd_half) crossFlag[*it][2]=1;
  }*/
 /* 
  for (int j = 1; j < nmolecule+1; j++) {
    xx.clear();
    for (int i=0;i<nlocal;i++){
      molId = molecule[i];
      if (j == molId){
        std::vector<double> tmp;
        tmp.push_back(x[i][0]);
        tmp.push_back(x[i][1]);
        tmp.push_back(x[i][2]);
        xx.push_back(tmp);
      }
    }
    maxX[0]=domain->sublo[0];
    maxX[1]=domain->sublo[1];
    maxX[2]=domain->sublo[2];
    minX[0]=domain->subhi[0];
    minX[1]=domain->subhi[1];
    minX[2]=domain->subhi[2];
    for (int i=0;i<xx.size();i++){
      if (xx[i][0] > maxX[0]) maxX[0] = xx[i][0];
      if (xx[i][1] > maxX[1]) maxX[1] = xx[i][1];
      if (xx[i][2] > maxX[2]) maxX[2] = xx[i][2];
      if (xx[i][0] < minX[0]) minX[0] = xx[i][0];
      if (xx[i][1] < minX[1]) minX[1] = xx[i][1];
      if (xx[i][2] < minX[2]) minX[2] = xx[i][2];
    }
    MPI_Allreduce(&maxX[0],&gmax[0],3,MPI_DOUBLE,MPI_MAX,world);
    MPI_Allreduce(&minX[0],&gmin[0],3,MPI_DOUBLE,MPI_MIN,world);
    if ((gmax[0]-gmin[0])>domain->xprd_half) crossFlag[j][0]=1;
    if ((gmax[1]-gmin[1])>domain->yprd_half) crossFlag[j][1]=1;
    if ((gmax[2]-gmin[2])>domain->zprd_half) crossFlag[j][2]=1;
  }*/
}

void AngleRbc::positionShift(double *x, double *xnew, int *crossFlag){
  xnew[0]=x[0];
  xnew[1]=x[1];
  xnew[2]=x[2];
  if (crossFlag[0] == 1){
    double xprd_half = domain->xprd_half;  
    double xprd = domain->xprd;
    if (x[0] < domain->boxlo[0] + xprd_half) xnew[0] = x[0] + xprd;
  }
  if (crossFlag[1] == 1){
    double yprd_half = domain->yprd_half;  
    double yprd = domain->yprd;
    if (x[1] < domain->boxlo[1] + yprd_half) xnew[1] = x[1] + yprd;
  }
  if (crossFlag[2] == 1){
    double zprd_half = domain->zprd_half; 
    double zprd = domain->zprd;
    if (x[2] < domain->boxlo[2] + zprd_half) xnew[2] = x[2] + zprd;
  }    
}
