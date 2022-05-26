/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS

AngleStyle(rbc,AngleRbc)

#else

#ifndef LMP_ANGLE_RBC_H
#define LMP_ANGLE_RBC_H

#include "stdio.h"
#include "angle.h"

namespace LAMMPS_NS {

class AngleRbc : public Angle {
 public:
  AngleRbc(class LAMMPS *);
  virtual ~AngleRbc();
  virtual void init_style();//calculate equilibrium angle
  virtual void compute(int, int);
  virtual void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, int, int, int);
  void computeAreaVol(double *, double *, int **);
  void computeAreaVol(double *, double *);//optimized projection method
  void minimum_image_xyz(double &dx, double &dy, double &dz, int &xflag, int &yflag, int &zflag);
  void check_crossing(int **crossFlag);
  void positionShift(double *, double *, int *);

 protected:
  double *Cq,*q,*ka,*A0t,*kv,*V0t,*kd,*A0;
  int nmolecule; 
  double *At,*Vt,*at,*vt;//global and local vairables
  int **crossFlag;
  double **MAXxyz, **maxxyz,**MINxyz, **minxyz;
  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for angle coefficients

Self-explanatory.  Check the input script or data file.

*/
