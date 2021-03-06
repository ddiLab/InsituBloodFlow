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

#ifdef FIX_CLASS

FixStyle(activate/platelet,FixActivatePlatelet)

#else

#ifndef LMP_FIX_ActivatePlatelet_H
#define LMP_FIX_ActivatePlatelet_H

#include "fix.h"

namespace LAMMPS_NS {

class FixActivatePlatelet : public Fix {
 public:
  double **fexternal;
  double rate;

  FixActivatePlatelet(class LAMMPS *, int, char **);
  ~FixActivatePlatelet();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void min_post_force(int);
  void end_of_step();
  //double compute_scalar();
  //void set_energy(double eng);

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  double releaseRate(){return rate;};
  //typedef void (*FnPtr)(void *, bigint, int, tagint *, double **, double **);
  //void set_callback(FnPtr, void *);

 private:
  //int mode,ncall,napply;
  double xi0, m0, tau, alpha, n;
  /*  FnPtr callback;
  void *ptr_caller;
  double user_energy;*/
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix ActivatePlatelet callback function not set

This must be done by an ActivatePlatelet program in order to use this fix.

*/
