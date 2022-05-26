#ifndef IBM_LBM_2D_H
#define IBM_LBM_2D_H
#include "lammpsWrapper.h"

namespace plb {
    
  template<typename T, template<typename U> class Descriptor>
  void interpolateVelocity2D(MultiBlockLattice2D<T,Descriptor> &lattice,
                           MultiTensorField2D<T,2> &velocity,
                           LammpsWrapper &wrapper);

  template<typename T, template<typename U> class Descriptor>
  void spreadForce2D(MultiBlockLattice2D<T,Descriptor> &lattice,  
                   LammpsWrapper &wrapper );
  
}; /* namespace plb */

#include "ibm2D.hh"

#endif /* IBDATAEXCHANGEWRAPPERS_H_LBDEM */
