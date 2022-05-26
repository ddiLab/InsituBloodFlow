#pragma once
#include <mpi.h>
#include <string>
#include "palabos3D.h"
#include "palabos3D.hh"
using namespace plb; 
namespace Bridge
{
  void Initialize(MPI_Comm world, const std::string& config_file);
  void SetData(double **x, long ntimestep, int nghost, 
               int nlocal, int **anglelist, int nanglelist,
	           TensorField3D<double, 3> velocityDoubleArray, 
	           TensorField3D<double, 3> vorticityDoubleArray, 
	           ScalarField3D<double> velocityNormDoubleArray,
             int nx, int ny, int nz, Box3D domainBox, plint envelopeWidth); 
  void Analyze(long ntimestep); 
  void Finalize();
}

