#pragma once
#include "palabos3D.h"
#include "palabos3D.hh"
#include <DataAdaptor.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkImageData.h>
using namespace plb; 
namespace senseiLP 
{
class LPDataAdaptor : public sensei::DataAdaptor
{
public:

  static LPDataAdaptor* New();

  senseiTypeMacro(LPDataAdaptor, sensei::DataAdaptor);

  void Initialize();

  void AddLAMMPSData(double **x, long ntimestep, int nghost, 
                     int nlocal, int **anglelist, int nanglelist);

   void AddPalabosData(vtkDoubleArray *velocityDoubleArray,
                      vtkDoubleArray *vorticityDoubleArray,
                      vtkDoubleArray *velocityNormDoubleArray,
                      int nx, int ny, int nz, Box3D domainBox, plint envelopeWidth);
// SENSEI API (Virtual functions overridden from sensei/DataAdaptor.h)
  int GetNumberOfMeshes(unsigned int &numMeshes) override;

  int GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &metadata) override;

  int GetMesh(const std::string &meshName, bool structureOnly, vtkDataObject *&mesh) override ;

  int GetMesh(const std::string &meshName, bool structureOnly, vtkCompositeDataSet *&mesh) override;

  int AddGhostNodesArray(vtkDataObject* mesh, const std::string &meshName) override;

  int AddGhostCellsArray(vtkDataObject* mesh, const std::string &meshName) override;

  int AddArray(vtkDataObject* mesh, const std::string &meshName, int association, const std::string &arrayName) override;

  int AddArrays(vtkDataObject* mesh, const std::string &meshName, int association, const std::vector<std::string> &arrayName) override;

  int ReleaseData() override;

protected:

LPDataAdaptor();

~LPDataAdaptor();

private:

  struct DInternals;
  DInternals* Internals;
};

}
