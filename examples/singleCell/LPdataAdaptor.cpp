#include "LPdataAdaptor.h" 
#include "palabos3D.h"
#include "palabos3D.hh"

#include "Error.h"
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>

#include "Bridge.h"
#include <iostream>
using namespace std;
using namespace plb; 

static
vtkUnsignedCharArray *newGhostCellsArray(plb::Box3D domain, int ng, int gnx, int gny, int gnz)
{
    // This sim is always 3D.
    int imin = domain.x0;
    int jmin = domain.y0;
    int kmin = domain.z0;
    int imax = domain.x1;
    int jmax = domain.y1;
    int kmax = domain.z1;
    int nx = domain.getNx()-1;
    int ny = domain.getNy()-1;
    int nz = domain.getNz()-1;
    int nxny = nx*ny;
    int ncells = nx*ny*nz;
    cout << "X Extents (domain.getNX()): " << domain.getNx() << endl;
    vtkUnsignedCharArray *g = vtkUnsignedCharArray::New();
    g->SetNumberOfTuples(ncells);
    memset(g->GetVoidPointer(0), 0, sizeof(unsigned char) * ncells);
    g->SetName("vtkTestType");
    
    unsigned char *gptr = (unsigned char *)g->GetVoidPointer(0);
    unsigned char ghost = 1;
    unsigned char external = 1;
   
    //I Low********************************************
    if(imin < 0) 
    {
      // Set the lowest I faces to external surface.
      for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ny; ++j)
      for(int i = 0; i < ng; ++i)
          gptr[k * nxny + j*nx + i] = external;
    }
    else
    { 
      // Set the low I faces to ghosts.
      for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ny; ++j)
      for(int i = 0; i < ng; ++i)
          gptr[k * nxny + j*nx + i] = ghost;
    }   
    //*************************************************
    
    //I High*******************************************
    if(imax > gnx)
    {
      // Set the highest I faces to external.
      for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ny; ++j)
      for(int i = nx-ng; i < nx; ++i)
          gptr[k * nxny + j*nx + i] = external;
    }
    else{
      // Set the high I faces to ghosts.
      for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ny; ++j)
      for(int i = nx-ng+1; i < nx; ++i)
          gptr[k * nxny + j*nx + i] = ghost;
    }
    //**************************************************
    
    //J Low*********************************************
    if(jmin < 0 )
    {
      // Set the lowest J faces to external.
      for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ng; ++j)
      for(int i = 0; i < nx; ++i)
          gptr[k * nxny + j*nx + i] = external;
    }
    else
    {
      // Set the low J faces to ghosts.
      for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ng; ++j)
      for(int i = 0; i < nx; ++i)
          gptr[k * nxny + j*nx + i] = ghost;
    }
    //**************************************************

    //J High********************************************
    if(jmax > gny)
    {
      // Set the highest J faces to external.
      for(int k = 0; k < nz; ++k)
      for(int j = ny-ng; j < ny; ++j)
      for(int i = 0; i < nx; ++i)
          gptr[k * nxny + j*nx + i] = external;
    }
    else
    {
        // Set the high J faces to ghosts.
      for(int k = 0; k < nz; ++k)
      for(int j = ny-ng+1; j < ny; ++j)
      for(int i = 0; i < nx; ++i)
          gptr[k * nxny + j*nx + i] = ghost;
    }
    //***************************************************

    //K Low**********************************************
    // Set the low K faces to ghosts.
    if(kmin < 0 )
    {
      // Set the lowest K faces to external.
      for(int k = 0; k < ng; ++k)
      for(int j = 0; j < ny; ++j)
      for(int i = 0; i < nx; ++i)
          gptr[k * nxny + j*nx + i] = external;
    }
    else
    {
      // Set the low K faces to ghosts.
      for(int k = 0; k < ng; ++k)
      for(int j = 0; j < ny; ++j)
      for(int i = 0; i < nx; ++i)
          gptr[k * nxny + j*nx + i] = ghost;
    }
    //****************************************************

    //K High**********************************************
    // Set the high K faces to ghosts.
    if(kmax > gnz)
    {
      // Set the high K faces to ghosts.
      for(int k = nz-ng; k < nz; ++k)
      for(int j = 0; j < ny; ++j)
      for(int i = 0; i < nx; ++i)
          gptr[k * nxny + j*nx + i] = external;
    }
    else
    {
      // Set the high K faces to ghosts.
      for(int k = nz-ng+1; k < nz; ++k)
      for(int j = 0; j < ny; ++j)
      for(int i = 0; i < nx; ++i)
          gptr[k * nxny + j*nx + i] = ghost;
    }    
    //****************************************************
    
    //gptr[0 * nxny + 0*nx + 0] = surface;
    //gptr[1 * nxny + 1*nx + 0] = surface;
    //gptr[2 * nxny + 2*nx + 0] = surface;
    
    return g;
}


namespace senseiLP
{
  //----------------------------------------------------------------------
  struct LPDataAdaptor::DInternals
  {
    vtkSmartPointer<vtkMultiBlockDataSet> mesh;
    vtkSmartPointer<vtkDoubleArray> AtomPositions;
    vtkSmartPointer<vtkIntArray> AtomTypes;
    vtkSmartPointer<vtkIntArray> AtomIDs;
    vtkSmartPointer<vtkCellArray> vertices;
  
  // -------- PALABOS ---------------------
    vtkDoubleArray *pb_velocityDoubleArray;
    vtkDoubleArray *pb_vorticityDoubleArray; 
    vtkDoubleArray *pb_velocityNormDoubleArray;
  // --------------------------------------
    long NumBlocks;
    int nlocal, nghost, nanglelist;
    double **x;
    int *type;
    int **anglelist;
    int *id;
    int pb_nx, pb_ny, pb_nz;
    Box3D domainBox; //XXX added for domainBox 2/23/22
    int envelopeWidth;
    
  };
  //----------------------------------------------------------------------
  senseiNewMacro(LPDataAdaptor);
  //----------------------------------------------------------------------
  LPDataAdaptor::LPDataAdaptor() :
    Internals(new LPDataAdaptor::DInternals())
  {
  }
  //----------------------------------------------------------------------
  LPDataAdaptor::~LPDataAdaptor()
  {
    delete this->Internals;
  }
  //----------------------------------------------------------------------
  void LPDataAdaptor::Initialize()
  {
    this->ReleaseData();//ReleaseData must be correctly defined!!
  }
  //----------------------------------------------------------------------
  void LPDataAdaptor::AddLAMMPSData(double **x, long ntimestep, int nghost, 
                                    int nlocal, int **anglelist, int nanglelist)
  {
    

    DInternals& internals = (*this->Internals);
    
    if(!internals.AtomPositions)
    {
      internals.AtomPositions = vtkSmartPointer<vtkDoubleArray>::New();
    }
    
  // atom coordinates
    if (internals.AtomPositions)
    {
      long nvals = nlocal+nghost;

      
      internals.AtomPositions->SetNumberOfComponents(3);
      internals.AtomPositions->SetArray(*x, nvals*3, 1);
      internals.AtomPositions->SetName("positions");
      
      internals.x = x;
    }  
    else
    {
      SENSEI_ERROR("Error. Internal AtomPositions structure not initialized")
    }

    

  // anglelists
    internals.anglelist = anglelist;
    internals.nanglelist = nanglelist;

  // number of atoms
    internals.nlocal = nlocal;
    internals.nghost = nghost;

  
  // timestep
    this->SetDataTimeStep(ntimestep);//XXX 4/20/22 There is a SetDataTimeStep in the Analyze() function in Bridge.cpp. Is this redundant?

  }
  //---------------------------------------------------------------------------
  void LPDataAdaptor::AddPalabosData(vtkDoubleArray *velocityDoubleArray,
                vtkDoubleArray *vorticityDoubleArray,
                vtkDoubleArray *velocityNormDoubleArray,
                        int nx, int ny, int nz, Box3D domainBox, plint envelopeWidth) 
  {
    DInternals& internals = (*this->Internals);

    internals.pb_velocityDoubleArray = velocityDoubleArray;
    internals.pb_vorticityDoubleArray = vorticityDoubleArray; 
    internals.pb_velocityNormDoubleArray = velocityNormDoubleArray;
    internals.pb_velocityDoubleArray->SetName("velocity");
    internals.pb_vorticityDoubleArray->SetName("vorticity"); 
    internals.pb_velocityNormDoubleArray->SetName("velocityNorm");
      
    internals.pb_nx = nx;
    internals.pb_ny = ny;
    internals.pb_nz = nz;
    internals.domainBox = domainBox;//XXX domainBox 2/23/22
    internals.envelopeWidth = envelopeWidth;
  }   
  //----------------------------------------------------------------------
  int LPDataAdaptor::GetNumberOfMeshes(unsigned int &numMeshes)
  {
    numMeshes = 2;
    return 0;
  }
  //----------------------------------------------------------------------
  int LPDataAdaptor::GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &metadata) 
  {
    //cout << "ID NUMBER: " << id << endl;
    //cout << "Calling GetMeshMetaData" << endl;
    int rank, nRanks;
    
    int nx = this->Internals->pb_nx;
    int ny = this->Internals->pb_ny;
    int nz = this->Internals->pb_nz; 

    //XXX Added for domainBox 2/23/22********	
    Box3D domainBox = this->Internals->domainBox;
    int nlx = domainBox.getNx(); 
    int nly = domainBox.getNy();
    int nlz = domainBox.getNz();
    plb::Array<plint, 6> localExtents = domainBox.to_plbArray();//XXX look at palabos/src/core/geometry3D.h for documentation
    //***************************************

    MPI_Comm_rank(this->GetCommunicator(), &rank);
    MPI_Comm_size(this->GetCommunicator(), &nRanks); 	

    if (id == 0) // id == 0 is cells
    {
      //cout << "GetMeshMetaData Cells Test" << endl;
      metadata->MeshName = "cells";
      metadata->MeshType = VTK_MULTIBLOCK_DATA_SET; //VTK_POLY_DATA;
      metadata->BlockType = VTK_POLY_DATA;
      metadata->CoordinateType = VTK_DOUBLE;
      metadata->NumBlocks = nRanks;
      metadata->NumBlocksLocal = {1};
      metadata->NumGhostCells = 0;
      metadata->NumArrays = 0;
      metadata->StaticMesh = 0;  

      if (metadata->Flags.BlockExtentsSet())
      {
        //SENSEI_WARNING("lammps data adaptor. Flags.BlockExtentsSet()")
        
        // There should be no extent for a PolyData, but ADIOS2 needs this
        std::array<int,6> ext = { 0, 0, 0, 0, 0, 0};
        metadata->Extent = std::move(ext);
        metadata->BlockExtents.reserve(1);	// One block per rank
        metadata->BlockExtents.emplace_back(std::move(ext));
      }
      
      if (metadata->Flags.BlockDecompSet())
      {
        metadata->BlockOwner.push_back(rank);
        metadata->BlockIds.push_back(rank);
      }
    
      //We use nanglelist for BlockNumCells because it give the number of triangles on a given processor
      metadata->BlockNumCells.push_back(this->Internals->nanglelist);
      metadata->BlockNumPoints.push_back(this->Internals->nlocal + this->Internals->nghost );
      metadata->BlockCellArraySize.push_back(0);
    }
    else if(id == 1) // id == 1 is fluid
    {
      //cout << "GetMeshMetaData Fluid Test" << endl;
      metadata->MeshName = "fluid"; 
      metadata->MeshType = VTK_MULTIBLOCK_DATA_SET;
      metadata->BlockType= VTK_IMAGE_DATA; 
      metadata->CoordinateType = VTK_DOUBLE;
      metadata->NumBlocks = nRanks;
      metadata->NumBlocksLocal = {1};
      metadata->NumGhostCells = this->Internals->envelopeWidth;  
      metadata->NumArrays=3;
      metadata->ArrayName = {"velocity","vorticity","velocityNorm"};
      metadata->ArrayComponents = {3, 3, 1}; 
      metadata->ArrayType = {VTK_DOUBLE, VTK_DOUBLE, VTK_DOUBLE};
      metadata->ArrayCentering = {vtkDataObject::POINT, vtkDataObject::POINT, vtkDataObject::POINT};
      metadata->StaticMesh = 1; 

      if (metadata->Flags.BlockDecompSet())
      {
        metadata->BlockOwner.push_back(rank);
        metadata->BlockIds.push_back(rank);
      }

      if (metadata->Flags.BlockExtentsSet())
      {
        //SENSEI_WARNING("lammps data adaptor. Flags.BlockExtentsSet()")
        std::array<int,6> ext = { 0, nx, 0, ny, 0, nz };
        std::array<int,6> blockext = {localExtents[0], localExtents[1], localExtents[2], localExtents[3], localExtents[4], localExtents[5]}; //XXX Changes 2/23/22
        cout << " CHECK DOMAIN: " << localExtents[0] << " " << localExtents[1] << " " << localExtents[2] << " " << localExtents[3] << " " << localExtents[4] << " " << localExtents[5] << endl;
        metadata->Extent = std::move(ext);
        metadata->BlockExtents.reserve(1);	// One block per rank
        metadata->BlockExtents.emplace_back(std::move(blockext)); //XXX We have to figure out the local numbers for block ext
      }

      metadata->BlockNumCells.push_back((nlx-1) * (nly-1) * (nlz-1)); //(nlx * nly * nlz * 3);//XXX Changed 2/23/22
      metadata->BlockNumPoints.push_back((nlx) * (nly) * (nlz)); //(nlx * nly * nlz * 3); //XXX Changed 2/23/22
      metadata->BlockCellArraySize.push_back(0); 
    }
    else
    {
      SENSEI_ERROR("MeshMetaData Error: id value does not exist")
    }
    return 0;
  }
  //----------------------------------------------------------------------
  int LPDataAdaptor::GetMesh(const std::string &meshName, bool structureOnly, vtkDataObject *&mesh)
  {
    int rank, size; 
    MPI_Comm_rank(this->GetCommunicator(), &rank);
    MPI_Comm_size(this->GetCommunicator(), &size);
    mesh = nullptr;
    //cout << "Calling GetMesh" << endl;
    if(meshName == "cells")
    {  
      DInternals& internals = (*this->Internals);

      vtkPolyData *pd = vtkPolyData::New();
      vtkMultiBlockDataSet *mb = vtkMultiBlockDataSet::New();

      if(!structureOnly)
      {
        vtkPoints *pts = vtkPoints::New();
        //pts->SetNumberOfPoints(internals.nlocal+internals.nghost);
        pts->SetData(internals.AtomPositions);
        vtkCellArray *Triangles = vtkCellArray::New();

        for (int i = 0 ; i < internals.nanglelist ; i++)
        {
          //cout << rank << " : Triangle Test :" << i << endl;
          vtkTriangle *Triangle = vtkTriangle::New();
          Triangle->GetPointIds()->SetId(0, internals.anglelist[i][0]);
          Triangle->GetPointIds()->SetId(1, internals.anglelist[i][1]);
          Triangle->GetPointIds()->SetId(2, internals.anglelist[i][2]);
          Triangles->InsertNextCell(Triangle);
	  
	  Triangle->Delete();
        }
        pd->SetPoints(pts);
        pd->SetPolys(Triangles);
	pts->Delete();
	Triangles->Delete();
      }
      
      pd->SetVerts( internals.vertices ); //XXX Does this do anything? vertices gets created in DInternals but doesn't get set as anything

      mb->SetNumberOfBlocks(size);
      mb->SetBlock(rank,pd);
 
      mesh = mb;
      pd->Delete();
    }

    else if(meshName == "fluid")
    {
      
      DInternals& internals = (*this->Internals); 
      //XXX added 2/24/22**********************
      Box3D domainBox = this->Internals->domainBox;
      int nlx = domainBox.getNx(); 
      int nly = domainBox.getNy();
      int nlz = domainBox.getNz();
      //***************************************

      vtkMultiBlockDataSet *mbfluid = vtkMultiBlockDataSet::New();
      
      //cout << "Inside get mesh " << meshName << endl;

      vtkImageData *FluidImageData = vtkImageData::New();
      FluidImageData->SetDimensions(nlx, nly, nlz); //XXX Changed on 2/24/22 (+2 because of extra layer on each side of vtk file)
      FluidImageData->SetExtent(domainBox.x0, domainBox.x1, domainBox.y0, domainBox.y1, domainBox.z0, domainBox.z1);

      //cout << internals.pb_nx << "," << internals.pb_ny << "," << internals.pb_nz << endl;
      /*
      FluidImageData->GetPointData()->AddArray(internals.pb_velocityDoubleArray);
      internals.pb_velocityDoubleArray->SetName("velocity");

      FluidImageData->GetPointData()->AddArray(internals.pb_vorticityDoubleArray);
      internals.pb_vorticityDoubleArray->SetName("vorticity"); 

      FluidImageData->GetPointData()->AddArray(internals.pb_velocityNormDoubleArray);
      internals.pb_velocityNormDoubleArray->SetName("velocityNorm");
      */

      mbfluid->SetNumberOfBlocks(size);
      mbfluid->SetBlock(rank,FluidImageData);
      FluidImageData->Delete();

      mesh = mbfluid;
    }
    else
    {
      SENSEI_ERROR("No mesh \"" << meshName << "\"")
      return -1;
    }

    return 0;
  }
  //----------------------------------------------------------------------
  int LPDataAdaptor::GetMesh(const std::string &meshName, bool structureOnly, vtkCompositeDataSet *&mesh)
  {
    return 0;
  }
  //----------------------------------------------------------------------
  int LPDataAdaptor::AddGhostNodesArray(vtkDataObject* mesh, const std::string &meshName)
  {
   
    return 0;
  }
  //----------------------------------------------------------------------
  int LPDataAdaptor::AddGhostCellsArray(vtkDataObject* mesh, const std::string &meshName)
  {
    //cout << " TESTING NODE ARRAY" << endl;
    
    int rank;
    MPI_Comm_rank(this->GetCommunicator(), &rank);
    if(meshName == "fluid")
    {
      DInternals& internals = (*this->Internals); 
      
      vtkMultiBlockDataSet *mbfluid = dynamic_cast<vtkMultiBlockDataSet*>(mesh); 
      if(!mbfluid)   
      {
        SENSEI_ERROR("unexpected mesh type "<< (mesh ? mesh->GetClassName() : "nullptr"))
        return -1;
      }
      vtkImageData *FluidImageData = (vtkImageData*)mbfluid->GetBlock(rank);
      if(!FluidImageData)
      {
        SENSEI_ERROR("Cannot Get Block in LPDataAdaptor::AddArray")
        return -1;
      }
      
      vtkDataSetAttributes *dsa = FluidImageData->GetAttributes(vtkDataObject::CELL);

      vtkUnsignedCharArray *ga = newGhostCellsArray(this->Internals->domainBox, this->Internals->envelopeWidth , this->Internals->pb_nx, this->Internals->pb_ny, this->Internals->pb_nz);//XXX Fix this to envelopeWidth
      dsa->AddArray(ga);
      ga->Delete(); 
    }
    
    return 0;
  }
  //----------------------------------------------------------------------
  int LPDataAdaptor::AddArray(vtkDataObject* mesh, const std::string &meshName,
      int association, const std::string &arrayName)
  {
    //cout << "meshname: " << meshName<< "  ArrayName: " << arrayName << endl;
    int rank;
    MPI_Comm_rank(this->GetCommunicator(), &rank);
    if(meshName == "fluid")
    {
      DInternals& internals = (*this->Internals); 
      /*
      vtkImageData *velocity = dynamic_cast<vtkImageData*>(mesh); //XXX PROBABLY NEEDS TO BE vtkmultiblockdataset
      vtkImageData *vorticity = dynamic_cast<vtkImageData*>(mesh); //XXX extract image data and then add arrays
      vtkImageData *velocityNorm = dynamic_cast<vtkImageData*>(mesh); //XXX Once again, this should be one line, not three
      velocity->GetPointData()->AddArray(internals.pb_velocityDoubleArray);
      vorticity->GetPointData()->AddArray(internals.pb_vorticityDoubleArray); 
      velocityNorm->GetPointData()->AddArray(internals.pb_velocityNormDoubleArray);
      */
      vtkMultiBlockDataSet *mbfluid = dynamic_cast<vtkMultiBlockDataSet*>(mesh); 
      if(!mbfluid)   
      {
        SENSEI_ERROR("unexpected mesh type "<< (mesh ? mesh->GetClassName() : "nullptr"))
        return -1;
      }
      vtkImageData *FluidImageData = (vtkImageData*)mbfluid->GetBlock(rank);
      if(!FluidImageData)
      {
        SENSEI_ERROR("Cannot Get Block in LPDataAdaptor::AddArray")
        return -1;
      }
      if(arrayName == "velocity")
      {
        FluidImageData->GetPointData()->AddArray(internals.pb_velocityDoubleArray);
      }
      else if(arrayName == "vorticity")
      {
        FluidImageData->GetPointData()->AddArray(internals.pb_vorticityDoubleArray);
      }
      else if(arrayName == "velocityNorm")
      {
        FluidImageData->GetPointData()->AddArray(internals.pb_velocityNormDoubleArray);
      }
      else
      {
        SENSEI_ERROR("Array name for Palabos AddArray does not exist in LPDataAdaptor::AddArray")
        return -1;
      }
       
    }
    return 0;
   
  }
  //----------------------------------------------------------------------
  int LPDataAdaptor::AddArrays(vtkDataObject* mesh, const std::string &meshName, int association, const std::vector<std::string> &arrayName)
  {
    return 0;
  }
  //----------------------------------------------------------------------
  int LPDataAdaptor::ReleaseData() 
  {
    //DInternals& internals = (*this->Internals); //XXX Might want to set if check depending on if you have solid or fluid data
    //internals.AtomPositions = NULL;
    //internals.anglelist = NULL;
    //internals.nanglelist = NULL;
    //internals.pb_velocityDoubleArray = NULL;
    //internals.pb_vorticityDoubleArray = NULL;
    //internals.pb_velocityNormDoubleArray = NULL;
    
    return 0;
  }

}
