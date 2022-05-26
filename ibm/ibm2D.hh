
#ifndef IBM_LBM_2D_HH
#define IBM_LBM_2D_HH

#include "atom.h"

namespace plb {
  
  template<typename T>
  void weight(T r, std::vector<T> & w){
      T q = sqrt(1 + 4*r*(1-r));
      w[0] = (3 - 2*r - q)/8.0;
      w[1] = (3 - 2*r + q)/8.0;
      w[2] = (1 + 2*r + q)/8.0;
      w[3] = (1 + 2*r - q)/8.0;
  }

  template<typename T>
  class Interpolation2D: public BoxProcessingFunctional2D_T<T,2>{
    public:
      Interpolation2D(LammpsWrapper &wrapper_):wrapper(wrapper_){}
      virtual void process(Box2D domain, TensorField2D<T,2> &velocity){
        Dot2D offset = velocity.getLocation();
        plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        T rx,ry,rz,wgt;
        Array<T,2> us(0.,0.);
        Array<T,2> uf;
        T **x = wrapper.lmp->atom->x;
        T **v = wrapper.lmp->atom->v;
        plint nlocal = wrapper.lmp->atom->nlocal;
        std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
        for (plint iS=0; iS<nlocal; iS++){
          xl = floor(x[iS][0]); 
          yl = floor(x[iS][1]); 
          //zl = floor(x[iS][2]);
          rx = x[iS][0] - xl;
          ry = x[iS][1] - yl;
          //rx = x[iS][2] - zl;
          weight<T>(rx,wx);
          weight<T>(ry,wy);
          //weight(rz,wz);
          us[0] = us[1] =0.0;
          for (ii=0;ii<4;ii++ )
            for (jj=0;jj<4;jj++ ){
                ix = xl-1 + ii - offset.x ;
                iy = yl-1 + jj - offset.y ;
                //iz = zl-1 + kk - offset.z ;
                //uf = velocity.get(ix,iy,iz);
                if (ix > domain.x1 || ix < domain.x0) continue;
                if (iy > domain.y1 || iy < domain.y0) continue;
                uf = velocity.get(ix,iy);
                wgt = wx[ii]*wy[jj];
                us[0] += wgt*uf[0];
                us[1] += wgt*uf[1];
                //us[2] += wgt*uf[2];
              }
          v[iS][0]=us[0];
          v[iS][1]=us[1];
          //v[iS][2]=us[2];
        }
      }
      virtual Interpolation2D<T> * clone() const{
        return new Interpolation2D(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::nothing; 
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;
  };


  template<typename T, template<typename U> class Descriptor>
  void interpolateVelocity2D(MultiBlockLattice2D<T,Descriptor> &lattice, MultiTensorField2D<T,2> &velocity, LammpsWrapper &wrapper)
  {/*
    plint rank = global::mpi().getRank();

    // this relies on the fact that there is exactly one block on each lattice
    plint iBlock = lattice.getLocalInfo().getBlocks()[0];
    std::map<plint,Box3D> blockmap = lattice.getSparseBlockStructure().getBulks();
    Box3D localBB = blockmap[iBlock];

    plint nx=lattice.getNx(), ny=lattice.getNy(), nz=lattice.getNz();
    plint nPart = wrapper.lmp->atom->nlocal;*/

    applyProcessingFunctional(new Interpolation2D<T>(wrapper), velocity.getBoundingBox(),velocity); 
  }

  template<typename T, template<typename U> class Descriptor>
  class Spreading2D: public BoxProcessingFunctional2D_L<T,Descriptor>{
      public:
      Spreading2D(LammpsWrapper &wrapper_):wrapper(wrapper_){}
      virtual void process(Box2D domain, BlockLattice2D<T,Descriptor> &lattice){
        Dot2D offset = lattice.getLocation();
        plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        T rx,ry,rz,wgt;
        Array<T,2> ff(0.,0.);
        //Array<T,2> fs;
        T **x = wrapper.lmp->atom->x;
        T **f = wrapper.lmp->atom->f;
        plint nlocal = wrapper.lmp->atom->nlocal;
        std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
        for (plint iS=0; iS<nlocal; iS++){
          xl = floor(x[iS][0]); 
          yl = floor(x[iS][1]); 
          //zl = floor(x[iS][2]);
          rx = x[iS][0] - xl;
          ry = x[iS][1] - yl;
          //rx = x[iS][2] - zl;
          weight<T>(rx,wx);
          weight<T>(ry,wy);
          //weight(rz,wz);
          //us[0] = us[1] =0.0;
          for (ii=0;ii<4;ii++ )
            for (jj=0;jj<4;jj++ ){
                ix = xl-1 + ii - offset.x ;
                iy = yl-1 + jj - offset.y ;
                //iz = zl-1 + kk - offset.z ;
                //uf = velocity.get(ix,iy,iz);
                if (ix > domain.x1 || ix < domain.x0) continue;
                if (iy > domain.y1 || iy < domain.y0) continue;
                Cell<T,Descriptor>& cell  = lattice.get(ix,iy);
                T *ff=cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
                wgt = wx[ii]*wy[jj];
                ff[0] += wgt*f[iS][0]; 
                ff[1] += wgt*f[iS][1]; 
                cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,ff );
                //us[0] += wgt*uf[0];
                //us[1] += wgt*uf[1];
                //us[2] += wgt*uf[2];
              }
          //v[iS][0]=us[0];
          //v[iS][1]=us[1];
          //v[iS][2]=us[2];
        }
      }
      virtual Spreading2D<T,Descriptor> * clone() const{
        return new Spreading2D(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::staticVariables; 
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;

  };

  template<typename T, template<typename U> class Descriptor>
  void spreadForce2D(MultiBlockLattice2D<T,Descriptor> &lattice,
                   LammpsWrapper &wrapper ){
    applyProcessingFunctional(new Spreading2D<T,Descriptor>(wrapper), lattice.getBoundingBox(),lattice); 
  }
  

/*
  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LiggghtsCouplingWrapper &wrapper,
                           PhysUnits3D<T> const &units,
                           bool initVelFlag)
  {
    std::vector<plint> dummyExcludeType;
    setSpheresOnLattice(lattice,wrapper,units,dummyExcludeType,initVelFlag);
  }



  template<typename T, template<typename U> class Descriptor>
  void getForcesFromLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                            LiggghtsCouplingWrapper &wrapper,
                            PhysUnits3D<T> const &units)
  {
    // debug stuff
    plint r = global::mpi().getRank();

    static std::vector<T> force,torque;
    static typename ParticleData<T>::ParticleDataArrayVector x_lb;

    plint const nPart = wrapper.lmp->atom->nlocal + wrapper.lmp->atom->nghost;
    plint const n_force = nPart*3;

    if(nPart == 0) return; // no particles - no work

    if(nPart > x_lb.size()){
      for(plint iPart=0;iPart<x_lb.size();iPart++){
        x_lb[iPart][0] = units.getLbPosition(wrapper.lmp->atom->x[iPart][0]);
        x_lb[iPart][1] = units.getLbPosition(wrapper.lmp->atom->x[iPart][1]);
        x_lb[iPart][2] = units.getLbPosition(wrapper.lmp->atom->x[iPart][2]);
      }
      for(plint iPart = x_lb.size();iPart < nPart; iPart++)
        x_lb.push_back( Array<T,3>( units.getLbPosition(wrapper.lmp->atom->x[iPart][0]),
                                    units.getLbPosition(wrapper.lmp->atom->x[iPart][1]),
                                    units.getLbPosition(wrapper.lmp->atom->x[iPart][2]) ) );
    } else{
      for(plint iPart=0;iPart<nPart;iPart++){
        x_lb[iPart][0] = units.getLbPosition(wrapper.lmp->atom->x[iPart][0]);
        x_lb[iPart][1] = units.getLbPosition(wrapper.lmp->atom->x[iPart][1]);
        x_lb[iPart][2] = units.getLbPosition(wrapper.lmp->atom->x[iPart][2]);
      }
    }

    if(n_force > force.size()){
      for(plint i=0;i<force.size();i++){
        force[i] = 0;
        torque[i] = 0;
      }
      for(plint i=force.size();i<n_force;i++){
        force.push_back(0.);
        torque.push_back(0.);
      }
    } else {
      for(plint i=0;i<n_force;i++){
        force[i] = 0;
        torque[i] = 0;
      }
    }

    SumForceTorque3D<T,Descriptor> *sft = new SumForceTorque3D<T,Descriptor>(x_lb,
                                                                             &force.front(),&torque.front(),
                                                                             wrapper);
    
    // this relies on the fact that there is exactly one block on each processor
    plint iBlock = lattice.getLocalInfo().getBlocks()[0];
    std::map<plint,Box3D> blockmap = lattice.getSparseBlockStructure().getBulks();
    Box3D localBB = blockmap[iBlock];
    applyProcessingFunctional(sft,localBB, lattice);


    LAMMPS_NS::FixLbCouplingOnetoone 
      *couplingFix 
      = dynamic_cast<LAMMPS_NS::FixLbCouplingOnetoone*>
      (wrapper.lmp->modify->find_fix_style("couple/lb/onetoone",0));

    double **f_liggghts = couplingFix->get_force_ptr();
    double **t_liggghts = couplingFix->get_torque_ptr();

    for(plint iPart=0;iPart<nPart;iPart++)
      for(plint j=0;j<3;j++){
        f_liggghts[iPart][j] = 0;
        t_liggghts[iPart][j] = 0;
      }


    for(plint iPart=0;iPart<nPart;iPart++){
      int tag = wrapper.lmp->atom->tag[iPart];
      int liggghts_ind = wrapper.lmp->atom->map(tag);

      for(plint j=0;j<3;j++){
        f_liggghts[liggghts_ind][j] += units.getPhysForce(force[3*iPart+j]);
        t_liggghts[liggghts_ind][j] += units.getPhysTorque(torque[3*iPart+j]);
      }
    }
    couplingFix->comm_force_torque();
  }
*/
}; /* namespace plb */

#endif /* IBDATAEXCHANGEWRAPPERS_HH_LBDEM */
