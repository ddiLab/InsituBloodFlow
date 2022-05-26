FROM ubuntu:18.04

SHELL ["/bin/bash", "-c"]

# directory structure for the build
ENV BASEDIR=/SeniorDesign
ENV SRCDIR=$BASEDIR/src
ENV INSTALLDIR=$BASEDIR/install
ENV BUILDDIR=$BASEDIR/build

RUN mkdir $BASEDIR
RUN mkdir $SRCDIR
RUN mkdir $INSTALLDIR
RUN mkdir $BUILDDIR

# ubuntu packages
RUN apt-get update 
RUN apt-get install -y \
    git \
    build-essential \
    autoconf \
    libtool \
    libmpich-dev \
    libssl-dev \
    wget \
    pkg-config \
    python3-dev \
    python3-numpy \
    libosmesa6-dev \
    libgl1-mesa-dev \
    libtbb-dev

# install cmake binary from Kitware
WORKDIR ${INSTALLDIR}
RUN wget https://github.com/Kitware/CMake/releases/download/v3.21.3/cmake-3.21.3-linux-x86_64.tar.gz && \
    tar xvzf cmake-3.21.3-linux-x86_64.tar.gz && \
    mv cmake-3.21.3-linux-x86_64 cmake && \
    rm cmake-3.21.3-linux-x86_64.tar.gz
ENV PATH=${INSTALLDIR}/cmake/bin:$PATH

##########################################
###         ParaView 5.9.1             ###
##########################################

WORKDIR ${SRCDIR}

# clone superbuild repo
RUN git clone --recursive https://gitlab.kitware.com/paraview/paraview.git \ 
 && cd paraview \
 && git checkout v5.9.1 \
 && git submodule update --init --recursive

# configure paraview
RUN cmake -B ${BUILDDIR}/paraview -S ${SRCDIR}/paraview \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/paraview \
  -DCMAKE_CXX_COMPILER=g++ \
  -DCMAKE_C_COMPILER=gcc \
  -DPARAVIEW_USE_PYTHON=ON \ 
  -DPARAVIEW_USE_MPI=ON \
  -DPARAVIEW_USE_QT=OFF \ 
  -DVTK_SMP_IMPLEMENTATION_TYPE=TBB \
  -DPARAVIEW_BUILD_EDITION=CATALYST_RENDERING \
  -DVTK_USE_X=OFF \
  -DVTK_OPENGL_HAS_OSMESA=ON \
  -DOSMESA_INCLUDE_DIR=/usr/include/GL/ \
  -DOSMESA_LIBRARY=/usr/lib/x86_64-linux-gnu/libOSMesa.so

# build and install
WORKDIR ${BUILDDIR}/paraview
 
run make -j8 install

##########################################
###            ADIOS 2.7.1             ###
##########################################

# clone ADIOS2
WORKDIR ${SRCDIR}
RUN git clone https://github.com/ornladios/ADIOS2.git && \
    cd ADIOS2 && \
    git checkout v2.7.1

# config ADIOS2
WORKDIR ${BASEDIR}
RUN cmake -S ${SRCDIR}/ADIOS2 -B ${BUILDDIR}/adios \
 -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/adios \
 -DADIOS2_USE_Fortran=OFF \
 -DADIOS2_BUILD_EXAMPLES=OFF \
 -DCMAKE_BUILD_TYPE=Release 

# build ADIOS2
WORKDIR ${BUILDDIR}/adios
RUN make -j8

# install ADIOS2
WORKDIR ${BUILDDIR}/adios
RUN make -j8 install

##########################################
###            SENSEI                  ###
##########################################

# clone SENSEI
WORKDIR ${SRCDIR}
# this line clones develop branch
RUN git clone https://github.com/SENSEI-insitu/SENSEI.git && \   
    cd SENSEI && \
    git checkout v3.2.2 

# config SENSEI
WORKDIR ${BASEDIR}
ENV ParaView_DIR=${INSTALLDIR}/paraview/lib/cmake/paraview-5.9
ENV ADIOS2_DIR=${INSTALLDIR}/adios/lib/cmake/adios2

RUN cmake -S ${SRCDIR}/SENSEI -B ${BUILDDIR}/sensei \
  -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/sensei \
  -DParaView_DIR=${ParaView_DIR} \
  -DENABLE_CATALYST=ON \
  -DENABLE_VTK_IO=ON \
  -DENABLE_LAMMPS=OFF \
  -DENABLE_MANDELBROT=OFF \
  -DENABLE_OSCILLATORS=OFF \
  -DENABLE_ADIOS2=OFF \
  -DADIOS2_DIR=${ADIOS2_DIR} 

# build SENSEI
WORKDIR ${BUILDDIR}/sensei
RUN make -j8

# install SENSEI
WORKDIR ${BUILDDIR}/sensei
RUN make -j8 install

###################################
###     BloodFlow repo          ###
###################################

WORKDIR ${SRCDIR}
RUN git clone https://github.com/Mectr0/BloodFlow.git

###################################
###          Palabos            ###
###################################

WORKDIR ${SRCDIR}
RUN git clone https://gitlab.com/unigespc/palabos.git && \
    cd palabos && \
    git checkout e498e8ad7f24fd7ff87313670db7873703c1fd3f

###################################
###          LAMMPS             ###
###################################

WORKDIR ${SRCDIR}
RUN git clone https://github.com/lammps/lammps.git && \
    cd lammps && \
    git checkout e960674cea38515ae3749218c314a9e1a3c6c140

###################################
###       Patch LAMMPS          ###
###################################

WORKDIR ${SRCDIR}
RUN cp ./BloodFlow/rbc/bond_wlc_pow.* ./lammps/src \
 && cp ./BloodFlow/rbc/angle_rbc.* ./lammps/src \
 && cp ./BloodFlow/rbc/dihedral_bend.* ./lammps/src \
 && cp ./BloodFlow/rbc/fix* ./lammps/src

###################################
###       Build LAMMPS          ###
###################################

WORKDIR ${SRCDIR}/lammps/src
RUN make yes-MOLECULE \
 && make yes-MC \
 && make mpi mode=lib -j 8

###################################
###      Application            ###
###################################

# configure
WORKDIR ${BASEDIR}
RUN cmake -S ${SRCDIR}/BloodFlow/examples/singleCell \
          -B ${SRCDIR}/BloodFlow/examples/singleCell/build \
          -DSENSEI_DIR=${INSTALLDIR}/sensei/lib/cmake \
          -DPALABOS_ROOT=${SRCDIR}/palabos \
          -DBLOODFLOW_ROOT=${SRCDIR}/BloodFlow \
          -DLAMMPS_DIR=${SRCDIR}/lammps \
          -DCMAKE_C_COMPILER=mpicc \
          -DCMAKE_CXX_STANDARD=11

# build
WORKDIR ${SRCDIR}/BloodFlow/examples/singleCell/build
RUN make -j 8

############################
####    Final settings   ###
############################

ENV LD_LIBRARY_PATH=${INSTALLDIR}/vtk/lib:${INSTALLDIR}/adios/lib:$LD_LIBRARY_PATH

