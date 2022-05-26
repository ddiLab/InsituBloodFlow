#!/bin/sh

echo "Starting Setup"
echo "Leaving BloodFlow"
cd ..
BFCDIR=$PWD # Directory Containing BloodFlow
echo "Making Directories"
mkdir insituBloodFlow
cd insituBloodFlow
mkdir install
mkdir src
mkdir build 

INSTALLDIR="$PWD/install"
BUILDDIR="$PWD/build"
SRCDIR="$PWD/src"
BASEDIR=$PWD

echo "\nchanging directory to SRC"
cd $SRCDIR
echo "\ncurrent directory is:" $PWD
echo "\ninstalling ParaView:"
git clone --recursive https://gitlab.kitware.com/paraview/paraview.git
cd paraview
git checkout v5.9.1
git submodule update --init --recursive

echo "\nchanging directory to build"
cd $BUILDDIR
echo "\nmaking paraview build directory"
mkdir paraview

echo "\nsetting ParaView install options:"
cmake -B $BUILDDIR/paraview -S $SRCDIR/paraview -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$INSTALLDIR/paraview -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DPARAVIEW_USE_PYTHON=ON -DPARAVIEW_USE_MPI=ON -DPARAVIEW_USE_QT=OFF -DVTK_SMP_IMPLEMENTATION_TYPE=TBB -DPARAVIEW_BUILD_EDITION=CATALYST_RENDERING -DVTK_USE_X=OFF -DVTK_OPENGL_HAS_OSMESA=ON -DOSMESA_INCLUDE_DIR=/usr/include/GL/ -DOSMESA_LIBRARY=/usr/lib/x86_64-linux-gnu/libOSMesa.so

echo "\nchanging directory to build/paraview"
cd paraview
echo "installing ParaView"
make -j8 install
echo "\ndone installing ParaView"

echo "\nchanging directory to src"
cd $SRCDIR
echo "\ndownloading ADIOS2:"
git clone https://github.com/ornladios/ADIOS2.git
cd ADIOS2
git checkout v2.7.1

echo "changing directory to /build/"
cd $BUILDDIR
echo "making ADIOS2 directory:"
mkdir ADIOS2
echo "setting up ADIOS2 options:"
cmake -S $SRCDIR/ADIOS2 -B $BUILDDIR/ADIOS2 -DCMAKE_INSTALL_PREFIX=$INSTALLDIR/ADIOS2 -DADIOS2_USE_Fortran=OFF -DADIOS2_BUILD_EXAMPLES=OFF -DCMAKE_BUILD_TYPE=Release
echo "changing directory to build/ADIOS2:"
cd ADIOS2
echo "installing ADIOS2:"
make -j8
make -j8 install
echo "done installing ADIOS2"

echo "changing directory to src"
cd $SRCDIR
echo "downloading SENSEI:"
git clone https://github.com/SENSEI-insitu/SENSEI.git
echo "changing to SENSEI directory"
cd SENSEI
echo "checking out a working version of SENSEI"
git checkout v3.2.2

echo "changing directory to install:"
cd $INSTALLDIR
echo "making install/SENSEI directory:"
mkdir SENSEI
echo "changing directory to build:"
cd $BUILDDIR
echo "making build/SENSEI directory:"
mkdir SENSEI
cmake -S  $SRCDIR/SENSEI -B $BUILDDIR/SENSEI -DCMAKE_INSTALL_PREFIX=$INSTALLDIR/SENSEI -DParaView_DIR=$INSTALLDIR/paraview/lib/cmake/paraview-5.9 -DENABLE_CATALYST=ON -DENABLE_VTK_IO=ON -DENABLE_LAMMPS=OFF -DENABLE_MANDELBROT=OFF -DENABLE_OSCILLATORS=OFF -DENABLE_ADIOS2=OFF -DADIOS2_DIR=$INSTALLDIR/ADIOS2/lib/cmake/adios2
echo "changing directory to build/SENSEI"
cd SENSEI
echo "installing SENSEI"
make -j8 
make -j8 install

echo "done installing SENSEI"
echo "changing directory to base"

cd $BFCDIR
echo "moving /BloodFlow directory to /insituBloodFlow:"
mv BloodFlow/ insituBloodFlow/

cd $SRCDIR
echo "downloading Palabos:"
git clone https://gitlab.com/unigespc/palabos.git

echo "changing directory to src/Palabos"
cd palabos
echo "checking out the compatible version"
git checkout e498e8ad7f24fd7ff87313670db7873703c1fd3f
echo "changing directory back to src"
cd ..

echo "downloading LAMMPS:"
git clone https://github.com/lammps/lammps.git
echo "changing directory to src/lammps"
cd lammps
echo "checking out the compatible version"
git checkout e960674cea38515ae3749218c314a9e1a3c6c140
cd ..

echo "changing directory to insituBloodFlow/BloodFlow/rbc"
cd $BASEDIR/BloodFlow/rbc
echo "copying the necessary additions to lammps library:"
cp bond_wlc_pow.* $SRCDIR/lammps/src
cp angle_rbc.* $SRCDIR/lammps/src
cp dihedral_bend.* $SRCDIR/lammps/src
cp fix* $SRCDIR/lammps/src

echo "changing directory to src/lammps/src"
cd $SRCDIR/lammps/src
echo "adding MOLECULE package to LAMMPS:"
make yes-MOLECULE
echo "adding MC package to LAMMPS:"
make yes-MC
echo "compiling LAMMPS as a library:"
make mpi mode=lib -j 8
echo "done compiling LAMMPS"

echo "changing directory to BloodFlow/examples/singleCell"
cd $BASEDIR/BloodFlow/examples/singleCell/
echo "making build directory:"
mkdir build
echo "changing directory to base"
cd $BASEDIR 
cmake -S $BASEDIR/BloodFlow/examples/singleCell -B $BASEDIR/BloodFlow/examples/singleCell/build -DSENSEI_DIR=$INSTALLDIR/SENSEI/lib/cmake -DPALABOS_ROOT=$SRCDIR/palabos -DBLOODFLOW_ROOT=$BASEDIR/BloodFlow -DLAMMPS_DIR=$SRCDIR/lammps -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_STANDARD=11


cd $BASEDIR/BloodFlow/examples/singleCell/build
make -j8
cd ..
echo "running test simulation for single cell example:"
mpirun -n 4 cellFlow in.lmp4cell