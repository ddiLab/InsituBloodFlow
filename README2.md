# Using the BloodFlow Repository
Palabos + LAMMPS + SENSEI Integration for in-situ visualization is the goal for this repository.
Dr. Jifu Tan had successfully coupled Palabos and LAMMPS to simulate the flow of blood cells within plasma. 
BloodFlow has several directories, each with a specific purpose. These are listed below:
##  Quick start embolism example for Palabos and LAMMPS coupling (no SENSEI involved) NOTE: To run on Cooley rather than your personnal computer, there are several changes listed at the bottom
1. In your home directory, make a directory named *src* and change directories to *src*. \
`cd ~` \
`mkdir src` \
`cd ~/src` 
2. Clone Palabos into your *src* directory. We include a command that changes the version to the one known to work for this project. 
`git clone https://gitlab.com/unigespc/palabos.git` \
`git checkout e498e8ad7f24fd7ff87313670db7873703c1fd3f` (update to the version we know works)
3. Also clone LAMMPS to the *src* directory. \
`git clone https://github.com/lammps/lammps.git` \
`git checkout e960674cea38515ae3749218c314a9e1a3c6c140` (update to the version we know works) 
4. Make sure you have cmake (version 3.17.3), make (GNU Make 4.2.1), openmpi (Open MPI 4.0.3) installed. My tested version for each are listed in the parenthesis. Other versions should work, but these are given just in case.
5. Clone the *BloodFlow* repository. This can be done using either ssh or https. If you create a ssh key with your github account, ssh allow you to push and pull changes without having to type your username and password.
`cd ~` \
`git clone git@github.com:ddiLab/BloodFlow.git` (clone with ssh)\
`git clone https://github.com/ddiLab/BloodFlow.git` (clone with https)
6. To run both the embolism and singleCell examples, several files need to be copied from *BloodFlow/rbc* to *lammps/src*.\
`cd ~/BloodFlow/rbc` \
`cp angle* ~/src/lammps/src` \
`cp bond* ~/src/lammps/src` \
`cp dihedral* ~/src/lammps/src` \
`cp fix* ~/src/lammps/src` 

7. Go to the *src* directory that is found within the LAMMPS repository you just cloned and make sure the MOLECULE and MC packages are installed. \
`cd ~/src/lammps/src` \
`make yes-MOLECULE` \
`make yes-MC` 
8. Compile lammps as a library. Before doing so, you will have to append a `-std=c++11` flag in CCFLAGS in the ~Makefile.mpi in the *lammps/src/MAKE* folder \
`cd ~/src/lammps/src/MAKE` \
`vi Makefile.mpi` \
My CCFLAGS line looks like this: `CCFLAGS = -g -03 -std=c++11` \
Save the file. \
`cd ~/src/lammps/src` \
`make mpi mode=lib` 
9. Create a *build* directory in the *embolism* example directory and compile the embolism example. \
`cd ~/BloodFlow/examples/embolism` \
`mkdir build` \
`cd build` \
`cmake ..` \
`make -j6` \
NOTE: `-jX` make running make command faster. This depends on how many processors you can run it on.
10. Run the Simulation.
`cd ~/BloodFlow/examples/embolism` \
`mpirun -np 4 embolism in.embolism` 
---
COOLEY CHANGES: 

- The embolism.sh executable file needs to have a directory path updated. This file is in *BloodFlow/examples/embolism*. \
Change the path to: `EMB_PATH=path/to/your/embolism` 

- Several paths need to be changed in *cooley.cmake*, which is found in *BloodFlow/sites*. 

- When running on Cooley, the cmake command in step 9 should be the following: 
`cmake -C /path/to/BloodFlow/sites/cooley.cmake ..` 
- To run the simulation on Cooley, use the following command:
`qsub -n X -t X -A <ProjectName> ./embolism.sh` \
Values of X: \
-n : number of nodes \
-t : allocated time 

---

## Implementing SENSEI using singleCell example. NOTE: Do the embolism example first because this example uses the softwares built in that example. 

1. There is a file names *personal.cmake* in *BloodFlow/sites*. The first three lines need their paths corrected.
2. Make sure all these packages are installed: \
`sudo apt-get update` 
```bash
sudo apt-get install -y \
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
```
3. Clone and configure Paraview v5.9.1 superbuild in *~/src* directory. Send build files to a new directory *home/build* and set install location to a new directory *home/install*. \
`git clone --recursive https://gitlab.kitware.com/paraview/paraview.git` \
`cd paraview/`
`git checkout v5.9.1` \
`git submodule update --init --recursive`
```bash
cmake -B ~/build/paraview -S ~/src/paraview \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=~/install/paraview \
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
  ```
4. Build and install Paraview.
`cd ~build/paraview` \
`make -j8 install` 

5. Clone SENSEI in the *~/build* directory and then configure.
`git clone https://github.com/SENSEI-insitu/SENSEI.git` 
```bash
cmake -S ~/src/SENSEI -B ~/build/sensei \
  -DCMAKE_INSTALL_PREFIX=~/install/sensei \
  -DParaView_DIR=~/install/paraview/lib/cmake/paraview-5.9 \
  -DENABLE_CATALYST=ON \
  -DENABLE_VTK_IO=ON \
  -DENABLE_LAMMPS=OFF \
  -DENABLE_MANDELBROT=OFF \
  -DENABLE_OSCILLATORS=OFF \
  -DENABLE_ADIOS2=OFF \
  -DADIOS2_DIR=${ADIOS2_DIR} 
```
XXX: Do we need -DADIOS2_DIR command?

6. Build and install SENSEI \
`cd ~/build/sensei` \
`make -j8` \
`make -j8 install` 

7. Create a *build* directory in *BloodFlow/examples/singleCell/*, move into it and make comple the code. \
`cmake -DSENSEI_DIR=~/install/SENSEI/lib/cmake/ -C ~/BloodFlow/sites/personal.cmake ../` \
`make -j8` \
XXX: It should be install/SENSEI, but maybe build?
8. In the *singleCell* directory, run the simulation. \
`mpirun -n 4 cellFlow in.lmp4cell`


