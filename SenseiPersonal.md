# Building SENSEI on a personnal computer

1. Clone the SENSEI repo: `git clone https://github.com/SENSEI-insitu/SENSEI.git`

2. If you don't have cmake, ccmake, and make, install them.

3. SENSEI has some dependencies on ParaView, so we will need to install the ParaView (specifically v5.9.1).

   * Clone ParaView repo with: \
   `git clone --recursive https://gitlab.kitware.com/paraview/paraview.git` \
   `cd paraview` \
   `git checkout v5.9.1` \
   `git submodule update --init --recursive` \
   
   * configure with cmake: \
  `cmake -B ${BUILDDIR}/paraview -S ${SRCDIR}/paraview \
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
  -DOSMESA_LIBRARY=/usr/lib/x86_64-linux-gnu/libOSMesa.so`

   * Build and Install:\
   `make -j8 install`

4. configure Sensei using cmake: *Here my command:
`cmake -B /home/mectro/build/sensei -S /home/mectro/src/SENSEI -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_INSTALL_PREFIX=/home/mectro/install/SENSEI -DENABLE_SENSEI=ON -DENABLE_VTK_IO=ON -DENABLE_CATALYST=ON -DParaView_DIR="/home/mectro/install/lib/cmake/paraview-5.9/"`

5. Build and Install Sensei: \
`make -j8`
`make -j8 install`


