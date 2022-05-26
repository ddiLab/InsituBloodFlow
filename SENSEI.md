# Building SENSEI on Cooley and running Oscillators example
1. * Clone the repo in the home directory:
`git clone https://gitlab.kitware.com/sensei/sensei.git src`

2. * Use the correct configuration of cmake while in the home directory:
  `cmake -B ./build -S ./src -DCMAKE_CXX_COMPILER=`which g++` -DCMAKE_C_COMPILER=`which gcc` -DCMAKE_INSTALL_PREFIX=./install -DENABLE_SENSEI=ON -DENABLE_VTK_IO=ON -DENABLE_CATALYST=ON -DParaView_DIR="/lus/grand/projects/visualization/srizzi/BUILDS/paraview/build/install/lib/cmake/paraview-5.9/"`

3. * Build Sensei. There should be a build directory in your home directory that was created by the cmake command:
`cd build` \
`make -j8`

4. * Install. While in the build directory:
`make install`
5. Copy oscillator executable (found in build/bin) to build/miniapps.

6. Copy oscillator.xml in src/configs to build/miniapps

7. Copy sample.osc from src/miniapps/oscillators/inputs into build/miniapps

8. Enable PosthocIO and disable any libsim in oscillator.xml 

NOTE: I copied the required files into build/miniapps, but this is not necessary. You just need to use the correct paths otherwise.

9. in build/miniapps: `export LD_LIBRARY_PATH=/lus/grand/projects/visualization/srizzi/BUILDS/paraview/build/install/lib:$LD_LIBRARY_PATH` and then `mpirun -n 2 ./oscillator -b 2 -p 0 -f oscillator.xml sample.osc` 
