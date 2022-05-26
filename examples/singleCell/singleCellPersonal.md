#Running singleCell With Sensei Instructions (Personal Computer)

1. Sensei must be properly built on your computer

2. Go to `Bloodflow/rbc` and use the following commands: `cp * /path/to/lammps/src`

3. Make a build directory while in `BloodFlow/examples/singleCell`: `mkdir build`

4. There are pathes in Bloodflow/sites/personal.cmake that need to udated. *Maybe we can get this to be relative paths

5. While in the newly created build directory, run the following command (install directory is created when building Sensei): `cmake -DSENSEI_DIR=/path/to/install/lib/cmake/ -C /path/to/BloodFlow/sites/personal.cmake ../`\
Next run: `make -j8`

6. Return to `BloodFlow/examples/singleCell` and run: mpirun -n 1 cellFlow in.lmp4cell
