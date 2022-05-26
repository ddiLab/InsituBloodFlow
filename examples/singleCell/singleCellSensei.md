#Running singleCell With Sensei Instructions (on Cooley)

1. Sensei must be properly built on your Cooley Account.

2. `EMB_PATH` in `BloodFlow/examples/singleCell/singleCell.sh` needs to be changed to your correct path.

3. Go to `Bloodflow/rbc` and use the following commands: `cp bond_wlc_pow.* ../../lammps/src` \ `angle_rbc.* ../../lammps/src` \ `dihedral_bend.* ../../lammps/src`

4. Make a build directory while in `BloodFlow/examples/singleCell`: `mkdir build`

5. While in the newly created build directory, run the following command (install directory is created when building Sensei): `cmake -DSENSEI_DIR=/path/to/install/lib/cmake/ -C /path/to/BloodFlow/sites/cooley.cmake ../`\
   Next run: `make -j8`

6. Return to `BloodFlow/examples/singleCell` and run: `qsub -n 1 -t 10 -q debug ./singleCell.sh`

