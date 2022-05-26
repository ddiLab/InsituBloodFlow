# How to run the embolism code

The embolism example is tested with the following versions:
- Palabos (https://gitlab.com/unigespc/palabos): \
`git clone git@gitlab.com:unigespc/palabos.git` (clone with ssh)\
`git clone https://gitlab.com/unigespc/palabos.git` (clone with https)\
`git checkout e498e8ad7f24fd7ff87313670db7873703c1fd3f` (update to the version we know works)\
- LAMMPS  (https://github.com/lammps/lammps): \
`git clone git@github.com:lammps/lammps.git` (clone with ssh)\
`git clone https://github.com/lammps/lammps.git` (clone with https)\
`git checkout e960674cea38515ae3749218c314a9e1a3c6c140` (update to the version we know works)
- Make sure you have cmake (version 3.17.3), make (GNU Make 4.2.1), openmpi (Open MPI 4.0.3) installed. My tested version for each are listed in the parenthesis. Other versions may work but I haven't tested yet.
 
1. Clone the `ddilab/BloodFlow` repository. In `BloodFlow/sites/cooley.cmake`, there are a few file paths that need to be modified. This is only required for running embolism on Cooley\
`git clone git@github.com:ddiLab/BloodFlow.git` (clone with ssh)\
`git clone https://github.com/ddiLab/BloodFlow.git` (clone with https)\
2. Embolism needs several files to be copied from `BloodFlow/rbc` to `lammps/src`. These all begin with fix. \
When in `BloodFlow/rbc`: \
 `cp fix* path/to/lammps/src`
3. The embolism.sh executable file needs to have a directory path updated. This file is in `BloodFlow/examples/embolism`. \
`EMB_PATH=path/to/your/embolism` \
**NOTE:** This step is only needed for running embolism on Cooley.\
4. Go to `lammps/src` and make sure the MOLECULE and MC packages are installed : \
`make yes-<packagename>`\
5. Compile lammps as a library, go to `lammps/src`, type `make mpi mode=lib`. You will need to append `-std=c++11` flags in CCFLAGS in Makefile.mpi under lammps/src/MAKE folder.\
6. go to `examples/embolism` and make a directory called `build`. Move to this directory and use the command: `cmake -C /path/to/BloodFlow/sites/cooley.cmake ../` and then `make`.\
**NOTE:** For running on your personnal computer the command `cmake ..` should be used. \
7. A executable named embolism should now exist in the embolism directory. Finally, use the following command while in the embolism directory to submit the code to run on a Cooley node:
`qsub -n 1 -t 10 -A <ProjectName> ./embolism.sh` (My ProjectName is visualization)
The parameters of this command can be different,this website gives the overview for qsub parameters on cooley: [Cooley Overview](https://www.alcf.anl.gov/support-center/cooley/submitting-jobs-cooley)\
**NOTE:** For your personnal computer, use a command like this when in the embolism directory: `mpirun -np 4 embolism in.embolism`


