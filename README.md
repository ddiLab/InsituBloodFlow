# Using the BloodFlow Repository

Palabos + LAMMPS + SENSEI Integration for in-situ visualization is the goal for this repository.
Dr. Jifu Tan had successfully coupled Palabos and LAMMPS to simulate the flow of blood cells within plasma. 
BloodFlow has several directories and text files, each with a specific purpose. These are listed below

1. PALABOS.md : This text file includes lines of code within embolism.cpp that are of interest.
It attempts to explain code developed by Palabos that someone new to the software may not understand, even when being proficient in c++.

1. [Embolism Example](examples/embolism/README.md) lists instructions on how to build LAMMPS and Palabos as well as running the embolism example. NOTE: This example doesn't include integration of SENSEI

2. SENSEI.md : Instructions on how to build SENSEI on Cooley and how to run an example called oscillator\
   [SenseiPersonal.md](SenseiPersonal.md) gives documentation on building SENSEI on one's personal computer. This is required to run the singleCells example with SENSEI.

3. [Follow this link](examples/singleCell/singleCellPersonal.md) for instructions on how to run the singleCell simulation on a personal computer. NOTE: Going through the embolism example will help with building LAMMPS and Palabos, both of which are needed for this example


4. examples (directory) : Several examples exist here to begin getting familiar to running such code. 


5. ibm (directory) : Newer users won't need to touch this directory until they need to dig deeper into understanding how LAMMPS and Palabos are coupled.
All of the coupling code is located here. If the class for a called object in embolism.cpp can't be found in the palabos/src, it most likely is located here.

6. lmp2vtk (directory) : The example codes output a .xyz file for the particles. The code in here converts this file to .vtk, which can be visualized in Paraview. lmp2vtk.bash file has a list of commands to use that can easily convert .xyz to .vtk. 

7. rbc (directory) : Multiple files are stored in here that are needed for the example simulations. Instruction on what to do with them is located 
in the README file for each example. NOTE: Training example uses the same files as singleCell.

8. sites (directory) : cooley.cmake is basically a bash file called when using the cmake command on Cooley (listed in embolism's README). 
The paths listed in the file will need to be changed to reflect one's personnal paths. 
