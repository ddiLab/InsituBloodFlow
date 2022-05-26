 g++ -o lmp2vtk4s lmp2vtk4s.cpp #Creates the executable used in the next line 
 ./lmp2vtk4s -p dump.rbc.xyz -t in.cells -w vtk/cell -s cell -c 3 -z 1 #Descriptions of parameters are found in lmp2vtk4s.cpp
