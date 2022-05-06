rm -f airfoil.vtk
rm -f airfoil.o
g++ -o airfoil.o airfoil.cpp
./airfoil.o
paraview airfoil.vtk