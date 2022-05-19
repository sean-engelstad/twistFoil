//cpp file to make arrays
#include <airfoil.h>
#include <wing.h>
#include <iostream>

int main() {    

    //set the wing design variables
    Wing::DV wingDV;
    int nrib = Wing::NRIB;
    int nspar = Wing::NSPAR;
    wingDV.dihedral = 10.0;
    wingDV.aspect = 8.0;
    wingDV.area = 120.0;
    wingDV.taper = 0.5;
    wingDV.lesweep = 10.0;
    
    int nribsp = Wing::NRIB-1;
    for (int ispace = 0; ispace < nribsp; ispace++) {
        wingDV.ribFrac[ispace] = 1.0;
    }

    int nsparsp = Wing::NSPAR+1;
    for (int ispace = 0; ispace < nsparsp; ispace++) {
        wingDV.sparFrac[ispace] = 1.0;
    }


    Airfoil::DV myDV[2];
    myDV[0].toverc = 0.1;
    myDV[0].twistAOA = 20;
    myDV[0].phaseShift = 0.3;
    myDV[0].twistMag = 0.5;
    myDV[0].AOA = 3.0;
    myDV[1].toverc = 0.1;
    myDV[1].twistAOA = 20;
    myDV[1].phaseShift = 0.3;
    myDV[1].twistMag = 0.3;
    myDV[1].AOA = -5.0;

    for (int i = 0; i < 2; i++) {
        wingDV.airfoilDV[i] = myDV[1];
    }

    Wing mywing(wingDV);
    mywing.printToVtk();
    
    return 0;
}