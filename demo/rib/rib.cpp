//cpp file to make arrays
#include <airfoil.h>
#include <rib.h>
#include <iostream>

int main() {    

    Airfoil::DV myDV;
    myDV.toverc = 0.1;
    myDV.twistAOA = 20;
    myDV.phaseShift = 0.3;
    myDV.twistMag = 0.3;
    myDV.AOA = 10.0;

    Rib myrib(myDV);
    myrib.printToVtk();
    
    return 0;
}