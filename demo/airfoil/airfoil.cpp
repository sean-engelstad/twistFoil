//cpp file to make arrays
#include <airfoil.h>



int main() {    

    Airfoil::DV myDV;
    myDV.toverc = 0.1;
    myDV.twistAOA = 20;
    myDV.phaseShift = 0.3;
    myDV.twistMag = 0.3;
    myDV.AOA = 10.0;

    Airfoil airfoil(myDV);
    airfoil.printToVtk();
    
    return 0;
}