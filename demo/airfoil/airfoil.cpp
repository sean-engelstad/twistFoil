//cpp file to make arrays
#include <airfoil.h>
#include <string>

int main() {    

    Airfoil::DV myDV;
    myDV.toverc = 0.1;
    myDV.twistAOA = 20;
    myDV.phaseShift = 0.3;
    myDV.twistMag = 0.3;
    myDV.AOA = 10.0;

    Airfoil airfoil(myDV);
    std::string scalarField = "phaseShift";
    airfoil.printToVtk(scalarField);
    
    return 0;
}