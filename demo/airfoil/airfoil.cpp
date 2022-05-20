//cpp file to make arrays
#include <airfoil.h>
#include <string>
#include <mesh.h>
#include <random>
#include <iostream>

//function to test airfoil derivatives
void testDerivatives() {
    //initial airfoil
    Airfoil::DV DV0;
    DV0.toverc = 0.1;
    DV0.twistAOA = 20;
    DV0.phaseShift = 0.3;
    DV0.twistMag = 0.3;
    DV0.AOA = 10.0;
    DV0.chord = 1.0;
    Airfoil airfoil0(DV0);

    static const int NXI = Mesh::NXI;

    double h = 1.0e-5;

    Airfoil::DV p;
    p.toverc = 0.01;
    p.twistAOA = 0.2;
    p.phaseShift = 0.3;
    p.twistMag = 0.5;
    p.AOA = 1.0;
    p.chord = 0.1;

    //final airfoil with perturbations
    Airfoil::DV DVf;
    DVf.toverc = DV0.toverc + p.toverc * h;
    DVf.twistAOA = DV0.twistAOA + p.twistAOA * h;
    DVf.phaseShift = DV0.phaseShift + p.phaseShift * h;
    DVf.twistMag = DV0.twistMag + p.twistMag * h;
    DVf.AOA = DV0.AOA + p.AOA * h;
    DVf.chord = DV0.chord + p.chord * h;
    Airfoil airfoilf(DVf);

    //directional derivative gradf * ph approx_eq deltaf over x+ph - x
    double pertubation[2][NXI];
    double finiteDiff, dderiv;
    finiteDiff = 0.0;
    dderiv = 0.0;
    Airfoil::DV df;
    for (int d = 0; d < 2; d++) {
        for (int i = 0; i < NXI; i++) {
            pertubation[d][i] = (double) rand() / (double) RAND_MAX; //rand double from 0 to 1
            dderiv += airfoil0.testsens[d][i].toverc * pertubation[d][i] * p.toverc * h;
            dderiv += airfoil0.testsens[d][i].twistAOA * pertubation[d][i] * p.twistAOA * h;
            dderiv += airfoil0.testsens[d][i].phaseShift * pertubation[d][i] * p.phaseShift * h;
            dderiv += airfoil0.testsens[d][i].twistMag * pertubation[d][i] * p.twistMag * h;
            dderiv += airfoil0.testsens[d][i].AOA * pertubation[d][i] * p.AOA * h;
            dderiv += airfoil0.testsens[d][i].chord * pertubation[d][i] * p.chord * h;
            finiteDiff += (airfoilf.testpos[d][i] - airfoil0.testpos[d][i]) * pertubation[d][i];

            std::cout << "dderiv = " << dderiv << " finiteDiff = " << finiteDiff << std::endl;
        }
    }
    

    //print results for toverc sens
    std::cout << "dderiv from sensitivity = " << dderiv << std::endl;
    std::cout << "finite difference deriv = " << finiteDiff << std::endl;

    //finite difference after perturb by ph


}

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
    
    testDerivatives();

    return 0;
}