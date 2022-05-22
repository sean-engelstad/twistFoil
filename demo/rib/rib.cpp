//cpp file to make arrays
#include <airfoil.h>
#include <rib.h>
#include <mesh.h>
#include <iostream>
#include <random>

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
    Rib rib0(DV0);

    static const int NXI = Mesh::NXI;
    static const int NETA = Mesh::NETA;

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
    Rib ribf(DVf);

    //directional derivative gradf * ph approx_eq deltaf over x+ph - x
    double pertubation[2][NXI][NETA];
    double finiteDiff, dderiv;
    finiteDiff = 0.0;
    dderiv = 0.0;
    Airfoil::DV df;
    for (int d = 0; d < 2; d++) {
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NETA; j++) {
                pertubation[d][i][j] = (double) rand() / (double) RAND_MAX; //rand double from 0 to 1
                dderiv += rib0.sens[d][i][j].toverc * pertubation[d][i][j] * p.toverc * h;
                dderiv += rib0.sens[d][i][j].twistAOA * pertubation[d][i][j] * p.twistAOA * h;
                dderiv += rib0.sens[d][i][j].phaseShift * pertubation[d][i][j] * p.phaseShift * h;
                dderiv += rib0.sens[d][i][j].twistMag * pertubation[d][i][j] * p.twistMag * h;
                dderiv += rib0.sens[d][i][j].AOA * pertubation[d][i][j] * p.AOA * h;
                dderiv += rib0.sens[d][i][j].chord * pertubation[d][i][j] * p.chord * h;
                finiteDiff += (ribf.pos[d][i][j] - rib0.pos[d][i][j]) * pertubation[d][i][j];

                //std::cout << "dderiv = " << dderiv << " finiteDiff = " << finiteDiff << std::endl;
            }
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

    Rib myrib(myDV);
    myrib.printToVtk();

    testDerivatives();
    
    return 0;
}