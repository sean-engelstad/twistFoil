//cpp file to make arrays
#include <airfoil.h>
#include <wing.h>
#include <iostream>
#include <random>

double randNumber() {
    double randNumber = (double) rand() / (double) RAND_MAX;
    return randNumber;
}

//function to test airfoil derivatives
void testDerivatives() {
    int NXI = Mesh::NXI;
    int NGAMMA = Mesh::NGAMMA;

    // //make initial wing
    Wing::DV wingDV[2];
    Wing::DV p;
    Wing* wing[2];
    double h;
    //make random perturbation DV vector
    p.area = randNumber();
    p.aspect = randNumber();
    p.dihedral = randNumber();
    p.lesweep = randNumber();
    p.taper = randNumber();
    for (int d = 0; d < 2; d++) {
        p.airfoilDV[d].toverc = randNumber();
        p.airfoilDV[d].twistAOA = randNumber();
        p.airfoilDV[d].twistMag = randNumber();
        p.airfoilDV[d].phaseShift = randNumber();
        p.airfoilDV[d].AOA = randNumber();
    }

    for (int c = 0; c < 2; c++) {
        if (c == 0) {
            h = 0.0;
        } else {
            h = 1.0e-5;
        }
        wingDV[c].dihedral = 10.0 + h * p.dihedral;
        wingDV[c].aspect = 8.0+ h * p.aspect;
        wingDV[c].area = 120.0+ h * p.area;
        wingDV[c].taper = 0.5+h * p.taper;
        wingDV[c].lesweep = 10.0+ h * p.lesweep;
        
        for (int d = 0; d < 2; d++) {
            wingDV[c].airfoilDV[d].toverc = 0.1+ h * p.airfoilDV[d].toverc;
            wingDV[c].airfoilDV[d].twistAOA = 20+ h * p.airfoilDV[d].twistAOA;
            wingDV[c].airfoilDV[d].phaseShift = 0.3+ h * p.airfoilDV[d].phaseShift;
            wingDV[c].airfoilDV[d].twistMag = 0.5+ h * p.airfoilDV[d].twistMag;
            wingDV[c].airfoilDV[d].AOA = 3.0+ h * p.airfoilDV[d].AOA;
        }

        //make the new wing
        wing[c] = new Wing(wingDV[c]);
    }  

    //std::cout << wing[1]->airfoilDV[0].toverc << " this value\n";

    //directional derivative gradf * ph approx_eq deltaf over x+ph - x
    double pertubation;
    double finiteDiff, dderiv;
    finiteDiff = 0.0;
    dderiv = 0.0;
    Airfoil::DV df;
    for (int d = 0; d < 3; d++) {
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                pertubation = 1.0 * randNumber(); //rand double from 0 to 1
                //spanwise design variables
                dderiv += wing[0]->OMLsens[d][i][j].area * pertubation * p.area * h;
                dderiv += wing[0]->OMLsens[d][i][j].aspect * pertubation * p.aspect * h;
                dderiv += wing[0]->OMLsens[d][i][j].dihedral * pertubation * p.dihedral * h;
                dderiv += wing[0]->OMLsens[d][i][j].taper * pertubation * p.taper * h;
                dderiv += wing[0]->OMLsens[d][i][j].lesweep * pertubation * p.lesweep * h;
                //airfoil design variables
                for (int d2 = 0; d2 < 2; d2++) {
                    dderiv += wing[0]->OMLsens[d][i][j].airfoilDV[d2].toverc * pertubation * p.airfoilDV[d2].toverc * h;
                    dderiv += wing[0]->OMLsens[d][i][j].airfoilDV[d2].twistAOA * pertubation * p.airfoilDV[d2].twistAOA * h;
                    dderiv += wing[0]->OMLsens[d][i][j].airfoilDV[d2].phaseShift * pertubation * p.airfoilDV[d2].phaseShift * h;
                    dderiv += wing[0]->OMLsens[d][i][j].airfoilDV[d2].twistMag * pertubation * p.airfoilDV[d2].twistMag * h;
                    dderiv += wing[0]->OMLsens[d][i][j].airfoilDV[d2].AOA * pertubation * p.airfoilDV[d2].AOA * h;
                }
                finiteDiff += (wing[1]->OMLpos[d][i][j] - wing[0]->OMLpos[d][i][j]) * pertubation;
                //std::cout << "(d,i,j) = (" << d << "," << i << "," << j << ") " << finiteDiff << "\n";
                //std::cout << "span1, span2 = " << span1 << ", " << span2 << std::endl;
                //std::cout << "dderiv = " << dderiv << " finiteDiff = " << wing[1]->testpos[d][i][j] << std::endl;
            }
        }
    }
    

    //print results for toverc sens
    std::cout << "dderiv from sensitivity = " << dderiv << std::endl;
    std::cout << "finite difference deriv = " << finiteDiff << std::endl;

    //finite difference after perturb by ph


}

int main() {    

    //set the wing design variables
    Wing::DV wingDV;    
    wingDV.dihedral = 10.0;
    wingDV.aspect = 8.0;
    wingDV.area = 120.0;
    wingDV.taper = 0.5;
    wingDV.lesweep = 10.0;
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
        wingDV.airfoilDV[i] = myDV[i];
    }

    Wing mywing(wingDV);
    mywing.printToVtk();

    //testDerivatives();
    
    return 0;
}