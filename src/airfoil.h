#ifndef AIRFOIL_H
#define AIRFOIL_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

class Airfoil {
public:
    //mesh settings
    //static const float PI = 3.14159265358979323846;
    static const int NXI = 200;
    static const int NETA = 1;
    static const int NGAMMA = 1;
    static const int NPTS = NXI * NETA * NGAMMA;

    //design variables
    float toverc, twistAOA, phaseShift, twistMag;

    //parametric curve data
    float xi[NXI];
    float sym[2][NXI];
    float tmp[2][NXI];
    float pos[2][NXI];

    Airfoil(float tc, float aoa, float beta, float mag) {
        toverc = tc;
        twistAOA = aoa;
        phaseShift = beta;
        twistMag = mag;

        //generate airfoil
        //
        //loop over each point in the spline
        for (int i = 0; i < NXI; i++) {
            //parametric curve parameter xi for x(xi), y(xi)
            xi[i] = i * 2.0 * M_PI / (NXI-1);

            //initialize symmetric airfoil
            sym[0][i] = 0.5*(cos(xi[i]) - 1);
            sym[1][i] = 2.0/3.0 * toverc * sin(xi[i]) * sin(xi[i]/2.0);
            std::cout << "step " << i << " x,y= " << sym[0][i] << ", " << sym[1][i] << "\n";

            //rotate airfoil
        }  
    }

    void printToVtk() {
        using namespace std;
        //print airfoil to vtk file        
        string sp = " ";
        string dataType = "Float64";

        //make vtk files
        ofstream myfile;
        myfile.open("airfoil.vtk");
        myfile << "# vtk DataFile Version 2.0\n";
        myfile << "Airfoil Example\n";
        myfile << "ASCII\n";

        //make a structured grid
        myfile << "DATASET STRUCTURED_GRID\n";
        myfile << "DIMENSIONS " << NXI << sp << NETA << sp << NGAMMA << "\n";

        //print out the points of the structured grid
        myfile << "POINTS " << NPTS << sp << dataType << "\n";

        //loop over xi the airfoil curve
        for (int i = 0; i < NXI; i++) {
            myfile << sym[0][i] << sp << sym[1][i] << sp << 0.0;
            myfile << "\n";
        }

        string scalarName = "xi";
        myfile << "POINT_DATA " << NPTS << "\n";
        myfile << "SCALARS " << scalarName << " Float64 1\n";
        myfile << "LOOKUP_TABLE default\n";
        for (int i = 0; i < NXI; i++) {
            myfile << xi[i];
            if (i < NXI-1) {
                myfile << "\n";
            }
        }

        myfile.close();
    }
};

#endif //AIRFOIL_H