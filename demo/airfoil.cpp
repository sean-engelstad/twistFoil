//cpp file to make arrays
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

int main() {
    //temp variables
    float x, x1, x2, x3;
    float y, y1, y2, y3;
    const float PI = 3.14159265358979323846;

    //airfoil design parameters
    float toverc = 0.1;
    float twistAOA = 10;
    float beta = 0.3;
    float twistMag = 0.3;
    float delta = 0.1;

    float twistAOA_rad = twistAOA * PI/180;

    //mesh spline settings
    const int NFOIL = 3;
    const int NX = 1000;
    const int NY = 10;
    const int NZ = 1;
    const int NPTS = NX*NY*NZ;
    float xi, eta, y1eta;
    struct Airfoil {
        float xi[NFOIL][NX][NY];
        float eta[NFOIL][NX][NY];
        float x[NFOIL][NX][NY];
        float y[NFOIL][NX][NY];
    } airfoil;

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
    myfile << "DIMENSIONS " << NX << sp << NY << sp << NZ << "\n";

    //choose which airfoil to plot
    int plotInd = 2;

    //print out the points of the structured grid
    myfile << "POINTS " << NPTS << sp << dataType << "\n";
    for (int i = 0; i < NX; i++) {
        xi = 2.0*PI*i/(NX-1);
        //cout << "xi = " << xi << " i = " << i << "\n";
        for (int j = 0; j < NY; j++) {
            //paramters xi,eta for mesh
            
            if (NY != 1) {
                //eta = 1+delta - 2*delta*j/(NY-1);
                eta = 1.0*j/(NY-1);
            } else {
                eta = 1.0;
            }
            airfoil.xi[0][i][j] = xi;
            airfoil.eta[0][i][j] = eta;

            //initial symmetric airfoil
            airfoil.x[0][i][j] = 0.5*(cos(xi)-1);
            airfoil.y[0][i][j] = eta * 2/3 * toverc * sin(xi) * sin(xi/2);
            x1 = airfoil.x[0][i][j];
            y1 = airfoil.y[0][i][j]; 

            //apply rotation by twistAOA
            airfoil.x[1][i][j] = cos(twistAOA_rad) * x1 + sin(twistAOA_rad) * y1;
            airfoil.y[1][i][j] = -sin(twistAOA_rad) * x1 + cos(twistAOA_rad) * y1;
            x2 = airfoil.x[1][i][j];
            y2 = airfoil.y[1][i][j];

            //twist down
            airfoil.x[2][i][j] = airfoil.x[1][i][j];
            airfoil.y[2][i][j] = airfoil.y[1][i][j] * (1 - twistMag* (1-cos(xi-beta*sin(0.5*xi))));

            //cout << "xi = " << xi << " eta = " << eta << "\n";

            //print airfoil to vtk file
            x = airfoil.x[plotInd][i][j];
            y = airfoil.y[plotInd][i][j];
            myfile << x << sp << y << sp << 0.0;
            myfile << "\n";
        }
    }

    string scalarName = "eta";
    myfile << "POINT_DATA " << NPTS << "\n";
    myfile << "SCALARS " << scalarName << " Float64 1\n";
    myfile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            myfile << airfoil.eta[0][i][j];
            if (i+j < NX+NY-2) {
                if (j < NY-1) {
                    myfile << sp;
                } else {
                  myfile << "\n";  
                }
            }
        }
    }

    myfile.close();
    return 0;

}