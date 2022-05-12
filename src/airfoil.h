#ifndef AIRFOIL_H
#define AIRFOIL_H

#include <mesh.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

class Airfoil {
public:
    static const int NXI = Mesh::NXI;
    static constexpr double KS_rho = 100000.0;

    struct DV {
        double toverc;
        double twistAOA;
        double phaseShift;
        double twistMag;
        double AOA;
        double chord = 1.0;
    };

    //design variables
    double toverc, twistAOA, phaseShift, twistMag, AOA, chord;

    double AOA_rad, twistAOA_rad;

    DV dummyDV;

    //parametric curve data
    double xi[NXI];
    double sym[2][NXI];
    double tmp[2][NXI];
    double pos[2][NXI];
    double dpos[2][NXI];
    double dotProd[NXI];
    double minVec[NXI];

    Airfoil() {
        dummyDV.toverc = 0.1;
        dummyDV.twistAOA = 0.0;
        dummyDV.phaseShift = 0.0;
        dummyDV.twistMag = 0.0;
        dummyDV.AOA = 0.0;
        dummyDV.chord = 1.0;
        //Airfoil(dummyDV);
    }

    Airfoil(DV airfoilDV) {
        toverc = airfoilDV.toverc;
        twistAOA = airfoilDV.twistAOA;
        phaseShift = airfoilDV.phaseShift;
        twistMag = airfoilDV.twistMag;
        AOA = airfoilDV.AOA; // in degrees
        chord = airfoilDV.chord;
        std::cout << "airfoil chord = " << chord << std::endl;
        //toverc *= chord;

        twistAOA_rad = twistAOA * M_PI / 180.0;
        AOA_rad = AOA * M_PI / 180.0;

        //generate symmetric airfoil, rotate and twist
        for (int i = 0; i < NXI; i++) {
            //parametric curve parameter xi for x(xi), y(xi)
            xi[i] = i * 2.0 * M_PI / (NXI-1);

            //initialize symmetric airfoil
            sym[0][i] = 0.5*(cos(xi[i]) - 1);
            sym[1][i] = 2.0/3.0 * toverc * sin(xi[i]) * sin(xi[i]/2.0);
            //std::cout << "step " << i << " x,y= " << sym[0][i] << ", " << sym[1][i] << "\n";

            //rotate airfoil
            tmp[0][i] = cos(twistAOA_rad) * sym[0][i] + sin(twistAOA_rad) * sym[1][i];
            tmp[1][i] = -sin(twistAOA_rad) * sym[0][i] + cos(twistAOA_rad) * sym[1][i];

            //twist downwards
            double argument = xi[i] - phaseShift *  sin(0.5*xi[i]);
            double LEactivation = 1 - cos(argument);
            tmp[1][i] = tmp[1][i] * (1-twistMag*LEactivation);
        }  

        //determine the leading edge and rotate down and normalize

        //determine the (dot product pos*dpos)^2
        double dxi = 2.0 * M_PI / (NXI-1);
        for (int i = 0; i < NXI; i++) {

            //compute position differentials, dpos
            for (int direc = 0; direc < 2; direc++) {
                if (i == 0) {
                    dpos[direc][0] = tmp[direc][1] - tmp[direc][0];
                } else if (i == NXI) {
                    dpos[direc][NXI] = tmp[direc][NXI] - tmp[direc][NXI-1];
                } else {
                    dpos[direc][i] = 0.5*tmp[direc][i+1] - 0.5*tmp[direc][i-1];
                }
                dpos[direc][i] /= dxi;
            }

            //std::cout << "ind = " << i << " dx = " << dpos[0][i] << " dy = " << dpos[1][i] << std::endl;
            
            //compute dot product squared (x*dx + y*dy)^2
            dotProd[i] = 0;
            for (int direc = 0; direc < 2; direc++) {
                dotProd[i] = dotProd[i] + tmp[direc][i] * dpos[direc][i];
            }
            dotProd[i] = pow(dotProd[i], 2);

            //std::cout << "(ind, dotProd^2) = (" << i << ", " << dotProd[i] << std::endl;
        }

        //use min dot product squared on left half to get leading edge
        int lowind = int(1.0/4 * NXI);
        int uppind = int(3.0/4*NXI);
        double min_vec_sum = 0.0;
        for (int i = 0; i < NXI; i++) {
            if (i < lowind || i > uppind) {
                minVec[i] = 0.0;
            } else {
                minVec[i] = exp(-KS_rho * dotProd[i]);
                min_vec_sum += minVec[i];
            }
        }

        for (int i = 0; i < NXI; i++) {
            minVec[i] = minVec[i] / min_vec_sum;
            //std::cout << "(ind, KSminvec) = (" << i << ", " << minVec[i] << ")" << std::endl;
        }

        //compute x and y coordinates of the leading edge using minvec
        double pos_LE[2] = {0.0, 0.0};
        for (int i = 0; i < NXI; i++) {
            for (int direc = 0; direc < 2; direc++) {
                pos_LE[direc] += tmp[direc][i] * minVec[i];
            }
        }

        //compute the angle of attack in radians, with x_LE negative and y_LE positive
        double AOA_LE = atan(-1.0 * pos_LE[1]/pos_LE[0]);

        //std::cout << "LE: (x,y,AOA) = (" << pos_LE[0] << ", " << pos_LE[1] << ", " << AOA_LE << ")" << std::endl;
        
        //rotate airfoil back to zero AOA
        for (int i = 0; i < NXI; i++) {
            pos[0][i] = tmp[0][i];
            pos[1][i] = tmp[1][i];

            pos[0][i] = cos(AOA_LE) * pos[0][i] - sin(AOA_LE) * pos[1][i];
            pos[1][i] = sin(AOA_LE) * pos[0][i] + cos(AOA_LE) * pos[1][i];
        }

        double pos2_LE[2];

        pos2_LE[0] = cos(AOA_LE) * pos_LE[0] - sin(AOA_LE) * pos_LE[1];
        pos2_LE[1] = sin(AOA_LE) * pos_LE[0] + cos(AOA_LE) * pos_LE[1];

        //std::cout << "LE: (x,y,AOA) = (" << pos2_LE[0] << ", " << pos2_LE[1] << ", " << AOA_LE << ")" << std::endl;

        //normalize airfoil in x-direction again
        for (int i = 0; i < NXI; i++) {
            pos[0][i] /= abs(pos2_LE[0]);
            pos[0][i] *= chord;
            pos[1][i] *= chord;
            pos[0][i] += chord;
        }

        //now shift it 1 to the left again and apply the AOA and shift it back
        double AOA_rad = AOA * M_PI / 180;
        for (int i = 0; i < NXI; i++) {
            pos[0][i] -= chord;

            pos[0][i] = cos(AOA_rad) * pos[0][i] + sin(AOA_rad) * pos[1][i];
            pos[1][i] = -sin(AOA_rad) * pos[0][i] + cos(AOA_rad) * pos[1][i];

            pos[0][i] += chord;
        }
    }

    void printToVtk() {
        using namespace std;
        //print airfoil to vtk file        
        string sp = " ";
        string dataType = "double64";

        //make vtk files
        ofstream myfile;
        myfile.open("airfoil.vtk");
        myfile << "# vtk DataFile Version 2.0\n";
        myfile << "Airfoil Example\n";
        myfile << "ASCII\n";

        //make a structured grid
        myfile << "DATASET STRUCTURED_GRID\n";
        myfile << "DIMENSIONS " << NXI << sp << 1 << sp << 1 << "\n";

        //print out the points of the structured grid
        myfile << "POINTS " << NXI << sp << dataType << "\n";

        //loop over xi the airfoil curve
        for (int i = 0; i < NXI; i++) {
            myfile << pos[0][i] << sp << pos[1][i] << sp << 0.0;
            myfile << "\n";
        }

        string scalarName = "xi";
        myfile << "POINT_DATA " << NXI << "\n";
        myfile << "SCALARS " << scalarName << " double64 1\n";
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