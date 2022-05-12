#ifndef RIB_H
#define RIB_H

#include <mesh.h>
#include <airfoil.h>
#include <iostream>
#include <fstream>

class Rib {
public:
    static const int NXI = Mesh::NXI;
    static const int NETA = Mesh::NETA;
    static const int NPTS = NXI * NETA;

    Airfoil airfoil;
    double eta[NETA];
    double camberLine[2][NXI];
    double pos[2][NXI][NETA];

    Rib(Airfoil::DV designVariables) {
        airfoil = Airfoil(designVariables);

        for (int j = 0; j < NETA;  j++) {
            eta[j] = 1.0 * j / (NETA - 1);
        }

        for (int i = 0; i < NXI; i++) {
            int oppInd = NXI-1-i;

            for (int d = 0; d < 2; d++) {
                
                //make the camberline
                camberLine[d][i] = 0.5*(airfoil.pos[d][i] + airfoil.pos[d][oppInd]);

                //make the rib coordinates
                for (int j = 0; j < NETA; j++) {
                    pos[d][i][j] = camberLine[d][i] * (1-eta[j]) + airfoil.pos[d][i] * eta[j];
                }
            }
            
        }
    }

    void printToVtk() {
        using namespace std;
        //print airfoil to vtk file        
        string sp = " ";
        string dataType = "double64";

        //make vtk files
        ofstream myfile;
        myfile.open("rib.vtk");
        myfile << "# vtk DataFile Version 2.0\n";
        myfile << "Airfoil Example\n";
        myfile << "ASCII\n";

        //make a structured grid
        myfile << "DATASET STRUCTURED_GRID\n";
        myfile << "DIMENSIONS " << NXI << sp << NETA << sp << 1 << "\n";

        //print out the points of the structured grid
        myfile << "POINTS " << NPTS << sp << dataType << "\n";

        //loop over xi the airfoil curve
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NETA; j++) {
                myfile << pos[0][i][j] << sp << pos[1][i][j] << sp << 0.0;
                myfile << "\n";
            }
            
        }

        string scalarName = "eta";
        myfile << "POINT_DATA " << NPTS << "\n";
        myfile << "SCALARS " << scalarName << " double64 1\n";
        myfile << "LOOKUP_TABLE default\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NETA; j++) {
                myfile << eta[j];
                myfile << sp;
            }
            myfile << "\n";
        }

        myfile.close();
    }


};

#endif