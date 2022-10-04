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

    //airfoil and rib curves
    Airfoil airfoil;
    double eta[NETA];
    double camberLine[2][NXI];
    double pos[2][NXI][NETA];

    //airfoil and rib curve mesh sensitivities
    Airfoil::DV camberline_sens[2][NXI];
    Airfoil::DV sens[2][NXI][NETA];

    Rib(Airfoil::DV designVariables) {
        //make airfoil curve
        airfoil = Airfoil(designVariables);
        Rib::init();
    }

    void init() {

        

        for (int j = 0; j < NETA;  j++) {

            //make the eta coordinate for rib interior
            eta[j] = 1.0 * j / (NETA - 1);

        }

        for (int i = 0; i < NXI; i++) {

            //get the index on opposite side of the airfoil
            int oppInd = NXI-1-i;

            for (int d = 0; d < 2; d++) {
                
                //make the camberline
                camberLine[d][i] = 0.5*(airfoil.pos[d][i] + airfoil.pos[d][oppInd]);

                //determine camber line sensitivity
                camberline_sens[d][i].toverc = 0.5*(airfoil.sens[d][i].toverc + airfoil.sens[d][oppInd].toverc);
                camberline_sens[d][i].twistAOA = 0.5*(airfoil.sens[d][i].twistAOA + airfoil.sens[d][oppInd].twistAOA);
                camberline_sens[d][i].phaseShift = 0.5*(airfoil.sens[d][i].phaseShift + airfoil.sens[d][oppInd].phaseShift);
                camberline_sens[d][i].twistMag = 0.5*(airfoil.sens[d][i].twistMag + airfoil.sens[d][oppInd].twistMag);
                camberline_sens[d][i].AOA = 0.5*(airfoil.sens[d][i].AOA + airfoil.sens[d][oppInd].AOA);
                camberline_sens[d][i].chord = 0.5*(airfoil.sens[d][i].chord + airfoil.sens[d][oppInd].chord);

                //set the camber

                //make the rib coordinates
                for (int j = 0; j < NETA; j++) {

                    //linear interpolate from camberline to outer airofoil curve via eta coordinate
                    pos[d][i][j] = camberLine[d][i] * (1-eta[j]) + airfoil.pos[d][i] * eta[j];

                    //sensitivity of rib positions, from linear interpolation by eta
                    sens[d][i][j].toverc = camberline_sens[d][i].toverc * (1-eta[j]) + airfoil.sens[d][i].toverc * eta[j];
                    sens[d][i][j].twistAOA = camberline_sens[d][i].twistAOA * (1-eta[j]) + airfoil.sens[d][i].twistAOA * eta[j];
                    sens[d][i][j].phaseShift = camberline_sens[d][i].phaseShift * (1-eta[j]) + airfoil.sens[d][i].phaseShift * eta[j];
                    sens[d][i][j].twistMag = camberline_sens[d][i].twistMag * (1-eta[j]) + airfoil.sens[d][i].twistMag * eta[j];
                    sens[d][i][j].AOA = camberline_sens[d][i].AOA * (1-eta[j]) + airfoil.sens[d][i].AOA * eta[j];
                    sens[d][i][j].chord = camberline_sens[d][i].chord * (1-eta[j]) + airfoil.sens[d][i].chord * eta[j];
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
        //myfile << "DATASET STRUCTURED_GRID\n";
        myfile << "DATASET UNSTRUCTURED_GRID\n";
        //myfile << "DIMENSIONS " << NXI << sp << NETA << sp << 1 << "\n";

        //print out the points of the structured grid
        myfile << "POINTS " << NPTS << sp << dataType << "\n";

        //loop over xi the airfoil curve
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NETA; j++) {
                myfile << pos[0][i][j] << sp << pos[1][i][j] << sp << 0.0;
                myfile << "\n";
            }
            
        }

        //cell data (elements)
        
        //count the number of cells and the number of points
        double pointinds[NXI][NETA];
        int pointind = 0;
        int ncells = 0;
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NETA; j++) {
                pointinds[i][j] = pointind;
                pointind += 1;
                if (i < NXI-1 && j < NETA-1) {
                    ncells += 1;
                }
            }
        }
        int ncell_pts = 5 * ncells;

        myfile << "CELLS " << ncells << " " << ncell_pts << "\n";
        int LLpt, RLpt, LRpt, RRpt;
        int cell_types[ncells];
        int icell = 0;
        for (int i = 0; i < NXI-1; i++) {
            for (int j = 0; j < NETA-1; j++) {
                LLpt = pointinds[i][j];
                RLpt = pointinds[i+1][j];
                LRpt = pointinds[i][j+1];
                RRpt = pointinds[i+1][j+1];
                myfile << 4 << sp << LLpt << sp << RLpt << sp << RRpt << sp << LRpt << "\n";
                cell_types[icell] = 9;
                icell += 1;
            }
        }

        myfile << "CELL_TYPES " << ncells << "\n";
        for (int icell = 0; icell < ncells; icell++) {
            myfile << cell_types[icell] << "\n";
        }

        //Point data for all of the vectors
        myfile << "POINT_DATA " << NPTS << "\n";
        string scalarName = "TOVERC_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NETA; j++) {
                myfile << sens[0][i][j].toverc << sp << sens[1][i][j].toverc << sp << 0.0 << "\n";
            }
        }

        scalarName = "TWISTAOA_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NETA; j++) {
                myfile << sens[0][i][j].twistAOA << sp << sens[1][i][j].twistAOA << sp << 0.0 << "\n";
            }
        }

        scalarName = "PHASESHIFT_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NETA; j++) {
                myfile << sens[0][i][j].phaseShift << sp << sens[1][i][j].phaseShift << sp << 0.0 << "\n";
            }
        }

        scalarName = "TWISTMAG_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NETA; j++) {
                myfile << sens[0][i][j].twistMag << sp << sens[1][i][j].twistMag << sp << 0.0 << "\n";
            }
        }

        scalarName = "AOA_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NETA; j++) {
                myfile << sens[0][i][j].AOA << sp << sens[1][i][j].AOA << sp << 0.0 << "\n";
            }
        }

        scalarName = "CHORD_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NETA; j++) {
                myfile << sens[0][i][j].chord << sp << sens[1][i][j].chord << sp << 0.0 << "\n";
            }
            myfile << "\n";
        }

        myfile.close();
    }


};

#endif