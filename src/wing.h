#ifndef WING_H
#define WING_H

#include <mesh.h>
#include <airfoil.h>
#include <rib.h>
#include <iostream>
#include <fstream>
#include <cmath>

class Wing {
public:
    //number of mesh points in 
    static const int NXI = Mesh::NXI;
    static const int NETA = Mesh::NETA;
    static const int NGAMMA = Mesh::NGAMMA;

    // number of ribs and spars
    static const int NRIB = 16;
    static const int NSPAR = 2;

    //design variable inputs as structure
    struct DV {
        double dihedral, lesweep, taper, area, aspect;
        Airfoil::DV airfoilDV[2]; //chord and root design variables
        double ribFrac[NRIB-1]; //for rib design variables
        double sparFrac[NSPAR+1]; //for spar design variables
    };

    //design variable attributes from wing struct
    double dihedral, lesweep, taper, area, aspect;
    Airfoil::DV airfoilDV[2];
    Airfoil::DV tempDV;
    double ribFrac[NRIB-1];
    double sparFrac[NSPAR+1];

    //intermediate calculations
    double gamma[NGAMMA];
    double AOA[NGAMMA];
    double chord[NGAMMA];
    double pos_LE[3][NGAMMA];
    Airfoil airfoil[NGAMMA];
    //Rib rib[NRIB];
    //Spar spar[NSPAR];

    //OML, rib, and spar positions
    double OMLpos[3][NXI][NGAMMA];

    Wing(DV wingDV) {
        //copy the wing design variables
        dihedral = wingDV.dihedral;
        lesweep = wingDV.lesweep;
        taper = wingDV.taper;
        area = wingDV.area;
        aspect = wingDV.aspect;
        airfoilDV[0] = wingDV.airfoilDV[0];
        airfoilDV[1] = wingDV.airfoilDV[1];
        for (int isp = 0; isp < NRIB-1; isp++) {
            ribFrac[isp] = wingDV.ribFrac[isp];
        }
        for (int isp = 0; isp < NSPAR+1; isp++) {
            sparFrac[isp] = wingDV.sparFrac[isp];
        }

        //solve related quantities
        double cmean = pow(area/aspect, 0.5);
        double span = cmean * aspect;
        double sspan = span * 0.5;
        double croot = 2*cmean / (1+taper);
        double ctip = croot * taper;
        double xroot = 0.0;
        double yroot = 0.0;
        double zroot = 0.0;
        double xtip = sspan * tan(lesweep * M_PI/180);
        double ytip = sspan * tan(dihedral * M_PI/180);
        double ztip = -1.0 * sspan;

        //make the airfoil objects at each spanwise station with loft
        //gamma parameter from 0 to 1 along spanwise direction
        for (int k = 0; k < NGAMMA; k++) {
            gamma[k] = 1.0 * k / (NGAMMA - 1);

            //compute airfoil design variables at the spanwise station
            tempDV.toverc = airfoilDV[0].toverc * (1-gamma[k]) + airfoilDV[1].toverc * gamma[k];
            tempDV.twistAOA = airfoilDV[0].twistAOA * (1-gamma[k]) + airfoilDV[1].twistAOA * gamma[k];
            tempDV.phaseShift = airfoilDV[0].phaseShift * (1-gamma[k]) + airfoilDV[1].phaseShift * gamma[k];
            tempDV.twistMag = airfoilDV[0].twistMag * (1-gamma[k]) + airfoilDV[1].twistMag * gamma[k];
            tempDV.AOA = airfoilDV[0].AOA * (1-gamma[k]) + airfoilDV[1].AOA * gamma[k];
            tempDV.chord = croot * (1-gamma[k]) + ctip * gamma[k];

            //construct the airfoil at that spanwise station
            airfoil[k] = Airfoil(tempDV);

            //get the leading edge positions
            pos_LE[0][k] = xroot * (1-gamma[k]) + xtip * gamma[k];
            pos_LE[1][k] = yroot * (1-gamma[k]) + ytip * gamma[k];
            pos_LE[2][k] = zroot * (1-gamma[k]) + ztip * gamma[k];

            //get chord length at each station
            chord[k] = croot * (1-gamma[k]) + ctip * gamma[k];

            //make the OML positions
            for (int i = 0; i < NXI; i++) {
                for (int d = 0; d < 3; d++) {
                    OMLpos[d][i][k] = pos_LE[d][k];
                }
                for (int d = 0; d < 2; d++) {
                    OMLpos[d][i][k] += airfoil[k].pos[d][i]; 
                }
            }
        }
    }

    void printToVtk() {
        using namespace std;
        //print airfoil to vtk file        
        string sp = " ";
        string dataType = "double64";
        
        //print OML to vtk
        //make vtk files
        ofstream myfile;
        myfile.open("wing.vtk");
        myfile << "# vtk DataFile Version 2.0\n";
        myfile << "Wing Example\n";
        myfile << "ASCII\n";

        //make a structured grid
        myfile << "DATASET STRUCTURED_GRID\n";
        myfile << "DIMENSIONS " << NXI << sp << 1 << sp << NGAMMA << "\n";

        //print out the points of the structured grid
        int NPTS_OML = NXI * NGAMMA;
        myfile << "POINTS " << NPTS_OML << sp << dataType << "\n";

        //loop over xi the airfoil curve
        for (int i = 0; i < NXI; i++) {
            for (int k = 0; k < NGAMMA; k++) {
                myfile << OMLpos[0][i][k] << sp << OMLpos[1][i][k] << sp << OMLpos[2][i][k];
                myfile << "\n";
            }
            
        }

        string scalarName = "GAMMA";
        myfile << "POINT_DATA " << NPTS_OML << "\n";
        myfile << "SCALARS " << scalarName << " double64 1\n";
        myfile << "LOOKUP_TABLE default\n";
        for (int i = 0; i < NXI; i++) {
            for (int k = 0; k < NGAMMA; k++) {
                myfile << gamma[k];
                myfile << endl;
            }
        }

        myfile.close();
    }
}; //end of wing class
#endif