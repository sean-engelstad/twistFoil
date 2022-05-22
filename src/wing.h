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
    //static const int NRIB = 16;
    //static const int NSPAR = 2;

    //design variable inputs as structure
    struct DV {
        double area = 0.0;
        double aspect = 0.0;
        double dihedral = 0.0;
        double lesweep = 0.0;
        double taper = 0.0;
        Airfoil::DV airfoilDV[2]; //chord and root design variables
        //double ribFrac[NRIB-1]; //for rib design variables
        //double sparFrac[NSPAR+1]; //for spar design variables
    };

    //design variable attributes from wing struct
    double dihedral, lesweep, taper, area, aspect;
    Airfoil::DV airfoilDV[2];
    Airfoil::DV tempDV;
    //double ribFrac[NRIB-1];
    //double sparFrac[NSPAR+1];

    //intermediate calculations
    double gamma[NGAMMA];
    double AOA[NGAMMA];
    double chord[NGAMMA];
    double pos_LE[3][NGAMMA];
    Airfoil airfoil[NGAMMA];
    //Rib rib[NRIB];
    //Spar spar[NSPAR];

    //OML, rib, and spar positions
    double testpos[3][NXI][NGAMMA];
    double OMLpos[3][NXI][NGAMMA];

    //position sensitivities
    DV testsens[3][NXI][NGAMMA];
    DV OMLsens[3][NXI][NGAMMA];

    Wing() {
        dihedral = 0.0;
    }

    Wing(DV wingDV) {
        //copy the wing design variables
        dihedral = wingDV.dihedral;
        lesweep = wingDV.lesweep;
        taper = wingDV.taper;
        area = wingDV.area;
        aspect = wingDV.aspect;
        airfoilDV[0] = wingDV.airfoilDV[0];
        airfoilDV[1] = wingDV.airfoilDV[1];

        //solve related quantities
        double cmean = pow(area/aspect, 0.5);
        double span = cmean * aspect;
        double sspan = span * 0.5;
        double croot = 2.0*cmean / (1+taper);
        double ctip = croot * taper;
        double xroot = 0.0;
        double yroot = 0.0;
        double zroot = 0.0;
        double xtip = sspan * tan(lesweep * M_PI/180);
        double ytip = sspan * tan(dihedral * M_PI/180);
        double ztip = -1.0 * sspan;
        double deg_rad = M_PI / 180.0;

        //solve intermediate quantity sensitivities
        DV cmean_sens;
        cmean_sens.area = 1.0/2.0/pow(area*aspect, 0.5);
        cmean_sens.aspect = -0.5 * pow(area,0.5) * pow(aspect, -1.5);
        DV span_sens;
        span_sens.area = cmean_sens.area * aspect;
        span_sens.aspect = 0.5 * pow(area/aspect, 0.5);
        DV sspan_sens;
        sspan_sens.area = 0.5 * span_sens.area;
        sspan_sens.aspect = 0.5 * span_sens.aspect;
        DV croot_sens;
        croot_sens.area = 2.0/(1+taper) * cmean_sens.area;
        croot_sens.aspect = 2.0/(1+taper) * cmean_sens.aspect;
        croot_sens.taper = -croot/(1+taper);
        DV ctip_sens;
        ctip_sens.area = croot_sens.area * taper;
        ctip_sens.aspect = croot_sens.aspect * taper;
        ctip_sens.taper = croot_sens.taper * taper + croot;
        DV xtip_sens, ytip_sens, ztip_sens;
        xtip_sens.area = sspan_sens.area * tan(lesweep * M_PI/180);
        xtip_sens.aspect = sspan_sens.aspect * tan(lesweep * M_PI/180);
        xtip_sens.lesweep = sspan * 1.0 / cos(lesweep * deg_rad) / cos(lesweep * M_PI/180) * deg_rad;
        ytip_sens.area = sspan_sens.area * tan(dihedral * M_PI/180);
        ytip_sens.aspect = sspan_sens.aspect * tan(dihedral * M_PI/180);
        ytip_sens.dihedral = sspan * 1.0 /cos(dihedral * deg_rad)/cos(dihedral * M_PI/180) * deg_rad;
        ztip_sens.area = -1.0 * sspan_sens.area;
        ztip_sens.aspect = -1.0 * sspan_sens.aspect;

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

            //temp DV sensitivity
            DV tempDV_toverc_sens;
            tempDV_toverc_sens.airfoilDV[0].toverc = 1-gamma[k];
            tempDV_toverc_sens.airfoilDV[1].toverc = gamma[k];
            DV tempDV_twistAOA_sens;
            tempDV_twistAOA_sens.airfoilDV[0].twistAOA = 1-gamma[k];
            tempDV_twistAOA_sens.airfoilDV[1].twistAOA = gamma[k];
            DV tempDV_phaseShift_sens;
            tempDV_phaseShift_sens.airfoilDV[0].phaseShift = 1-gamma[k];
            tempDV_phaseShift_sens.airfoilDV[1].phaseShift = gamma[k];
            DV tempDV_twistMag_sens;
            tempDV_twistMag_sens.airfoilDV[0].twistMag = 1-gamma[k];
            tempDV_twistMag_sens.airfoilDV[1].twistMag = gamma[k];
            DV tempDV_AOA_sens;
            tempDV_AOA_sens.airfoilDV[0].AOA = 1-gamma[k];
            tempDV_AOA_sens.airfoilDV[1].AOA = gamma[k];
            DV tempDV_chord_sens;
            tempDV_chord_sens.area = croot_sens.area * (1-gamma[k]) + ctip_sens.area * gamma[k];
            tempDV_chord_sens.aspect = croot_sens.aspect * (1-gamma[k]) + ctip_sens.aspect * gamma[k];
            tempDV_chord_sens.taper = croot_sens.taper * (1-gamma[k]) + ctip_sens.taper * gamma[k];

            //construct the airfoil at that spanwise station
            airfoil[k] = Airfoil(tempDV);

            //airfoil position sensitivities with respect to wing design variables
            DV airfoil_pos_sens[2][NXI];
            for (int d = 0; d < 2; d++) {
                for (int i = 0; i < NXI; i++) {
                    //airfoilDV toverc chain rule product
                    airfoil_pos_sens[d][i].airfoilDV[0].toverc = airfoil[k].sens[d][i].toverc * tempDV_toverc_sens.airfoilDV[0].toverc;
                    airfoil_pos_sens[d][i].airfoilDV[1].toverc = airfoil[k].sens[d][i].toverc * tempDV_toverc_sens.airfoilDV[1].toverc;

                    //airfoilDV twistAOA chain rule product
                    airfoil_pos_sens[d][i].airfoilDV[0].twistAOA = airfoil[k].sens[d][i].twistAOA * tempDV_twistAOA_sens.airfoilDV[0].twistAOA;
                    airfoil_pos_sens[d][i].airfoilDV[1].twistAOA = airfoil[k].sens[d][i].twistAOA * tempDV_twistAOA_sens.airfoilDV[1].twistAOA;

                    //airfoilDV phaseShift chain rule product
                    airfoil_pos_sens[d][i].airfoilDV[0].phaseShift = airfoil[k].sens[d][i].phaseShift * tempDV_phaseShift_sens.airfoilDV[0].phaseShift;
                    airfoil_pos_sens[d][i].airfoilDV[1].phaseShift = airfoil[k].sens[d][i].phaseShift * tempDV_phaseShift_sens.airfoilDV[1].phaseShift;

                    //airfoilDV twistMag chain rule product
                    airfoil_pos_sens[d][i].airfoilDV[0].twistMag = airfoil[k].sens[d][i].twistMag * tempDV_twistMag_sens.airfoilDV[0].twistMag;
                    airfoil_pos_sens[d][i].airfoilDV[1].twistMag = airfoil[k].sens[d][i].twistMag * tempDV_twistMag_sens.airfoilDV[1].twistMag;
                    
                    //airfoilDV AOA chain rule product
                    airfoil_pos_sens[d][i].airfoilDV[0].AOA = airfoil[k].sens[d][i].AOA * tempDV_AOA_sens.airfoilDV[0].AOA;
                    airfoil_pos_sens[d][i].airfoilDV[1].AOA = airfoil[k].sens[d][i].AOA * tempDV_AOA_sens.airfoilDV[1].AOA;
                    
                    //airfoilDV chord chain rule product
                    airfoil_pos_sens[d][i].area = airfoil[k].sens[d][i].chord * tempDV_chord_sens.area;
                    airfoil_pos_sens[d][i].aspect = airfoil[k].sens[d][i].chord * tempDV_chord_sens.aspect;
                    airfoil_pos_sens[d][i].taper = airfoil[k].sens[d][i].chord * tempDV_chord_sens.taper;
                }
            }

            //get the leading edge positions
            pos_LE[0][k] = xroot * (1-gamma[k]) + xtip * gamma[k];
            pos_LE[1][k] = yroot * (1-gamma[k]) + ytip * gamma[k];
            pos_LE[2][k] = zroot * (1-gamma[k]) + ztip * gamma[k];

            // //get leading edge sensitivities
            DV pos_LE_sens[3];
            pos_LE_sens[0].area = xtip_sens.area * gamma[k];
            pos_LE_sens[0].aspect = xtip_sens.aspect * gamma[k];
            pos_LE_sens[0].lesweep = xtip_sens.lesweep * gamma[k];
            //std::cout << pos_LE_sens[0].taper << " duh " << "\n";
            pos_LE_sens[1].area = ytip_sens.area * gamma[k];
            pos_LE_sens[1].aspect = ytip_sens.aspect * gamma[k];
            pos_LE_sens[1].dihedral = ytip_sens.dihedral * gamma[k];
            pos_LE_sens[2].area = ztip_sens.area * gamma[k];
            pos_LE_sens[2].aspect = ztip_sens.aspect * gamma[k];
            pos_LE_sens[2].dihedral = 0.0;
            pos_LE_sens[2].taper = 0.0;
            pos_LE_sens[2].lesweep = 0.0;

            //make the OML positions
            for (int i = 0; i < NXI; i++) {
                for (int d = 0; d < 3; d++) {
                    
                    //store OML position as leading edge point first
                    OMLpos[d][i][k] = pos_LE[d][k];

                    // //initial OML LE sensitivity
                    OMLsens[d][i][k].area = pos_LE_sens[d].area;
                    OMLsens[d][i][k].aspect = pos_LE_sens[d].aspect;
                    OMLsens[d][i][k].lesweep = pos_LE_sens[d].lesweep;
                    OMLsens[d][i][k].dihedral = pos_LE_sens[d].dihedral;
                    OMLsens[d][i][k].taper = pos_LE_sens[d].taper;

                }
                for (int d = 0; d < 2; d++) {

                    //add airfoil coordinates to OML
                    OMLpos[d][i][k] += airfoil[k].pos[d][i]; 

                    //add airfoil position sensitivities
                    OMLsens[d][i][k].area += airfoil_pos_sens[d][i].area;
                    OMLsens[d][i][k].aspect += airfoil_pos_sens[d][i].aspect;
                    OMLsens[d][i][k].taper += airfoil_pos_sens[d][i].taper;
                    for (int d2 = 0; d2 < 2; d2++) {
                        OMLsens[d][i][k].airfoilDV[d2].toverc += airfoil_pos_sens[d][i].airfoilDV[d2].toverc;
                        OMLsens[d][i][k].airfoilDV[d2].twistAOA += airfoil_pos_sens[d][i].airfoilDV[d2].twistAOA;
                        OMLsens[d][i][k].airfoilDV[d2].phaseShift += airfoil_pos_sens[d][i].airfoilDV[d2].phaseShift;
                        OMLsens[d][i][k].airfoilDV[d2].twistMag += airfoil_pos_sens[d][i].airfoilDV[d2].twistMag;
                        OMLsens[d][i][k].airfoilDV[d2].AOA += airfoil_pos_sens[d][i].airfoilDV[d2].AOA;
                        OMLsens[d][i][k].airfoilDV[d2].chord = 0.0;
                    }
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

        //make an unstructured grid even though it is really structured
        myfile << "DATASET UNSTRUCTURED_GRID\n";

        //print out the points of the structured grid
        int NPTS_OML = NXI * NGAMMA;
        myfile << "POINTS " << NPTS_OML << sp << dataType << "\n";

        //print out all xyz coordinates for each point
        for (int i = 0; i < NXI; i++) {
            for (int k = 0; k < NGAMMA; k++) {
                myfile << OMLpos[0][i][k] << sp << OMLpos[1][i][k] << sp << OMLpos[2][i][k] << "\n";
            }
        }

        //count the number of cells and the number of points
        double pointinds[NXI][NGAMMA];
        int pointind = 0;
        int ncells = 0;
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                pointinds[i][j] = pointind;
                pointind += 1;
                if (i < NXI-1 && j < NGAMMA-1) {
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
            for (int j = 0; j < NGAMMA-1; j++) {
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
        myfile << "POINT_DATA " << NPTS_OML << "\n";
        string scalarName = "AREA_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].area << sp << OMLsens[1][i][j].area << sp << OMLsens[2][i][j].area << "\n";
            }
        }

        scalarName = "ASPECT_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].aspect << sp << OMLsens[1][i][j].aspect << sp << OMLsens[2][i][j].aspect << "\n";
            }
        }

        scalarName = "DIHEDRAL_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].dihedral << sp << OMLsens[1][i][j].dihedral << sp << OMLsens[2][i][j].dihedral << "\n";
            }
        }

        scalarName = "TAPER_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].taper << sp << OMLsens[1][i][j].taper << sp << OMLsens[2][i][j].taper << "\n";
            }
        }

        scalarName = "LESWEEP_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].lesweep << sp << OMLsens[1][i][j].lesweep << sp << OMLsens[2][i][j].lesweep << "\n";
            }
        }

        scalarName = "ROOT_TOVERC_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].airfoilDV[0].toverc << sp << OMLsens[1][i][j].airfoilDV[0].toverc << sp << OMLsens[2][i][j].airfoilDV[0].toverc << "\n";
            }
        }

        scalarName = "ROOT_TWISTAOA_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].airfoilDV[0].twistAOA << sp << OMLsens[1][i][j].airfoilDV[0].twistAOA << sp << OMLsens[2][i][j].airfoilDV[0].twistAOA << "\n";
            }
        }

        scalarName = "ROOT_TWISTMAG_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].airfoilDV[0].twistMag << sp << OMLsens[1][i][j].airfoilDV[0].twistMag << sp << OMLsens[2][i][j].airfoilDV[0].twistMag << "\n";
            }
        }

        scalarName = "ROOT_AOA_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].airfoilDV[0].AOA << sp << OMLsens[1][i][j].airfoilDV[0].AOA << sp << OMLsens[2][i][j].airfoilDV[0].AOA << "\n";
            }
        }

        scalarName = "ROOT_PHASESHIFT_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].airfoilDV[0].phaseShift << sp << OMLsens[1][i][j].airfoilDV[0].phaseShift << sp << OMLsens[2][i][j].airfoilDV[0].phaseShift << "\n";
            }
        }

        scalarName = "TIP_TOVERC_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].airfoilDV[1].toverc << sp << OMLsens[1][i][j].airfoilDV[1].toverc << sp << OMLsens[2][i][j].airfoilDV[1].toverc << "\n";
            }
        }

        scalarName = "TIP_TWISTAOA_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].airfoilDV[1].twistAOA << sp << OMLsens[1][i][j].airfoilDV[1].twistAOA << sp << OMLsens[2][i][j].airfoilDV[1].twistAOA << "\n";
            }
        }

        scalarName = "TIP_TWISTMAG_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].airfoilDV[1].twistMag << sp << OMLsens[1][i][j].airfoilDV[1].twistMag << sp << OMLsens[2][i][j].airfoilDV[1].twistMag << "\n";
            }
        }

        scalarName = "TIP_AOA_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].airfoilDV[1].AOA << sp << OMLsens[1][i][j].airfoilDV[1].AOA << sp << OMLsens[2][i][j].airfoilDV[1].AOA << "\n";
            }
        }

        scalarName = "TIP_PHASESHIFT_SENS";
        myfile << "VECTORS " << scalarName << " double64\n";
        for (int i = 0; i < NXI; i++) {
            for (int j = 0; j < NGAMMA; j++) {
                myfile << OMLsens[0][i][j].airfoilDV[1].phaseShift << sp << OMLsens[1][i][j].airfoilDV[1].phaseShift << sp << OMLsens[2][i][j].airfoilDV[1].phaseShift << "\n";
            }
        }

        myfile.close();
    }
}; //end of wing class
#endif