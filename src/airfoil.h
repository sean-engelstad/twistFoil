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

    double AOA_rad, twistAOA_rad;

    //parametric curve data
    double xi[NXI];
    double sym[2][NXI];
    double tmp[2][NXI];
    double pos[2][NXI];
    double dpos[2][NXI];
    double dotProd[NXI];
    double minVec[NXI];

    //sensitivity data
    DV sens[2][NXI];

    //derivative testing data
    double testpos[2][NXI];
    DV testsens[2][NXI];

    Airfoil (double vector[]) {
        Airfoil::DV airfoilDV;
        airfoilDV.toverc = vector[0];
        airfoilDV.twistAOA = vector[1];
        airfoilDV.phaseShift = vector[2];
        airfoilDV.twistMag = vector[3];
        airfoilDV.AOA = vector[4];
        airfoilDV.chord = vector[5];

        Airfoil::init(airfoilDV);
    }

    Airfoil(DV airfoilDV) {
        Airfoil::init(airfoilDV);
    }

    void init(DV airfoilDV) {
        //initialize airfoil sensitivity
        for (int d = 0; d < 2; d++) {
            for (int i = 0; i < NXI; i++) {
                sens[d][i].AOA = 0.0;
                sens[d][i].chord = 0.0;
                sens[d][i].phaseShift = 0.0;
                sens[d][i].toverc = 0.0;
                sens[d][i].twistAOA = 0.0;
                sens[d][i].twistMag = 0.0;
            }
        }
        
        //convert angles from deg to radians
        twistAOA_rad = airfoilDV.twistAOA * M_PI / 180.0;
        AOA_rad = airfoilDV.AOA * M_PI / 180.0;

        double deg_rad_sens = M_PI / 180.0;

        //generate symmetric airfoil, rotate and twist
        for (int i = 0; i < NXI; i++) {
            //parametric curve parameter xi for x(xi), y(xi)
            xi[i] = i * 2.0 * M_PI / (NXI-1);

            //initialize symmetric airfoil
            sym[0][i] = 0.5*(cos(xi[i]) - 1);
            sym[1][i] = 2.0/3.0 * airfoilDV.toverc * sin(xi[i]) * sin(xi[i]/2.0);
            //std::cout << "step " << i << " x,y= " << sym[0][i] << ", " << sym[1][i] << "\n";
            
            //record airfoil sensitivity to toverc in y direction
            sens[0][i].toverc = 0.0;
            sens[1][i].toverc = 2.0/3.0 * sin(xi[i]) * sin(xi[i]/2.0);

            //rotate airfoil
            tmp[0][i] = cos(twistAOA_rad) * sym[0][i] + sin(twistAOA_rad) * sym[1][i];
            tmp[1][i] = -sin(twistAOA_rad) * sym[0][i] + cos(twistAOA_rad) * sym[1][i];

            //record airfoil sensitivity to rotation and propagate toverc
            sens[0][i].toverc = sin(twistAOA_rad) * sens[1][i].toverc;
            sens[1][i].toverc = cos(twistAOA_rad) * sens[1][i].toverc;
            sens[0][i].twistAOA = -sin(twistAOA_rad) * deg_rad_sens * sym[0][i] + cos(twistAOA_rad) * deg_rad_sens * sym[1][i];
            sens[1][i].twistAOA = -cos(twistAOA_rad) * deg_rad_sens * sym[0][i] - sin(twistAOA_rad) * deg_rad_sens * sym[1][i];

            //twist downwards quantities
            double argument = xi[i] - airfoilDV.phaseShift *  sin(0.5*xi[i]);
            double LEactivation = 1 - cos(argument);
            double twistFunction = (1-airfoilDV.twistMag*LEactivation);

            //record new sensitivities to twist
            double arg_phaseshift_sens = -sin(0.5 * xi[i]);
            double LEactivation_phaseshift_sens = sin(argument) * arg_phaseshift_sens;
            double twistFunction_phaseshift_sens = -airfoilDV.twistMag * LEactivation_phaseshift_sens;
            double twistFunction_twistMag_sens = -LEactivation;
            sens[1][i].phaseShift = tmp[1][i] * twistFunction_phaseshift_sens;
            sens[1][i].twistMag = tmp[1][i] * twistFunction_twistMag_sens;

            //apply the twist downwards
            tmp[1][i] *= twistFunction;

            //update previous sensitivities during twist
            sens[1][i].toverc *= twistFunction;
            sens[1][i].twistAOA *= twistFunction;

        }  

        //initialize dsens which represents sensitivity of dpos
        DV dsens[2][NXI];
        for (int d = 0; d < 2; d++) {
            for (int i = 0; i < NXI; i++) {
                dsens[d][i].AOA = 0.0;
                dsens[d][i].chord = 0.0;
                dsens[d][i].phaseShift = 0.0;
                dsens[d][i].toverc = 0.0;
                dsens[d][i].twistAOA = 0.0;
                dsens[d][i].twistMag = 0.0;
            }
        }

        //determine the (dot product pos*dpos)^2
        double dxi = 2.0 * M_PI / (NXI-1);
        //initialize dot product sensitivity
        DV dotprod_sens[NXI];

        for (int i = 0; i < NXI; i++) {

            for (int direc = 0; direc < 2; direc++) {

                //compute position differentials, dpos
                if (i == 0) {
                    dpos[direc][0] = tmp[direc][1] - tmp[direc][0];
                } else if (i == NXI-1) {
                    dpos[direc][NXI-1] = tmp[direc][NXI-1] - tmp[direc][NXI-2];
                } else {
                    dpos[direc][i] = 0.5*tmp[direc][i+1] - 0.5*tmp[direc][i-1];
                }
                dpos[direc][i] /= dxi;

                //record differential position sensitivities, dsens
                //current list of sens args: toverc, twistAOA, phaseShift, twistMag
                if (i == 0) {
                    dsens[direc][0].toverc = sens[direc][1].toverc - sens[direc][0].toverc;
                    dsens[direc][0].twistAOA = sens[direc][1].twistAOA - sens[direc][0].twistAOA;
                    dsens[direc][0].phaseShift = sens[direc][1].phaseShift - sens[direc][0].phaseShift;
                    dsens[direc][0].twistMag = sens[direc][1].twistMag - sens[direc][0].twistMag;
                } else if (i == NXI-1) {
                    dsens[direc][NXI-1].toverc = sens[direc][NXI-1].toverc - sens[direc][NXI-2].toverc;
                    dsens[direc][NXI-1].twistAOA = sens[direc][NXI-1].twistAOA - sens[direc][NXI-2].twistAOA;
                    dsens[direc][NXI-1].phaseShift = sens[direc][NXI-1].phaseShift - sens[direc][NXI-2].phaseShift;
                    dsens[direc][NXI-1].twistMag = sens[direc][NXI-1].twistMag - sens[direc][NXI-2].twistMag;
                } else {
                    dsens[direc][i].toverc = 0.5 * sens[direc][i+1].toverc - 0.5 * sens[direc][i-1].toverc;
                    dsens[direc][i].twistAOA = 0.5 * sens[direc][i+1].twistAOA - 0.5 * sens[direc][i-1].twistAOA;
                    dsens[direc][i].phaseShift = 0.5 * sens[direc][i+1].phaseShift - 0.5 * sens[direc][i-1].phaseShift;
                    dsens[direc][i].twistMag = 0.5 * sens[direc][i+1].twistMag - 0.5 * sens[direc][i-1].twistMag;
                }
                dsens[direc][i].toverc /= dxi;
                dsens[direc][i].twistAOA /= dxi;
                dsens[direc][i].phaseShift /= dxi;
                dsens[direc][i].twistMag /= dxi;
            }

            //std::cout << "ind = " << i << " dx = " << dpos[0][i] << " dy = " << dpos[1][i] << std::endl;
            
            //compute dot product squared (x*dx + y*dy)^2
            dotProd[i] = 0.0;
            for (int direc = 0; direc < 2; direc++) {
                dotProd[i] += tmp[direc][i] * dpos[direc][i];
            }
            dotProd[i] *= dotProd[i]; //square it

            //initialize dot product sensitivity
            dotprod_sens[i].AOA = 0.0;
            dotprod_sens[i].chord = 0.0;
            dotprod_sens[i].phaseShift = 0.0;
            dotprod_sens[i].toverc = 0.0;
            dotprod_sens[i].twistAOA = 0.0;
            dotprod_sens[i].twistMag = 0.0;

            //compute dot product sensitivity
            //current list of sens args: toverc, twistAOA, phaseShift, twistMag
            for (int direc = 0; direc < 2; direc++) {
                dotprod_sens[i].toverc += sens[direc][i].toverc * dpos[direc][i] + tmp[direc][i] * dsens[direc][i].toverc;
                dotprod_sens[i].twistAOA += sens[direc][i].twistAOA * dpos[direc][i] + tmp[direc][i] * dsens[direc][i].twistAOA;
                dotprod_sens[i].phaseShift += sens[direc][i].phaseShift * dpos[direc][i] + tmp[direc][i] * dsens[direc][i].phaseShift;
                dotprod_sens[i].twistMag += sens[direc][i].twistMag * dpos[direc][i] + tmp[direc][i] * dsens[direc][i].twistMag;
            }
            dotprod_sens[i].toverc *= 2*dotProd[i];
            dotprod_sens[i].twistAOA *= 2*dotProd[i];
            dotprod_sens[i].phaseShift *= 2*dotProd[i];
            dotprod_sens[i].twistMag *= 2*dotProd[i];
            //std::cout << "(ind, dotProd^2) = (" << i << ", " << dotProd[i] << std::endl;

            //current state of derivative verification, stuck here and dotProd sens not matching
            testpos[0][i] = 1.0 * dotProd[i];
            testsens[0][i].toverc = 1.0 * dotprod_sens[i].toverc;
            testsens[0][i].twistAOA = 1.0 * dotprod_sens[i].twistAOA;
            testsens[0][i].AOA = 1.0 * dotprod_sens[i].AOA;
            testsens[0][i].chord = 1.0 * dotprod_sens[i].chord;
            testsens[0][i].phaseShift = 1.0 * dotprod_sens[i].phaseShift;
            testsens[0][i].twistMag = 1.0 * dotprod_sens[i].twistMag;
            testpos[1][i] = 0.0;
            testsens[1][i].toverc = 0.0;
            testsens[1][i].twistAOA = 0.0;
            testsens[1][i].AOA = 0.0;
            testsens[1][i].chord = 0.0;
            testsens[1][i].phaseShift = 0.0;
            testsens[1][i].twistMag = 0.0;
        }

        // for (int d = 0; d < 2; d++) {
        //     testpos[d][i] = 1.0 * dpos[d][i];
        //     testsens[d][i].toverc = 1.0 * dsens[d][i].toverc;
        //     testsens[d][i].twistAOA = 1.0 * dsens[d][i].twistAOA;
        //     testsens[d][i].AOA = 1.0 * dsens[d][i].AOA;
        //     testsens[d][i].chord = 1.0 * dsens[d][i].chord;
        //     testsens[d][i].phaseShift = 1.0 * dsens[d][i].phaseShift;
        //     testsens[d][i].twistMag = 1.0 * dsens[d][i].twistMag;
        //     //std::cout << "direction " << d << " point " << i << " sensitivity " << sens[d][i].toverc << std::endl;
        // }

        //use min dot product squared on left half to get leading edge
        int lowind = int(1.0/4 * NXI);
        int uppind = int(3.0/4*NXI);
        double min_vec_sum = 0.0;

        //initialize minvec su[m sensitivity
        DV minVec_sens[NXI];
        DV minVec_sum_sens;
        minVec_sum_sens.toverc = 0.0;
        minVec_sum_sens.twistAOA = 0.0;
        minVec_sum_sens.phaseShift = 0.0;
        minVec_sum_sens.twistMag = 0.0;

        for (int i = 0; i < NXI; i++) {

            //compute initial minvec before normalization
            if (i < lowind || i > uppind) {
                minVec[i] = 0.0;
            } else {
                minVec[i] = exp(-KS_rho * dotProd[i]);
                min_vec_sum += minVec[i];
            }

            //initialize minvec sensitivities
            //current list of sens args: toverc, twistAOA, phaseShift, twistMag
            minVec_sens[i].toverc = 0.0;
            minVec_sens[i].twistAOA = 0.0;
            minVec_sens[i].phaseShift = 0.0;
            minVec_sens[i].twistMag = 0.0;

            //initialize minvec sum sens

            //compute min vec sensitivities before normalization
            if (i < lowind || i > uppind) {
                minVec_sens[i].toverc = 0.0;
                minVec_sens[i].twistAOA = 0.0;
                minVec_sens[i].phaseShift = 0.0;
                minVec_sens[i].twistMag = 0.0;
            } else {
                //update minvec sens
                minVec_sens[i].toverc = exp(-KS_rho * dotProd[i]) * dotprod_sens[i].toverc;
                minVec_sens[i].twistAOA = exp(-KS_rho * dotProd[i]) * dotprod_sens[i].twistAOA;
                minVec_sens[i].phaseShift = exp(-KS_rho * dotProd[i]) * dotprod_sens[i].phaseShift;
                minVec_sens[i].twistMag = exp(-KS_rho * dotProd[i]) * dotprod_sens[i].twistMag;

                //update minvec sum sens
                minVec_sum_sens.toverc += minVec_sens[i].toverc;
                minVec_sum_sens.twistAOA += minVec_sens[i].twistAOA;
                minVec_sum_sens.phaseShift += minVec_sens[i].phaseShift;
                minVec_sum_sens.twistMag += minVec_sens[i].twistMag;
            }
        }

        for (int i = 0; i < NXI; i++) {

            //normalize minvec so total sum is 1
            minVec[i] = minVec[i] / min_vec_sum;

            //update minvec sensitivity with quotient rule
            //current list of sens args: toverc, twistAOA, phaseShift, twistMag
            minVec_sens[i].toverc = minVec_sens[i].toverc / min_vec_sum - minVec[i] * minVec_sum_sens.toverc / min_vec_sum / min_vec_sum;
            minVec_sens[i].twistAOA = minVec_sens[i].twistAOA / min_vec_sum - minVec[i] * minVec_sum_sens.twistAOA / min_vec_sum / min_vec_sum;
            minVec_sens[i].phaseShift = minVec_sens[i].phaseShift / min_vec_sum - minVec[i] * minVec_sum_sens.phaseShift / min_vec_sum / min_vec_sum;
            minVec_sens[i].twistMag = minVec_sens[i].twistMag / min_vec_sum - minVec[i] * minVec_sum_sens.twistMag / min_vec_sum / min_vec_sum;

            //std::cout << "(ind, KSminvec) = (" << i << ", " << minVec[i] << ")" << std::endl;
        }

        //compute x and y coordinates of the leading edge using minvec
        double pos_LE[2] = {0.0, 0.0};

        //initialize sensitivity of position of leading edge
        DV pos_LE_sens[2];
        for (int d = 0; d < 2; d++) {
            pos_LE_sens[d].toverc = 0.0;
            pos_LE_sens[d].twistAOA = 0.0;
            pos_LE_sens[d].phaseShift = 0.0;
            pos_LE_sens[d].twistMag = 0.0;
        }

        for (int i = 0; i < NXI; i++) {
            for (int direc = 0; direc < 2; direc++) {

                //dot product current position tmp with minVec unit vector
                pos_LE[direc] += tmp[direc][i] * minVec[i];

                //compute pos_LE sensitivities with product rule
                pos_LE_sens[direc].toverc += sens[direc][i].toverc * minVec[i] + tmp[direc][i] * minVec_sens[i].toverc;
                pos_LE_sens[direc].twistAOA += sens[direc][i].twistAOA * minVec[i] + tmp[direc][i] * minVec_sens[i].twistAOA;
                pos_LE_sens[direc].phaseShift += sens[direc][i].phaseShift * minVec[i] + tmp[direc][i] * minVec_sens[i].phaseShift;
                pos_LE_sens[direc].twistMag += sens[direc][i].twistMag * minVec[i] + tmp[direc][i] * minVec_sens[i].twistMag;
            }
        }

        //compute the angle of attack in radians, with x_LE negative and y_LE positive
        double argument = -1.0 * pos_LE[1]/pos_LE[0];
        double AOA_LE = atan(argument);

        //compute sensitivity of argument for AOA_LE
        DV argument_sens;
        argument_sens.toverc = -1.0 * pos_LE_sens[1].toverc / pos_LE[0] + pos_LE_sens[0].toverc * pos_LE[1] / pos_LE[0] / pos_LE[0];
        argument_sens.twistAOA = -1.0 * pos_LE_sens[1].twistAOA / pos_LE[0] + pos_LE_sens[0].twistAOA * pos_LE[1] / pos_LE[0] / pos_LE[0];
        argument_sens.phaseShift = -1.0 * pos_LE_sens[1].phaseShift / pos_LE[0] + pos_LE_sens[0].phaseShift * pos_LE[1] / pos_LE[0] / pos_LE[0];
        argument_sens.twistMag = -1.0 * pos_LE_sens[1].twistMag / pos_LE[0] + pos_LE_sens[0].twistMag * pos_LE[1] / pos_LE[0] / pos_LE[0];

        //compute sensitivity for AOA_LE
        DV AOA_LE_sens;
        AOA_LE_sens.toverc = 1.0/(1+argument*argument) * argument_sens.toverc;
        AOA_LE_sens.twistAOA = 1.0/(1+argument*argument) * argument_sens.twistAOA;
        AOA_LE_sens.phaseShift = 1.0/(1+argument*argument) * argument_sens.phaseShift;
        AOA_LE_sens.twistMag = 1.0/(1+argument*argument) * argument_sens.twistMag;


        //std::cout << "LE: (x,y,AOA) = (" << pos_LE[0] << ", " << pos_LE[1] << ", " << AOA_LE << ")" << std::endl;
        
        //rotate airfoil back to zero AOA
        for (int i = 0; i < NXI; i++) {

            //copy tmp positions to final airfoil positions
            pos[0][i] = tmp[0][i];
            pos[1][i] = tmp[1][i];

            //rotate airfoil back to zero angle of attack with AOA_LE
            pos[0][i] = cos(AOA_LE) * pos[0][i] - sin(AOA_LE) * pos[1][i];
            pos[1][i] = sin(AOA_LE) * pos[0][i] + cos(AOA_LE) * pos[1][i];

            //update sens with AOA_LE transformation
            sens[0][i].toverc = -sin(AOA_LE) * AOA_LE_sens.toverc * pos[0][i] + cos(AOA_LE) * sens[0][i].toverc -\
            cos(AOA_LE) * AOA_LE_sens.toverc * pos[1][i] - sin(AOA_LE) * sens[1][i].toverc;
            sens[1][i].toverc = cos(AOA_LE) * AOA_LE_sens.toverc * pos[0][i] + sin(AOA_LE) * sens[0][i].toverc -\
            sin(AOA_LE) * AOA_LE_sens.toverc * pos[1][i] + cos(AOA_LE) * sens[1][i].toverc;

            sens[0][i].twistAOA = -sin(AOA_LE) * AOA_LE_sens.twistAOA * pos[0][i] + cos(AOA_LE) * sens[0][i].twistAOA -\
            cos(AOA_LE) * AOA_LE_sens.twistAOA * pos[1][i] - sin(AOA_LE) * sens[1][i].twistAOA;
            sens[1][i].twistAOA = cos(AOA_LE) * AOA_LE_sens.twistAOA * pos[0][i] + sin(AOA_LE) * sens[0][i].twistAOA -\
            sin(AOA_LE) * AOA_LE_sens.twistAOA * pos[1][i] + cos(AOA_LE) * sens[1][i].twistAOA;

            sens[0][i].phaseShift = -sin(AOA_LE) * AOA_LE_sens.phaseShift * pos[0][i] + cos(AOA_LE) * sens[0][i].phaseShift -\
            cos(AOA_LE) * AOA_LE_sens.phaseShift * pos[1][i] - sin(AOA_LE) * sens[1][i].phaseShift;
            sens[1][i].phaseShift = cos(AOA_LE) * AOA_LE_sens.phaseShift * pos[0][i] + sin(AOA_LE) * sens[0][i].phaseShift -\
            sin(AOA_LE) * AOA_LE_sens.phaseShift * pos[1][i] + cos(AOA_LE) * sens[1][i].phaseShift;

            sens[0][i].twistMag = -sin(AOA_LE) * AOA_LE_sens.twistMag * pos[0][i] + cos(AOA_LE) * sens[0][i].twistMag -\
            cos(AOA_LE) * AOA_LE_sens.twistMag * pos[1][i] - sin(AOA_LE) * sens[1][i].twistMag;
            sens[1][i].twistMag = cos(AOA_LE) * AOA_LE_sens.twistMag * pos[0][i] + sin(AOA_LE) * sens[0][i].twistMag -\
            sin(AOA_LE) * AOA_LE_sens.twistMag * pos[1][i] + cos(AOA_LE) * sens[1][i].twistMag;
        }

        //compute pos_LE after rotation
        double pos2_LE[2];
        pos2_LE[0] = cos(AOA_LE) * pos_LE[0] - sin(AOA_LE) * pos_LE[1];
        pos2_LE[1] = sin(AOA_LE) * pos_LE[0] + cos(AOA_LE) * pos_LE[1];

        //compute sensitivities for updated pos_LE after rotation
        DV pos2_LE_sens[2];
        for (int d = 0; d < 2; d++) {
            pos2_LE_sens[d].toverc = 0.0;
            pos2_LE_sens[d].twistAOA = 0.0;
            pos2_LE_sens[d].phaseShift = 0.0;
            pos2_LE_sens[d].twistMag = 0.0;
        }

        pos2_LE_sens[0].toverc = cos(AOA_LE) * pos_LE_sens[0].toverc - sin(AOA_LE) * pos_LE_sens[1].toverc -\
        sin(AOA_LE) * AOA_LE_sens.toverc * pos_LE[0] - cos(AOA_LE) * AOA_LE_sens.toverc * pos_LE[1];
        pos2_LE_sens[1].toverc = sin(AOA_LE) * pos_LE_sens[0].toverc + cos(AOA_LE) * pos_LE_sens[1].toverc +\
        cos(AOA_LE) * AOA_LE_sens.toverc * pos_LE[0] - sin(AOA_LE) * AOA_LE_sens.toverc * pos_LE[1];

        pos2_LE_sens[0].twistAOA = cos(AOA_LE) * pos_LE_sens[0].twistAOA - sin(AOA_LE) * pos_LE_sens[1].twistAOA -\
        sin(AOA_LE) * AOA_LE_sens.twistAOA * pos_LE[0] - cos(AOA_LE) * AOA_LE_sens.twistAOA * pos_LE[1];
        pos2_LE_sens[1].twistAOA = sin(AOA_LE) * pos_LE_sens[0].twistAOA + cos(AOA_LE) * pos_LE_sens[1].twistAOA +\
        cos(AOA_LE) * AOA_LE_sens.twistAOA * pos_LE[0] - sin(AOA_LE) * AOA_LE_sens.twistAOA * pos_LE[1];

        pos2_LE_sens[0].phaseShift = cos(AOA_LE) * pos_LE_sens[0].phaseShift - sin(AOA_LE) * pos_LE_sens[1].phaseShift -\
        sin(AOA_LE) * AOA_LE_sens.phaseShift * pos_LE[0] - cos(AOA_LE) * AOA_LE_sens.phaseShift * pos_LE[1];
        pos2_LE_sens[1].phaseShift = sin(AOA_LE) * pos_LE_sens[0].phaseShift + cos(AOA_LE) * pos_LE_sens[1].phaseShift +\
        cos(AOA_LE) * AOA_LE_sens.phaseShift * pos_LE[0] - sin(AOA_LE) * AOA_LE_sens.phaseShift * pos_LE[1];

        pos2_LE_sens[0].twistMag = cos(AOA_LE) * pos_LE_sens[0].twistMag - sin(AOA_LE) * pos_LE_sens[1].twistMag -\
        sin(AOA_LE) * AOA_LE_sens.twistMag * pos_LE[0] - cos(AOA_LE) * AOA_LE_sens.twistMag * pos_LE[1];
        pos2_LE_sens[1].twistMag = sin(AOA_LE) * pos_LE_sens[0].twistMag + cos(AOA_LE) * pos_LE_sens[1].twistMag +\
        cos(AOA_LE) * AOA_LE_sens.twistMag * pos_LE[0] - sin(AOA_LE) * AOA_LE_sens.twistMag * pos_LE[1];


        //std::cout << "LE: (x,y,AOA) = (" << pos2_LE[0] << ", " << pos2_LE[1] << ", " << AOA_LE << ")" << std::endl;

        //normalize airfoil in x-direction again
        for (int i = 0; i < NXI; i++) {

            //normalize airfoil x position
            pos[0][i] /= -1.0 * pos2_LE[0];

            //update sensitivities for normalize by x position
            sens[0][i].toverc = sens[0][i].toverc * -1.0 / pos2_LE[0] + pos[0][i] * pos2_LE_sens[0].toverc / pos2_LE[0] / pos2_LE[0];
            sens[0][i].twistAOA = sens[0][i].twistAOA * -1.0 / pos2_LE[0] + pos[0][i] * pos2_LE_sens[0].twistAOA / pos2_LE[0] / pos2_LE[0];
            sens[0][i].phaseShift = sens[0][i].phaseShift * -1.0 / pos2_LE[0] + pos[0][i] * pos2_LE_sens[0].phaseShift / pos2_LE[0] / pos2_LE[0];
            sens[0][i].twistMag = sens[0][i].twistMag * -1.0 / pos2_LE[0] + pos[0][i] * pos2_LE_sens[0].twistMag / pos2_LE[0] / pos2_LE[0];

            for (int d = 0; d < 2; d++) {

                //scale up position by chord size
                pos[d][i] *= airfoilDV.chord;

                //update sensitivities for scaling up chord size
                sens[d][i].toverc *= airfoilDV.chord;
                sens[d][i].twistAOA *= airfoilDV.chord;
                sens[d][i].phaseShift *= airfoilDV.chord;
                sens[d][i].twistMag *= airfoilDV.chord;
                sens[d][i].chord = pos[d][i];
            }

            //shift airfoil so LE at origin
            pos[0][i] += airfoilDV.chord;

            //update sensitivities for chord shift
            sens[0][i].chord += 1.0;
        }

        //now shift it 1 to the left again and apply the AOA and shift it back
        double AOA_rad = airfoilDV.AOA * M_PI / 180;
        double AOA_rad_sens = M_PI / 180;

        for (int i = 0; i < NXI; i++) {

            //shift airfoil to left again
            pos[0][i] -= airfoilDV.chord;
            sens[0][i].chord -= 1.0;

            //rotate airfoil by normal AOA
            pos[0][i] = cos(AOA_rad) * pos[0][i] + sin(AOA_rad) * pos[1][i];
            pos[1][i] = -sin(AOA_rad) * pos[0][i] + cos(AOA_rad) * pos[1][i];

            //update the sensitivities for AOA rotation
            sens[0][i].toverc = cos(AOA_rad) * sens[0][i].toverc + sin(AOA_rad) * sens[0][i].toverc;
            sens[1][i].toverc = -sin(AOA_rad) * sens[0][i].toverc + cos(AOA_rad) * sens[1][i].toverc;
            
            sens[0][i].twistAOA = cos(AOA_rad) * sens[0][i].twistAOA + sin(AOA_rad) * sens[0][i].twistAOA;
            sens[1][i].twistAOA = -sin(AOA_rad) * sens[0][i].twistAOA + cos(AOA_rad) * sens[1][i].twistAOA;

            sens[0][i].phaseShift = cos(AOA_rad) * sens[0][i].phaseShift + sin(AOA_rad) * sens[0][i].phaseShift;
            sens[1][i].phaseShift = -sin(AOA_rad) * sens[0][i].phaseShift + cos(AOA_rad) * sens[1][i].phaseShift;

            sens[0][i].twistMag = cos(AOA_rad) * sens[0][i].twistMag + sin(AOA_rad) * sens[0][i].twistMag;
            sens[1][i].twistMag = -sin(AOA_rad) * sens[0][i].twistMag + cos(AOA_rad) * sens[1][i].twistMag;

            sens[0][i].chord = cos(AOA_rad) * sens[0][i].chord + sin(AOA_rad) * sens[0][i].chord;
            sens[1][i].chord = -sin(AOA_rad) * sens[0][i].chord + cos(AOA_rad) * sens[1][i].chord;

            sens[0][i].AOA = -sin(AOA_rad) * AOA_rad_sens * pos[0][i] + cos(AOA_rad) * AOA_rad_sens * pos[1][i];
            sens[1][i].AOA = -cos(AOA_rad) * AOA_rad_sens * pos[0][i] - sin(AOA_rad) * AOA_rad_sens * pos[1][i];

            //shift airfoil chord back so LE at origin
            pos[0][i] += airfoilDV.chord;
            sens[0][i].chord += 1.0;
        }
    }

    void printToVtk(std::string dvType) {
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

        string scalarName;
        double total_sens[NXI];
        if (dvType == "toverc") {
            scalarName = "TOVERC_SENS";
            for (int i = 0; i < NXI; i++) {
                total_sens[i] = pow( pow(sens[0][i].toverc,2) + pow(sens[1][i].toverc,2), 0.5);
                cout << "total_sens = " << total_sens[i] << std::endl;
            }
            myfile << "POINT_DATA " << NXI << "\n";
            myfile << "SCALARS " << scalarName << " double64 1\n";
            myfile << "LOOKUP_TABLE default\n";
            for (int i = 0; i < NXI; i++) {
                myfile << total_sens[i];
                myfile << "\n";
            }
        } else if (dvType == "twistAOA") {
            scalarName = "TWISTAOA_SENS";
            for (int i = 0; i < NXI; i++) {
                total_sens[i] = pow( pow(sens[0][i].twistAOA,2) + pow(sens[1][i].twistAOA,2), 0.5);
                cout << "total_sens = " << total_sens[i] << std::endl;
            }
            myfile << "POINT_DATA " << NXI << "\n";
            myfile << "SCALARS " << scalarName << " double64 1\n";
            myfile << "LOOKUP_TABLE default\n";
            for (int i = 0; i < NXI; i++) {
                myfile << total_sens[i];
                myfile << "\n";
            }
        } else if (dvType == "phaseShift") {
            scalarName = "PHASESHIFT_SENS";
            for (int i = 0; i < NXI; i++) {
                total_sens[i] = pow( pow(sens[0][i].phaseShift,2) + pow(sens[1][i].phaseShift,2), 0.5);
            }
            myfile << "POINT_DATA " << NXI << "\n";
            myfile << "SCALARS " << scalarName << " double64 1\n";
            myfile << "LOOKUP_TABLE default\n";
            for (int i = 0; i < NXI; i++) {
                myfile << total_sens[i];
                myfile << "\n";
            }
        } else if (dvType == "twistMag") {
            scalarName = "TWISTMAG_SENS";
            for (int i = 0; i < NXI; i++) {
                total_sens[i] = pow( pow(sens[0][i].twistMag,2) + pow(sens[1][i].twistMag,2), 0.5);
            }
            myfile << "POINT_DATA " << NXI << "\n";
            myfile << "SCALARS " << scalarName << " double64 1\n";
            myfile << "LOOKUP_TABLE default\n";
            for (int i = 0; i < NXI; i++) {
                myfile << total_sens[i];
                myfile << "\n";
            }
        } else if (dvType == "chord") {
            scalarName = "CHORD_SENS";
            for (int i = 0; i < NXI; i++) {
                total_sens[i] = pow( pow(sens[0][i].chord,2) + pow(sens[1][i].chord,2), 0.5);
            }
            myfile << "POINT_DATA " << NXI << "\n";
            myfile << "SCALARS " << scalarName << " double64 1\n";
            myfile << "LOOKUP_TABLE default\n";
            for (int i = 0; i < NXI; i++) {
                myfile << total_sens[i];
                myfile << "\n";
            }
        } else if (dvType == "AOA") {
            scalarName = "AOA_SENS";
            for (int i = 0; i < NXI; i++) {
                total_sens[i] = pow( pow(sens[0][i].AOA,2) + pow(sens[1][i].AOA,2), 0.5);
            }
            myfile << "POINT_DATA " << NXI << "\n";
            myfile << "SCALARS " << scalarName << " double64 1\n";
            myfile << "LOOKUP_TABLE default\n";
            for (int i = 0; i < NXI; i++) {
                myfile << total_sens[i];
                myfile << "\n";
            }
        }

        

        myfile.close();
    }
};

#endif //AIRFOIL_H