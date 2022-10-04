
/*
Georgia Tech SMDO Lab, October 2022
Sean Engelstad, Aerospace Engineering

Class to generate symmetric airfoils with chord length 1
with leading edge at (0,0) and trailing edge at (1,0).
There is one design parameter here - t/c thickness to chord ratio.
*/

#ifndef SYM_AIRFOIL_H
#define SYM_AIRFOIL_H

#include <string>
#include <array>

template <unsigned int N> 
class SymmetricAirfoil {
public:
    // constructor builds airfoil curve
    SymmetricAirfoil(double);

    // print out the points to terminal
    void printPoints();
    
private:
    // toverc of the symmetric airfoil
    double toverc{0.0};

    // array of points on the curve
    std::array<std::array<double, 3>, N> points;
};

#endif //SYM_AIRFOIL_H