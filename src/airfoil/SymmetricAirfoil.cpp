
#define _USE_MATH_DEFINES

#include <array>
#include <cmath>
#include <iostream>
#include "SymmetricAirfoil.h"
using namespace std;

template<unsigned int N>
SymmetricAirfoil<N>::SymmetricAirfoil(double tc) {
    // store the value toverc
    toverc = tc;

    // build the array of points
    for (size_t ipoint{0}; ipoint < N; ipoint++) {
        
        // set xi the angle of each point
        double xi = ipoint/(N-1) * M_PI;

        // set the x value of each point
        points[ipoint][0] = 0.5*(cos(xi) - 1.0);
        points[ipoint][1] = 2.0/3 * toverc * sin(xi) * sin(xi/2);
    }
}

template<unsigned int N>
void SymmetricAirfoil<N>::printPoints() {
    // loop over each point on the airfoil and print it
    for (auto const& point : points) {
        for (auto const& coord : point) {
            cout << coord << " ";
        }
        cout << endl;
    }
}
