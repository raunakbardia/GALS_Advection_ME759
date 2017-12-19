//
//  InitializeLevelSet.h
//
//
//  Created by Raunak Bardia on 12/18/17.
//
//

const double k = 0.25;              // Wavenumber
const double omega = 5.0;            // Rate of exponential growth for modes -- arbitrarily chosen here
const double U0 = 1.0;               // Velocity field
const double c = 0.01;                // Rate at which velocity vanishes near r=0   ---  ~(1-exp(c*r)) -->0 as r-->0
const double U1 = U0;
const double na = 0.1;
const double A = 0.0;
const double H = -1.0;

using namespace std;
double initialize(double x, double y){
    double k;
    k = A * y - x - H;
    return k;
}

double derivxinit(double x, double y){
    double k;
    k = -1;
    return k;
}

double derivyinit(double x, double y){
    double k;
    k = A;
    return k;
}

double derivxyinit(double x, double y){
    double k;
    k = 0;
    return k;
}
