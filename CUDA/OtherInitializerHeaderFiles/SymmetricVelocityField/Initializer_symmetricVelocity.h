//
//  Initializer.h
//
//
//  Created by Raunak Bardia on 10/10/14.
//
//

#ifndef _Initializer_h
#define _Initializer_h

#include <math.h>
#include <stdio.h>
#include <string.h>

const double pi = M_PI;
const double xo = .5;
const double yo = .5;
const double rcircle = .15;
const double xc = .5;
const double yc = .75;

using namespace std;

double Velx(double x, double y, double t, double T_period){
    double temp;
    temp = 0.5 * pi * cos(pi * t/T_period) * ((2*sin(2*pi*xc)*sin(2*pi*yc))*(x-xc)+2*(-cos(2*pi*xc)*(sin(pi*yc)*sin(pi*yc))+(sin(pi*xc)*sin(pi*xc))*cos(2*pi*yc))*(y-yc) );
    return temp;
}

double Vely(double x, double y,double t, double T_period){
    double temp;
    temp = 0.5 * pi * cos(pi * t/T_period) * ((-2*sin(2*pi*xc)*sin(2*pi*yc))*(y-yc)+2*(-cos(2*pi*xc)*(sin(pi*yc)*sin(pi*yc)) +(sin(pi*xc)*sin(pi*xc))*cos(2*pi*yc))*(x-xc) );
    return temp;
}

double gradUx(double x, double y, double t, double T_period){
    double temp;
    temp = 0.5*pi*cos(pi*t/T_period)*( ( 2*sin(2*pi*xc)*sin(2*pi*yc)));
    return temp;
}

double gradUy(double x, double y, double t, double T_period){
    double temp;
	temp = 0.5*pi*cos(pi*t/T_period)*(2*(-cos(2*pi*xc)*pow(sin(pi*yc),2)+pow(sin(pi*xc),2)*cos(2*pi*yc)) );
    return temp;
}

double gradVx(double x, double y, double t, double T_period){
    double temp;
    temp = gradUy(x,y,t,T_period);
    return temp;
}

double gradVy(double x, double y, double t, double T_period){
    double temp;
	temp = -gradUx(x,y,t,T_period);
    return temp;
}
double initialize(double x, double y){
    double k;
    k = exp(-(pow((x - xo),2) + pow((y - yo),2))) - exp(-pow(rcircle,2));
    return k;
}

double derivxinit(double x, double y){
    double k;
    k = -2*(x-xo)*exp(-(pow((x - xo),2) + pow((y - yo),2)));
    return k;
}

double derivyinit(double x, double y){
    double k;
    k = -2*(y-yo)*exp(-(pow((x - xo),2) + pow((y - yo),2)));
    return k;
}

double derivxyinit(double x, double y){
    double k;
    k = 4 * (y - yo) * (x - xo) * exp(-(pow((x - xo),2) + pow((y - yo),2)));
    return k;
}

#endif
