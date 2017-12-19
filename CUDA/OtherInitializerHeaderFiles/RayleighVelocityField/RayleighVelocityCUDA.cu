//
//  RayleighVelocityCUDA.cu
//
//
//  Created by Raunak Bardia on 12/18/17.
//
//

#ifndef _RayleighVelocityCUDA_cu
#define _RayleighVelocityCUDA_cu

__device__ __host__ double Velx(double x, double y, double t, double T_period){
    double temp;
    temp = -(1 - exp(c * pow(x,na))) * (U0 * cos(k * y) + U1 * cos(2 * k * y)) * exp(omega * t);
    return temp;
}

__device__ __host__ double Vely(double x, double y,double t, double T_period){
    double temp;
    temp = -(c * na * pow(x, na - 1))/(2 * k) * exp(c * pow(x,na) + omega * t) * (2 * U0 * sin(k * y) + U1 * sin(2 * k * y));
    return temp;
}

__device__ double gradUx(double x, double y, double t, double T_period){
    double temp;
    temp = c * na * pow(x, na - 1) * exp(c * pow(x,na) + omega * t) * (U0 * cos(k * y) + U1 * cos(2 * k * y));
    return temp;
}

__device__ double gradUy(double x, double y, double t, double T_period){
    double temp;
    temp = -k * exp(omega * t) * (exp(c * pow(x,na)) - 1) * sin(k * y) * (U0 + 4 * U1 * cos(k * y));
    return temp;
}

__device__ double gradVx(double x, double y, double t, double T_period){
    double temp;
    temp = -c/k * na * pow(x, na - 2) * (c * na * pow(x,na) + na - 1) * sin(k * y) * exp(c * pow(x, na) + omega * t) * (U0 + U1 * cos(k * y));
    return temp;
}

__device__ double gradVy(double x, double y, double t, double T_period){
    double temp;
    temp = -c * na * pow(x, na - 1) * exp(c * pow(x,na) + omega * t) * (U0 * cos(k * y) + U1 * cos(2 * k * y));
    return temp;
}
#endif
