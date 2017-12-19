//
//  OscillatingVelocityCUDA.cu
//
//
//  Created by Raunak Bardia on 12/18/17.
//
//

#ifndef _OscillatingVelocityCUDA_cu
#define _OscillatingVelocityCUDA_cu

__device__ __host__ double Velx(double x, double y, double t, double T_period){
    double temp;
    double signum = cos(pi * t/T_period);
    if(signum == 0)
        signum = 0;
    else if(signum < 0)
        signum = -1;
    else
        signum = 1;
    temp = (50.0 - x) * signum;
    return temp;
}

__device__ __host__ double Vely(double x, double y,double t, double T_period){
    double temp;
    double signum = cos(pi * t/T_period);
    if(signum == 0)
        signum = 0;
    else if(signum < 0)
        signum = -1;
    else
        signum = 1;
    temp = (y - 50.0) * signum;
    return temp;
}

__device__ double gradUx(double x, double y, double t, double T_period){
    double temp;
    double signum = cos(pi * t/T_period);
    if(signum == 0)
        signum = 0;
    else if(signum < 0)
        signum = -1;
    else
        signum = 1;
    temp = -1 * signum;
    return temp;
}

__device__ double gradUy(double x, double y, double t, double T_period){
    double temp;
    double signum = cos(pi * t/T_period);
    if(signum == 0)
        signum = 0;
    else if(signum < 0)
        signum = -1;
    else
        signum = 1;
    temp = 0;
    return temp;
}

__device__ double gradVx(double x, double y, double t, double T_period){
    double temp;
    double signum = cos(pi * t/T_period);
    if(signum == 0)
        signum = 0;
    else if(signum < 0)
        signum = -1;
    else
        signum = 1;
    temp = 0;
    return temp;
}

__device__ double gradVy(double x, double y, double t, double T_period){
    double temp;
    double signum = cos(pi * t/T_period);
    if(signum == 0)
        signum = 0;
    else if(signum < 0)
        signum = -1;
    else
        signum = 1;
    temp = signum;
    return temp;
}
#endif
