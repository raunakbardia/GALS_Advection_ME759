//
//  SymmetricVelocityCUDA.cu
//
//
//  Created by Raunak Bardia on 12/18/17.
//
//

#ifndef _SymmetricVelocityCUDA_cu
#define _SymmetricVelocityCUDA_cu

__device__ __host__ double Velx(double x, double y, double t, double T_period){
    double temp;
    temp = 0.5 * pi * cos(pi * t/T_period) * ((2*sin(2*pi*xc)*sin(2*pi*yc))*(x-xc)+2*(-cos(2*pi*xc)*(sin(pi*yc)*sin(pi*yc))+(sin(pi*xc)*sin(pi*xc))*cos(2*pi*yc))*(y-yc) );
    return temp;
}

__device__ __host__ double Vely(double x, double y,double t, double T_period){
    double temp;
    temp = 0.5 * pi * cos(pi * t/T_period) * ((-2*sin(2*pi*xc)*sin(2*pi*yc))*(y-yc)+2*(-cos(2*pi*xc)*(sin(pi*yc)*sin(pi*yc)) +(sin(pi*xc)*sin(pi*xc))*cos(2*pi*yc))*(x-xc) );
    return temp;
}

__device__ double gradUx(double x, double y, double t, double T_period){
    double temp;
    temp = 0.5*pi*cos(pi*t/T_period)*( ( 2*sin(2*pi*xc)*sin(2*pi*yc)));
    return temp;
}

__device__ double gradUy(double x, double y, double t, double T_period){
    double temp;
    temp = 0.5*pi*cos(pi*t/T_period)*(2*(-cos(2*pi*xc)*pow(sin(pi*yc),2)+pow(sin(pi*xc),2)*cos(2*pi*yc)) );
    return temp;
}

__device__ double gradVx(double x, double y, double t, double T_period){
    double temp;
    temp = gradUy(x,y,t,T_period);
    return temp;
}

__device__ double gradVy(double x, double y, double t, double T_period){
    double temp;
    temp = -gradUx(x,y,t,T_period);
    return temp;
}
#endif
