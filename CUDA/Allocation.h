//
//  Allocation.cpp
//
//
//  Created by Raunak Bardia on 12/01/17.
//
//
#include "InitializeLevelSet.h"
#include "VortexVelocityCUDA.cu"

using namespace std;

    void allocate_levelset_matrices(double *mphi, double *mpsix, double *mpsiy, double *mpsixy, double *x, double *y, unsigned int nx, unsigned int ny)
    {
        // master matrices of phi, psix, psiy, psixy which are updated at each time step & temporary matrices which is updated constantly throughout the cell iterations
        for(unsigned int i = 0; i < ny; i++){
            for(unsigned int j = 0; j < nx; j++){
    		unsigned long indexToWrite = i * nx + j;
                mphi[indexToWrite] = initialize(x[j], y[i]);
                mpsix[indexToWrite] = derivxinit(x[j], y[i]);
                mpsiy[indexToWrite] = derivyinit(x[j], y[i]);
                mpsixy[indexToWrite] = derivxyinit(x[j], y[i]);
            }
        }
    }
    
    void allocate_levelset_velocity(double *u, double *v, unsigned int nx, unsigned int ny, double *x, double *y, double time, double T_period)
    {
        for(unsigned int i = 0; i < ny; i++){
            for(unsigned int j = 0; j < nx; j++){
    		unsigned long indexToWrite = i * nx + j;
                u[indexToWrite] = Velx(x[j], y[i], time, T_period);
                v[indexToWrite] = Vely(x[j], y[i], time, T_period);
            }
        }
    }
    
    void gridnodes(double *x, double *y, double xlim1, double ylim1, double dx, double dy, unsigned int nx, unsigned int ny)
    {
        // Defining node points
        for(unsigned int i = 0; i < nx; i++)
            x[i] = xlim1 + dx * i;
        for(unsigned int i = 0; i < ny; i++)
            y[i] = ylim1 + dy * i;
        // Node point definition ends
    }
    
    void fileprint(double *mphi, double *mpsix, double *mpsiy, double *mpsixy, unsigned int nx, unsigned int ny,
            double *x, double *y, double time, double T_period)
    {
        unsigned long gridmemory = nx * ny * sizeof(double);
        double* u = (double*) malloc(gridmemory);
        double* v = (double*) malloc(gridmemory);
        allocate_levelset_velocity(u,v,nx,ny,x,y,time,T_period);
        
        // Opening new files
        ofstream phifile, psixfile, psiyfile, psixyfile, ufile, vfile;
        phifile.open("phi.txt", ios::out | ios::app);
        psixfile.open("psix.txt", ios::out | ios::app);
        psiyfile.open("psiy.txt", ios::out | ios::app);
        psixyfile.open("psixy.txt", ios::out | ios::app);
        ufile.open("Velocity_x.txt", ios::out | ios::app);
        vfile.open("Velocity_y.txt",ios::out | ios::app);
        
        for(unsigned int i = 0; i < ny; i++){
            
            // Initializing the first element of each row separately to avoid trailing commas in the file
            phifile << std::fixed << std::setprecision(10) << mphi[i * nx];
            psixfile << std::fixed << std::setprecision(10) << mpsix[i * nx];
            psiyfile << std::fixed << std::setprecision(10) << mpsiy[i * nx];
            psixyfile << std::fixed << std::setprecision(10) << mpsixy[i * nx];
            ufile << std::fixed << std::setprecision(10) << u[i * nx];
            vfile << std::fixed << std::setprecision(10) << v[i * nx];
            
            for(unsigned int j = 1; j < nx; j++){
    		unsigned long indexToWrite = i * nx + j;
                phifile << "," << std::fixed << std::setprecision(10) << mphi[indexToWrite];
                psixfile << "," << std::fixed << std::setprecision(10) << mpsix[indexToWrite];
                psiyfile << "," << std::fixed << std::setprecision(10) << mpsiy[indexToWrite];
                psixyfile << "," << std::fixed << std::setprecision(10) << mpsixy[indexToWrite];
                ufile << "," << std::fixed << std::setprecision(10) << u[indexToWrite];
                vfile << "," << std::fixed << std::setprecision(10) << v[indexToWrite];
            }
            phifile << '\n';
            psixfile << '\n';
            psiyfile << '\n';
            psixyfile << '\n';
            ufile << '\n';
            vfile << '\n';
        }
        
        phifile << '\n';           psixfile << '\n';             psiyfile << '\n';            psixyfile << '\n';             ufile << '\n';         vfile << '\n';
        phifile.close();       psixfile.close();         psiyfile.close();          psixyfile.close();         ufile.close();        vfile.close();
    }
