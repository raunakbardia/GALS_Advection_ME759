//
//  GALS - Before compiling the program, update the section of this program that is in the beginning of main and update all the functions in Initializer.h
//
//  Created by Raunak Bardia on 10/22/14.
//
// DISCLAIMER:
// Use the indexes carefully
// First index of array represents movement along y-direction because it represents rows
// Second index of array represents movement along x-direction because it represents columns
//
// Implementing GALS for a given initial level set function
// in a specified velocity field for a grid of cells
//
// Given -
// 1. Defining function at t = 0 which implies that phi and psi values are available for all node points at t = 0
// 2. Given velocity for the complete domain at all times
//
// All required data is stored in separate 2D matrices of phi, psix, psiy and psixy
// Boundary Condition grad(velocity).n > 0
//

// THIS IMPLEMENTATION WON'T WORK IF THE GRID IS SMALLER THAN (2 X 2)
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <vector>
#include "VariableDefinitions.h"
#include "Hermite.h"
#include "Initializer_vortex.h"    //UPDATE THIS HEADER FILE WITH THE RELEVANT FUNCTIONS
#include "TimeSteppingMethods.h"
#include "Allocation.h"
//#include "Bilinear.h"

using namespace std;
using namespace galsfunctions;

// y direction is the first index of array, x direction is the second index of array

int main(){
    
    /* UPDATE ALL THE FOLLOWING VALUES */
    double xlim1 = 0.0;                       //Lower limit on x-axis
    double xlim2 = 1.0;                      //Upper limit on x-axis
    int nx = 64 + 1;                         //Number of nodes in x-direction INCLUDING THE EXTREME VALUES
    
    double ylim1 = 0.0;                       //Lower limit on y-axis
    double ylim2 = 1.0;                     //Upper limit on y-axis
    int ny = 64 + 1;                        //Number of nodes INCLUDING THE EXTREME VALUES
    
    double dt = (1/100.0);                     //Length of time step
    double Tfinal = 2.0;                    //Total time period for the simulation
    int option = 2;                         //Option - if you need animation initialize at 1 else initialize at 2
    int printstep = 8;                      //How frequently do you want to store the images (every nth time step)
    char psischeme[] = "SuperConsistent";   //'SuperConsistent' or 'Heuns'
    char backtrace_scheme[] = "RK3" ;      //'Euler' or 'RK3'
    double T_period = 1.0;                  //Period of the velocity field
    
    //MAKE SURE THAT YOU HAVE ENOUGH MEMORY SPACE IF YOU ARE STORING A LOT OF TIME STEP VALUES BECAUSE IT STORES ACROSS GRID POINTS FOR EACH PRINTSTEP
    
    /* USER UPDATE OVER */
    double np = nx * ny; // Number of node points
    int cornerx, cornery;   //To store indexes of corners - Not sure where they are being used
    int n = Tfinal/dt; //Number of time steps
    if(option != 1)
        printstep = n;
    
    // Node Locations
    double dx = (xlim2 - xlim1)/(nx - 1);
    double dy = (ylim2 - ylim1)/(ny - 1);
    vectorarray x,y;
    gridnodes(x,y,xlim1,ylim1,dx,dy,nx,ny);
    
    // Initializing at t = 0
    gridarray mphi, mpsix, mpsiy, mpsixy; // level set matrix
    allocate_levelset_matrices(mphi,mpsix,mpsiy,mpsixy,x,y,nx,ny); //Initializing level set matrices
    
    gridarray tempphi(mphi), temppsix(mpsix), temppsiy(mpsiy), temppsixy(mpsixy);	// Making level set matrix copies for time loop use
    
    
// Removing existing files with these names if any
    int i = remove("phi.txt");
    i = remove("psix.txt");
    i = remove("psiy.txt");
    i = remove("psixy.txt");
    i = remove("details.txt");
    i = remove("Velocity_x.txt");
    i = remove("Velocity_y.txt");
    fileprint(mphi,mpsix,mpsiy,mpsixy,nx,ny,x,y,0.0,T_period);
    
    /*
     * Let the following represent one cell
     *
     * 2      3
     *      *
     *
     *
     *      *
     * 0      1
     *
     * the value of the loop for this cell varies from 0 -> 3 but the (i, j) coordinates that represent these points in an array change as (i,j), (i,j+1), (i+1,j), (i+1, j+1)
     *
     * Hence, tempindexes take care of these changes
     *
     */
    
    ofstream details;
    details.open("details.txt", ios::out | ios::app);
    details<< nx << "," << ny << "," << std::fixed << std::setprecision(10) << dx << "," << dy << "," << xlim1 << "," << xlim2 << "," << ylim1 << "," << ylim2 << "," << n << "," << dt << "," << printstep;
    details.close();
    
    // TIME STEPPING LOOP
    // If only the initial and final profiles are needed
    for(int t = 0; t < n; t++){
        
        gridarray tracker;
        tracker.resize(ny,vectorarray(nx,0.0));
        
        gridarray xadv, yadv;
        advection_point(x, y, xadv, yadv, t, dt, T_period, backtrace_scheme);
        
        for(int i = 0; i < ny-1; i++){  // loop for y - rows in the array
            for(int j = 0; j < nx-1; j++){  // loop for x - columns in the array
                
                double phi[4], psix[4], psiy[4], psixy[4];
                // Storing the four values for four nodes of each cell
                phi[0] = mphi[i][j];            psix[0] = mpsix[i][j];          psiy[0] = mpsiy[i][j];          psixy[0] = mpsixy[i][j];
                phi[1] = mphi[i][j+1];          psix[1] = mpsix[i][j+1];        psiy[1] = mpsiy[i][j+1];        psixy[1] = mpsixy[i][j+1];
                phi[2] = mphi[i+1][j];          psix[2] = mpsix[i+1][j];        psiy[2] = mpsiy[i+1][j];        psixy[2] = mpsixy[i+1][j];
                phi[3] = mphi[i+1][j+1];        psix[3] = mpsix[i+1][j+1];      psiy[3] = mpsiy[i+1][j+1];      psixy[3] = mpsixy[i+1][j+1];
                // Node value assignment ends
                
                // Storing the coordinates of first node of the working cell
                double xo = x[j], yo = y[i];
                int corner = -1;    // To save which node of the cell is the grid corner if a corner is encountered
                
                // calculate 4 nodal values on the cell
#include "update_node_data.h"
                
            }   // end of x-columns marching
        }   // end of y-rows marching
        
        //---------------------------------------------------------------------------------------------------------
        
        
        //---------------------------------------------------------------------------------------------------------
        // Update the mixed derivatives now for the remaining grid points
#include "update_mixed_derivatives.h"
        
        //---------------------------------------------------------------------------------------------------------
        
        
        //---------------------------------------------------------------------------------------------------------
        // Feeding values back to the master matrix
        
        mphi = tempphi;
        mpsix = temppsix;
        mpsiy = temppsiy;
        mpsixy = temppsixy;
        
        //---------------------------------------------------------------------------------------------------------
        
        
        //---------------------------------------------------------------------------------------------------------
        // Feeding phi, psix, psiy and psixy values in their respective files
        if((t+1) % printstep == 0)
            fileprint(mphi,mpsix,mpsiy,mpsixy,nx,ny,x,y,(t+1)*dt,T_period);
        
        cout<< t+1;
        cout<< " Time Step Completed" <<'\n';
        
        //---------------------------------------------------------------------------------------------------------
        
    }  // end of time marching loop
    
    return 0;
}
