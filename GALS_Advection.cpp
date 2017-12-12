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
#include <tuple>
#include "VariableDefinitions.h"
#include "Hermite.h"
#include "InitializeLevelSet.h"    //UPDATE THIS HEADER FILE WITH THE RELEVANT FUNCTIONS
#include "VortexVelocity.h"    //UPDATE THIS HEADER FILE WITH THE RELEVANT FUNCTIONS
#include "TimeSteppingMethods.h"
#include "Allocation.h"
#include "AdvectionPointCalcs.h"
#include "Constants.h"

using namespace std;
using namespace galsfunctions;

// y direction is the first index of array, x direction is the second index of array

int main(){

    /* UPDATE ALL THE FOLLOWING VALUES */
    double xlim1 = 0.0;                       //Lower limit on x-axis
    double xlim2 = 1.0;                      //Upper limit on x-axis
    unsigned int nx = 127 + 1;                         //Number of nodes in x-direction INCLUDING THE EXTREME VALUES

    double ylim1 = 0.0;                       //Lower limit on y-axis
    double ylim2 = 1.0;                     //Upper limit on y-axis
    unsigned int ny = 127 + 1;                        //Number of nodes INCLUDING THE EXTREME VALUES

    double dt = (1.0/128.0);                     //Length of time step
    double Tfinal = 1.0;                    //Total time period for the simulation
    unsigned int option = 1;                         //Option - if you need animation initialize at 1 else initialize at 2
    unsigned int printstep = 8;                      //How frequently do you want to store the images (every nth time step)
    char psischeme[] = "SuperConsistent";   //'SuperConsistent' or 'Heuns'
    char backtrace_scheme[] = "RK3" ;      //'Euler' or 'RK3'
    double T_period = 1.0;                  //Period of the velocity field

    //MAKE SURE THAT YOU HAVE ENOUGH MEMORY SPACE IF YOU ARE STORING A LOT OF TIME STEP VALUES BECAUSE IT STORES ACROSS GRID POINTS FOR EACH PRINTSTEP

    /* USER UPDATE OVER */
    unsigned int n = Tfinal/dt; //Number of time steps
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
    remove("phi.txt");
    remove("psix.txt");
    remove("psiy.txt");
    remove("psixy.txt");
    remove("details.txt");
    remove("Velocity_x.txt");
    remove("Velocity_y.txt");
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
    for(unsigned int t = 0; t < n; t++){

		intgridarray tracker;
    		tracker.resize(ny,intvectorarray(nx,0));
    
    		gridarray xadv, yadv;
    		xadv.resize(ny,vectorarray(nx,0.0));
    		yadv.resize(ny,vectorarray(nx,0.0));
        
		advection_point(x, y, xadv, yadv, t, dt, T_period, backtrace_scheme);
        
    		intgridarray cellx, celly;
    		cellx.resize(ny,intvectorarray(nx,0));
    		celly.resize(ny,intvectorarray(nx,0));
	
		find_advection_point_location(x, y, xadv, yadv, cellx, celly, tracker, xlim1, xlim2, ylim1, ylim2);

		update_levelset_data(x, y, xadv, yadv, cellx, celly, tracker, t, dt, 
				tempphi, temppsix, temppsiy, temppsixy, mphi, mpsix, mpsiy, mpsixy,
				psischeme,backtrace_scheme,T_period);	
        
        // Update the mixed derivatives now for the remaining grid points
        	update_mixed_derivatives(temppsix, temppsiy, temppsixy, nx, ny, dx, dy);

        //---------------------------------------------------------------------------------------------------------
        // Feeding values back to the master matrix

        mphi = tempphi;
        mpsix = temppsix;
        mpsiy = temppsiy;
        mpsixy = temppsixy;

        //---------------------------------------------------------------------------------------------------------
        // Feeding phi, psix, psiy and psixy values in their respective files
        if((t+1) % printstep == 0)
            fileprint(mphi,mpsix,mpsiy,mpsixy,nx,ny,x,y,(t+1)*dt,T_period);

        cout<< t+1;
        cout<< " Time Step Completed" <<'\n';

        //---------------------------------------------------------------------------------------------------------
        xadv.clear();
        yadv.clear();
        tracker.clear();
        cellx.clear();
        celly.clear();
    }  // end of time marching loop

    return 0;
}
