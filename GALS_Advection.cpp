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
//#include "Bilinear.h"
#include "Hermite.h"
#include "Initializer_vortex.h"    //UPDATE THIS HEADER FILE WITH THE RELEVANT FUNCTIONS

using namespace std;

// y direction is the first index of array, x direction is the second index of array

int main(){

    clock_t tStart = clock();

    struct  timeval start1;       struct  timeval end1;          unsigned  long diff1;
    struct  timeval start2;       struct  timeval end2;          unsigned  long diff2;
    struct  timeval start3;       struct  timeval end3;          unsigned  long diff3;
    struct  timeval start4;       struct  timeval end4;          unsigned  long diff4;
    struct  timeval start5;       struct  timeval end5;          unsigned  long diff5;
    struct  timeval start999;     struct  timeval end999;        unsigned  long diff999;

    /* UPDATE ALL THE FOLLOWING VALUES */
    double xlim1 = 0.0;                       //Lower limit on x-axis
    double xlim2 = 1.0;                      //Upper limit on x-axis
    int nx = 192 + 1;                         //Number of nodes in x-direction INCLUDING THE EXTREME VALUES

    double ylim1 = 0.0;                       //Lower limit on y-axis
    double ylim2 = 1.0;                     //Upper limit on y-axis
    int ny = 192 + 1;                        //Number of nodes INCLUDING THE EXTREME VALUES

    double dt = (1/100.0);                     //Length of time step
    double Tfinal = 0.01;                    //Total time period for the simulation
    int option = 2;                         //Option - if you need animation initialize at 1 else initialize at 2
    int printstep = 8;                      //How frequently do you want to store the images (every nth time step)
    char psischeme[] = "SuperConsistent";   //'SuperConsistent' or 'Heuns'
    char backtrace_scheme[] = "RK3" ;      //'Euler' or 'RK3'
    double T_period = 1.0;                  //Period of the velocity field

    //MAKE SURE THAT YOU HAVE ENOUGH MEMORY SPACE IF YOU ARE STORING A LOT OF TIME STEP VALUES BECAUSE IT STORES ACROSS GRID POINTS FOR EACH PRINTSTEP

    /* USER UPDATE OVER */

    int n = Tfinal/dt; //Number of time steps
    // If only the initial and final profiles are needed
    if(option != 1)
        printstep = n;


    /*
     Let the following represent one cell

     2      3
     *      *


     *      *
     0      1

     the value of the loop for this cell varies from 0 -> 3 but the (i, j) coordinates that represent these points in an array change as (i,j), (i,j+1), (i+1,j), (i+1, j+1)

     Hence, tempindexes take care of these changes

     */


    // --------------- define all required parameters ---------------------------------------------------------
    #include "parameter_setting.h"
    //---------------------------------------------------------------------------------------------------------


    // --------------- 1ST Part --------------------------------------------------------------------------------
    // Initializing at t = 0
        gettimeofday(&start1,NULL);

        #include "initialization.h"

        gettimeofday(&end1,NULL);
        diff1 = 1000000 * (end1.tv_sec-start1.tv_sec)+ end1.tv_usec-start1.tv_usec;
    //---------------------------------------------------------------------------------------------------------


    // -------------- 2ND Part --------------------------------------------------------------------------------
    // Giving an extra empty line for distinction between next time step
    myfile << '\n';           myfile1 << '\n';             myfile2 << '\n';            myfile3 << '\n';             uu << '\n';         vv << '\n';

    details<< nx << "," << ny << "," << std::fixed << std::setprecision(10) << dx << "," << dy << "," << xlim1 << "," << xlim2 << "," << ylim1 << "," << ylim2 << "," << n << "," << dt << "," << printstep;
    details.close();

    gettimeofday(&start2,NULL);

    // TIME STEPPING LOOP
    for(t = 0; t < n; t++){

        // Initializing tracker to zero at every time step
        for(int p = 0; p < ny; p ++){
            for(int q = 0; q<nx; q++){
                tracker[p][q] = 0;
            }
        }

        gettimeofday(&start999,NULL);

        for(int i = 0; i < ny-1; i++){  // loop for y - rows in the array
            for(int j = 0; j < nx-1; j++){  // loop for x - columns in the array

                // Storing the four values for four nodes of each cell
                phi[0] = mphi[i][j];            psix[0] = mpsix[i][j];          psiy[0] = mpsiy[i][j];          psixy[0] = mpsixy[i][j];
                phi[1] = mphi[i][j+1];          psix[1] = mpsix[i][j+1];        psiy[1] = mpsiy[i][j+1];        psixy[1] = mpsixy[i][j+1];
                phi[2] = mphi[i+1][j];          psix[2] = mpsix[i+1][j];        psiy[2] = mpsiy[i+1][j];        psixy[2] = mpsixy[i+1][j];
                phi[3] = mphi[i+1][j+1];        psix[3] = mpsix[i+1][j+1];      psiy[3] = mpsiy[i+1][j+1];      psixy[3] = mpsixy[i+1][j+1];
                // Node value assignment ends

                // Storing the coordinates of first node of the working cell
                xo = x[j];                yo = y[i];
                corner = -1;

                // calculate 4 nodal values on the cell
                #include "update_node_data.h"

            }   // end of x-columns marching
        }   // end of y-rows marching

        gettimeofday(&end999,NULL);
        diff999 = 1000000 * (end999.tv_sec - start999.tv_sec) + end999.tv_usec - start999.tv_usec;
    //---------------------------------------------------------------------------------------------------------


    //---------------------------------------------------------------------------------------------------------
    // Update the mixed derivatives now for the remaining grid points
        gettimeofday(&start3,NULL);

        #include "update_mixed_derivatives.h"

        gettimeofday(&end3,NULL);
        diff3 = 1000000 * (end3.tv_sec - start3.tv_sec) + end3.tv_usec - start3.tv_usec;
    //---------------------------------------------------------------------------------------------------------


    //---------------------------------------------------------------------------------------------------------
    // Feeding values back to the master matrix
        gettimeofday(&start4,NULL);

        for(int k = 0; k < ny; k++){
            for(int l = 0; l < nx; l++){
                mphi[k][l] = tempphi[k][l];
                mpsix[k][l] = temppsix[k][l];
                mpsiy[k][l] = temppsiy[k][l];
                mpsixy[k][l] = temppsixy[k][l];
            }
        }

        gettimeofday(&end4,NULL);
        diff4 = 1000000 * (end4.tv_sec - start4.tv_sec) + end4.tv_usec-start4.tv_usec;
    //---------------------------------------------------------------------------------------------------------


    //---------------------------------------------------------------------------------------------------------
    // Feeding phi, psix, psiy and psixy values in their respective files
        gettimeofday(&start5,NULL);

            #include "feeding_output.h"

            cout<< t+1;
            cout<< " Time Step Completed" <<'\n';

        gettimeofday(&end5,NULL);
        diff5 = 1000000 * (end5.tv_sec - start5.tv_sec) + end5.tv_usec - start5.tv_usec;

    //---------------------------------------------------------------------------------------------------------

    }  // end of time marching loop

    gettimeofday(&end2,NULL);
    diff2 = 1000000 * (end2.tv_sec - start2.tv_sec) + end2.tv_usec - start2.tv_usec;

    myfile.close();       myfile1.close();         myfile2.close();          myfile3.close();         uu.close();        vv.close();

    cout<<diff1*0.001<<endl;
    cout<<diff2*0.001<<endl;
    cout<<diff3*0.001<<endl;
    cout<<diff4*0.001<<endl;
    cout<<diff5*0.001<<endl;
    cout<<diff999*0.001<<endl;
//    double ttt = (clock() - tStart)/CLOCKS_PER_SEC;
//    cout<<'\n'<<"Time Taken: "<<ttt<<'\n'<< '\n';
    return 0;
}
