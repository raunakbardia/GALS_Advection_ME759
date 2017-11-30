// define all required parameters

#ifndef _PARAMETER_SETTING_H
#define _PARAMETER_SETTING_H

double xo, yo, rootpsix, rootpsiy; // {xo,yo} - first node in each cell; rootpsi* - temporary variables that stores updated psi*
    
    // Temporary variables to store the values of phi, psix, psiy and psixy at all the nodes of the cell being considered for the iteration
    double phi[4];
    double psix[4];
    double psiy[4];
    double psixy[4];
    int tempindex_x, tempindex_y;   // These variable select each node one by one in a cell
    int cornerx, cornery;   //To store indexes of corners
    int corner; // To save which node of the cell is the grid corner if a corner is encountered
    // Float constants for RK-3
    
    double c1 = (1/6.0);
    double c2 = (1/6.0);
    double c3 = (2/3.0);    //RK-3 constants
    
    // Defining node points
    double dx = (xlim2 - xlim1)/(nx - 1);
    double x[nx];
    for(int i = 0; i < nx; i++){
        x[i] = xlim1 + dx * i;
    }
    double dy = (ylim2 - ylim1)/(ny - 1);
    double y[ny];
    for(int i = 0; i < ny; i++){
        y[i] = ylim1 + dy * i;
    }
    // Node point definition ends
    
    int temp1, temp2;       // These temporary variables are updated when the advection point lies between the x and y limits of the cell
    int tracker[ny][nx];    // Tracker, tracks every node point and updates that particular index to 1 if that point has been updated
    
    // master matrices of phi, psix, psiy, psixy which are updated at each time step & temporary matrices which is updated constantly throughout the cell iterations
    double mphi[ny][nx], mpsix[ny][nx], mpsiy[ny][nx], mpsixy[ny][nx], tempphi[ny][nx], temppsix[ny][nx], temppsiy[ny][nx], temppsixy[ny][nx];
    
    double u[ny][nx], v[ny][nx];    //Storing the velocity at node points only for those time steps when printstep is initiated
    double d1, d2, d3, d4;  //Temporary variable to compute mixed derivatives
    double xadv, yadv;  //Advection point
    double np = nx * ny; // Number of node points
    
    // ux, vy - X velocity & Y velocity; ux1, vy1, ux2, vy2 - X & Y velocity of the xadv & yadv at separate RK steps
    double ux, vy, ux1, vy1, ux2, vy2;
    
    //FOR SUPERCONSISTENT METHOD & HEUNS
    double xr1,yr1,xr2,yr2,gradux1,graduy1,gradvx1,gradvy1,gradux2,graduy2,gradvx2,gradvy2,graduxnode,graduynode,gradvxnode,gradvynode,graduxroot,graduyroot,gradvxroot,gradvyroot;
    double gradxr1x,gradyr1x,gradxr2x,gradyr2x,gradxr1y,gradyr1y,gradxr2y,gradyr2y,gradxrx,gradxry,gradyrx,gradyry, graduxnew, graduynew, gradvxnew, gradvynew, temp;
    
    double t; // Time variable
    
    
    ofstream myfile, myfile1, myfile2, myfile3, details, uu, vv;
    // Removing existing files with these names if any
    int i = remove("phi.txt");
    i = remove("psix.txt");
    i = remove("psiy.txt");
    i = remove("psixy.txt");
    i = remove("details.txt");
    i = remove("Velocity_x.txt");
    i = remove("Velocity_y.txt");
    
    // Opening new files
    myfile.open("phi.txt", ios::out);
    myfile1.open("psix.txt", ios::out);
    myfile2.open("psiy.txt", ios::out);
    myfile3.open("psixy.txt", ios::out);
    details.open("details.txt", ios::out);
    uu.open("Velocity_x.txt", ios::out);
    vv.open("Velocity_y.txt",ios::out);

#endif
