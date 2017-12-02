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
    
    double d1, d2, d3, d4;  //Temporary variable to compute mixed derivatives
    double np = nx * ny; // Number of node points
    
    // Removing existing files with these names if any
    int i = remove("phi.txt");
    i = remove("psix.txt");
    i = remove("psiy.txt");
    i = remove("psixy.txt");
    i = remove("details.txt");
    i = remove("Velocity_x.txt");
    i = remove("Velocity_y.txt");
    
#endif
