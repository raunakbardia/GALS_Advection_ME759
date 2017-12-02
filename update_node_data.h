// define all required parameters

#ifndef _UPDATE_NODE_DATA_H
#define _UPDATE_NODE_DATA_H

for(int k = 0; k < 4; k++)
{
    int count = 0;  // Additional variable to make sure that the advection point is available for the relevant node point
    int temp1 = 0, temp2 = 0;       // These temporary variables are updated when the advection point lies between the x and y limits of the cell
    int tempindex_y = i + (-2 * pow(k,3) + 9 * pow(k,2) - 7 * k)/6;   // Temporary Indexes polynomial
    int tempindex_x = j + (2 * pow(k,3) - 9 * pow(k,2) + 10 * k)/3;
    
    if(tracker[tempindex_y][tempindex_x] == 0){ // If the node has already been updated we need not update it again
        
        
        if((xadv[tempindex_y][tempindex_x] >= x[j]) && (xadv[tempindex_y][tempindex_x] <= x[j+1]))
            temp1 = 1;
        if((yadv[tempindex_y][tempindex_x] >= y[i]) && (yadv[tempindex_y][tempindex_x] <= y[i+1]))
            temp2 = 1;
        
        //---------------------------------------------------------------------------------------------------------
        
        //count variable makes sure that the values are updated if they can be through this cell's hermite interpolant itself
        
        // If the advection point lies within the cell then we don't need to change the advection point
        if(temp1 && temp2){
            count = 1;
        }
        // If only the x coordinate is within the cell limits then only update the point if it is on y-boundary
        if(temp1 && !temp2){
            if((yadv[tempindex_y][tempindex_x] < ylim1) || (yadv[tempindex_y][tempindex_x] > ylim2)){
                count = 2;
            }
        }
        // If only y coordinate is within the cell limits then only update the point if it is on x-boundary
        if(!temp1 && temp2){
            if((xadv[tempindex_y][tempindex_x] < xlim1) || (xadv[tempindex_y][tempindex_x] > xlim2)){
                count = 3;
            }
        }
        // If both coordinates dont lie in the cell limits, then the value is not changed for phi and first derivatives
        if(!temp1 && !temp2){
            if(((tempindex_y == 0) || (tempindex_y == ny-1)) && ((tempindex_x == 0) || (tempindex_x == nx-1))){
                count = 4;
            }
        }
        //---------------------------------------------------------------------------------------------------------
        
        if(((tempindex_y == 0) || (tempindex_y==ny-1)) && ((tempindex_x == 0) || (tempindex_x==nx-1))){
            corner = k; // This variable is to check if we have a grid corner in the cell because the psixy value for corners is calculated differently
        }
        
        //This condition and variable has been defined so that these statements need to be written only once otherwise all these statements would have been repeated in all the above if's
        if(count != 0) {
#include "update_levelset_data.h"
        }
    }   // end of if(tracker[tempindex_y][tempindex_x] == 0) loop
}   // end of for(int k = 0; k < 4; k++)

#endif
