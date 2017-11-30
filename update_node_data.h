// define all required parameters

#ifndef _UPDATE_NODE_DATA_H
#define _UPDATE_NODE_DATA_H

for(int k = 0; k < 4; k++) 
    {                     
                    int count = 0;  // Additional variable to make sure that the advection point is available for the relevant node point
                    temp1 = 0;           temp2 = 0;
                    tempindex_y = i + (-2 * pow(k,3) + 9 * pow(k,2) - 7 * k)/6;   // Temporary Indexes polynomial
                    tempindex_x = j + (2 * pow(k,3) - 9 * pow(k,2) + 10 * k)/3;
                    
                    if(tracker[tempindex_y][tempindex_x] == 0){ // If the node has already been updated we need not update it again
                        
                        ux = Velx(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                        vy = Vely(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                        
                        // Advected Points - 1 Step
                        xr1 = x[tempindex_x] - ux * dt;
                        yr1 = y[tempindex_y] - vy * dt;

                        if(strcmp("RK3",backtrace_scheme) == 0){
                            
                            ux1 = Velx(xr1, yr1, t * dt, T_period);
                            vy1 = Vely(xr1, yr1, t * dt, T_period);
                            
                            // 2 Step
                            xr2 = x[tempindex_x] - ((dt)/2.0 * (0.5 * ux + 0.5 * ux1));
                            yr2 = y[tempindex_y] - ((dt)/2.0 * (0.5 * vy + 0.5 * vy1));
                            
                            ux2 = Velx(xr2, yr2,(t * dt + dt/2.0), T_period);
                            vy2 = Vely(xr2, yr2,(t * dt + dt/2.0), T_period);
                            
                            // 3 Step
                            xadv = x[tempindex_x] - (dt * (c1 * ux + c2 * ux1 + c3 * ux2));
                            yadv = y[tempindex_y] - (dt * (c1 * vy + c2 * vy1 + c3 * vy2));
                        }

                        else{
                            xadv = xr1;
                            yadv = yr1;
                        }
                      

                        if((xadv >= x[j]) && (xadv <= x[j+1]))
                            temp1 = 1;
                        if((yadv >= y[i]) && (yadv <= y[i+1]))
                            temp2 = 1;
         
                        //---------------------------------------------------------------------------------------------------------

                        //count variable makes sure that the values are updated if they can be through this cell's hermite interpolant itself   

                        // If the advection point lies within the cell then we don't need to change the advection point
                        if(temp1 && temp2){ 
                            count = 1; 
                        }
                        // If only the x coordinate is within the cell limits then only update the point if it is on y-boundary
                        if(temp1 && !temp2){
                            if((yadv < ylim1) || (yadv > ylim2)){
                                count = 2;
                            }
                        }                        
                        // If only y coordinate is within the cell limits then only update the point if it is on x-boundary
                        if(!temp1 && temp2){
                            if((xadv < xlim1) || (xadv > xlim2)){
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
