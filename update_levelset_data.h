// case selection

#ifndef _UPDATE_LEVELSET_DATA_H
#define _UPDATE_LEVELSET_DATA_H

        switch(count){
            case 1:{
                tempphi[tempindex_y][tempindex_x] = hp(phi, psix, psiy, psixy, xadv, yadv, xo, yo, dx, dy);
                rootpsix = hermx(phi,psix,psiy,psixy,xadv,yadv,xo,yo,dx, dy);
                rootpsiy = hermy(phi,psix,psiy,psixy,xadv,yadv,xo,yo,dx, dy);

                if(strcmp("Heuns",psischeme) == 0){
                                        
                // The following derivatives are the velocity gradients required for the third equation in the mehtod of characteristics for GRADIENTS
                    graduxnode = gradUx(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                    graduynode = gradUy(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                    gradvxnode = gradVx(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                    gradvynode = gradVy(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                                        
                    graduxroot = gradUx(xadv, yadv, t * dt, T_period);
                    graduyroot = gradUy(xadv, yadv, t * dt, T_period);
                    gradvxroot = gradVx(xadv, yadv, t * dt, T_period);
                    gradvyroot = gradVy(xadv, yadv, t * dt, T_period);
                                        
                // Heun's RK-2 method - first step - for x & y
                    temppsix[tempindex_y][tempindex_x] = rootpsix - dt * (graduxroot * rootpsix + gradvxroot * rootpsiy);
                    temppsiy[tempindex_y][tempindex_x] = rootpsiy - dt * (graduyroot * rootpsix + gradvyroot * rootpsiy);
                                        
                // Heun's RK-2 method - second step - for x & y
                    temppsix[tempindex_y][tempindex_x] = rootpsix - 0.5 * dt * (graduxroot * rootpsix + gradvxroot * rootpsiy + graduxnode * temppsix[tempindex_y][tempindex_x] + gradvxnode * temppsiy[tempindex_y][tempindex_x]);
                    temppsiy[tempindex_y][tempindex_x] = rootpsiy - 0.5 * dt * (graduyroot * rootpsix + gradvyroot * rootpsiy + graduynode * temppsix[tempindex_y][tempindex_x] + gradvynode * temppsiy[tempindex_y][tempindex_x]);
                }
                                    
                else if(strcmp("SuperConsistent",psischeme) == 0){                                        
                    graduxnode = gradUx(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                    graduynode = gradUy(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                    gradvxnode = gradVx(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                    gradvynode = gradVy(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                                        
                    gradxr1x = 1 - dt * graduxnode;
                    gradxr1y = 0 - dt * graduynode;
                    gradyr1x = 0 - dt * gradvxnode;
                    gradyr1y = 1 - dt * gradvynode;
                                        
                        if(strcmp("RK3",backtrace_scheme) == 0){
                                            
                            gradux1 = gradUx(xr1, yr1, t * dt, T_period);
                            graduy1 = gradUy(xr1, yr1, t * dt, T_period);
                            gradvx1 = gradVx(xr1, yr1, t * dt, T_period);
                            gradvy1 = gradVy(xr1, yr1, t * dt, T_period);
                                            
                            gradux2 = gradUx(xr2, yr2, (t * dt + dt/2.0), T_period);
                            graduy2 = gradUy(xr2, yr2, (t * dt + dt/2.0), T_period);
                            gradvx2 = gradVx(xr2, yr2, (t * dt + dt/2.0), T_period);
                            gradvy2 = gradVy(xr2, yr2, (t * dt + dt/2.0), T_period);
                                            
                            gradxr2x = 1 - dt/4.0 * (graduxnode + gradxr1x * gradux1 + gradyr1x * graduy1);
                            gradxr2y = 0 - dt/4.0 * (graduynode + gradxr1y * gradux1 + gradyr1y * graduy1);
                            gradyr2x = 0 - dt/4.0 * (gradvxnode + gradxr1x * gradvx1 + gradyr1x * gradvy1);
                            gradyr2y = 1 - dt/4.0 * (gradvynode + gradxr1y * gradvx1 + gradyr1y * gradvy1);
                                            
                            gradxrx = 1 - dt/6.0 * (graduxnode + gradxr1x * gradux1 + gradyr1x * graduy1 + 4 * (gradxr2x * gradux2 + gradyr2x * graduy2));
                            gradxry = 0 - dt/6.0 * (graduynode + gradxr1y * gradux1 + gradyr1y * graduy1 + 4 * (gradxr2y * gradux2 + gradyr2y * graduy2));
                            gradyrx = 0 - dt/6.0 * (gradvxnode + gradxr1x * gradvx1 + gradyr1x * gradvy1 + 4 * (gradxr2x * gradvx2 + gradyr2x * gradvy2));
                            gradyry = 1 - dt/6.0 * (gradvynode + gradxr1y * gradvx1 + gradyr1y * gradvy1 + 4 * (gradxr2y * gradvx2 + gradyr2y * gradvy2));
                        }

                        else {
                            gradxrx = gradxr1x;
                            gradxry = gradxr1y;
                            gradyrx = gradyr1x;
                            gradyry = gradyr1y;
                        }  // end of if(strcmp("RK3",backtrace_scheme)                                    
                                      
                        // Super consistent method -> psi(n+1) = grad(xroot,yroot) * psiroot(n)
                        temppsix[tempindex_y][tempindex_x] = (gradxrx * rootpsix + gradyrx * rootpsiy);
                        temppsiy[tempindex_y][tempindex_x] = (gradxry * rootpsix + gradyry * rootpsiy);
                }   // end of   else if(strcmp("SuperConsistent",psischeme) == 0)

                 tracker[tempindex_y][tempindex_x] = 1;
                 break;
            } // end of case 1

            case 2:{
                rootpsix = hermx(phi,psix,psiy,psixy,xadv, y[tempindex_y],xo,yo,dx, dy);
                rootpsiy = temppsiy[tempindex_y][tempindex_x];
                tempphi[tempindex_y][tempindex_x] = hp(phi, psix, psiy, psixy, xadv, y[tempindex_y], xo, yo, dx, dy) - dt * Vely(x[tempindex_x], y[tempindex_y], t * dt, T_period) * rootpsiy;
                
                // The following derivatives are the velocity gradients required for the third equation in the mehtod of characteristics for GRADIENTS                                  
                graduxnode = gradUx(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                graduynode = gradUy(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                gradvxnode = gradVx(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                gradvynode = gradVy(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                                    
                graduxnew = gradUx(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                graduynew = gradUy(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                gradvxnew = gradVx(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                gradvynew = gradVy(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                                    
                temp = rootpsix - dt * (gradvxnode * rootpsiy + Vely(x[tempindex_x], y[tempindex_y], t * dt, T_period) * temppsixy[tempindex_y][tempindex_x] + graduxnode * temppsix[tempindex_y][tempindex_x]);
                                    
                //HEUNS
                temppsix[tempindex_y][tempindex_x] = rootpsix - (dt/2.0) * (gradvxnode * rootpsiy + Vely(x[tempindex_x], y[tempindex_y], t * dt, T_period) * temppsixy[tempindex_y][tempindex_x] + graduxnode * temppsix[tempindex_y][tempindex_x] + gradvxnew * rootpsiy + Vely(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period) * temppsixy[tempindex_y][tempindex_x] + graduxnew * temp);
                tracker[tempindex_y][tempindex_x] = 2;
                break;
            }   // end of case 2

            case 3:{                                    
                rootpsix = temppsix[tempindex_y][tempindex_x];
                rootpsiy = hermy(phi,psix,psiy,psixy,x[tempindex_x],yadv,xo,yo,dx, dy);
                tempphi[tempindex_y][tempindex_x] = hp(phi, psix, psiy, psixy, x[tempindex_x],yadv, xo, yo, dx, dy) - dt * Velx(x[tempindex_x], y[tempindex_y], t * dt, T_period) * rootpsix;
                                    
                // The following derivatives are the velocity gradients required for the third equation in the mehtod of characteristics for GRADIENTS                                 
                graduxnode = gradUx(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                graduynode = gradUy(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                gradvxnode = gradVx(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                gradvynode = gradVy(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                                    
                graduxnew = gradUx(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                graduynew = gradUy(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                gradvxnew = gradVx(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                gradvynew = gradVy(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                                   
                temp = rootpsiy - dt * (graduynode * rootpsix + Velx(x[tempindex_x], y[tempindex_y], t * dt, T_period) * temppsixy[tempindex_y][tempindex_x] + gradvynode * temppsiy[tempindex_y][tempindex_x]);
                                    
                //HEUNS
                temppsiy[tempindex_y][tempindex_x] = rootpsiy - (dt/2.0) * (graduynode * rootpsix + Velx(x[tempindex_x], y[tempindex_y], t * dt, T_period) * temppsixy[tempindex_y][tempindex_x] + gradvynode * temppsiy[tempindex_y][tempindex_x] + graduynew * rootpsix + Velx(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period) * temppsixy[tempindex_y][tempindex_x] + gradvynew * temp);
                tracker[tempindex_y][tempindex_x] = 3;
                break;
            }   // end of case 3

            case 4:{
                tempphi[tempindex_y][tempindex_x] = tempphi[tempindex_y][tempindex_x] - dt * (Velx(x[tempindex_x], y[tempindex_y], t * dt, T_period) * temppsix[tempindex_y][tempindex_x] + Vely(x[tempindex_x], y[tempindex_y], t * dt, T_period) * temppsiy[tempindex_y][tempindex_x]);
                tracker[tempindex_y][tempindex_x] = 4;
                break;
            } //end of case4

            default:{break;}
        //---------------------------------------------------------------------------------------------------------
        }   // end of switch loop

      count = 0;

#endif
