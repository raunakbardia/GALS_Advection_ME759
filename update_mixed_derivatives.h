// Update the mixed derivatives now for the remaining grid points

#ifndef _UPDATE_MIXED_DERIVATIVES_H
#define _UPDATE_MIXED_DERIVATIVES_H

double dxy1, dxy2;
double d1, d2, d3, d4;  //Temporary variable to compute mixed derivatives


for(int k = 0; k < ny; k++){
    for(int l = 0; l < nx; l++){
        if ((k == 0 || k == ny - 1) && (l != 0 && l != nx - 1)){
            temppsixy[k][l] = (temppsiy[k][l+1] - temppsiy[k][l-1])/(2 * dx);
        }
        else
            if ((k != 0 && k != ny - 1) && (l == 0 || l == nx - 1)){
                temppsixy[k][l] = (temppsix[k+1][l] - temppsix[k-1][l])/(2 * dy);
            }
            else
                if((k == 0 || k == ny - 1) && (l == 0 || l == nx - 1)){
                    if(k == 0 && l == 0){
                        d1 = (temppsiy[0][1] - temppsiy[0][0])/dx;
                        d2 = (temppsix[1][0] - temppsix[0][0])/dy;
                        d3 = (temppsix[1][1] - temppsix[0][1])/dy;
                        d4 = (temppsiy[1][1] - temppsiy[1][0])/dx;
                        temppsixy[k][l] = 0.75 * (d1 + d2) - 0.25 * (d3 + d4);
                    }
                    else if(k == 0 && l == nx-1){
                        d1 = (temppsiy[0][nx-1] - temppsiy[0][nx-2])/dx;
                        d2 = (temppsix[1][nx-2] - temppsix[0][nx-2])/dy;
                        d3 = (temppsix[1][nx-1] - temppsix[0][nx-1])/dy;
                        d4 = (temppsiy[1][nx-1] - temppsiy[1][nx-2])/dx;
                        temppsixy[k][l] = 0.75 * (d1 + d3) - 0.25 * (d2 + d4);
                        
                    }
                    else if(k == ny-1 && l == 0){
                        d1 = (temppsiy[ny-2][1] - temppsiy[ny-2][0])/dx;
                        d2 = (temppsix[ny-1][0] - temppsix[ny-2][0])/dy;
                        d3 = (temppsix[ny-1][1] - temppsix[ny-2][1])/dy;
                        d4 = (temppsiy[ny-1][1] - temppsiy[ny-1][0])/dx;
                        temppsixy[k][l] = 0.75 * (d2 + d4) - 0.25 * (d3 + d1);
                        
                    }
                    else if(k == ny-1 && l == nx-1){
                        d1 = (temppsiy[ny-2][nx-1] - temppsiy[ny-2][nx-2])/dx;
                        d2 = (temppsix[ny-1][nx-2] - temppsix[ny-2][nx-2])/dy;
                        d3 = (temppsix[ny-1][nx-1] - temppsix[ny-2][nx-1])/dy;
                        d4 = (temppsiy[ny-1][nx-1] - temppsiy[ny-1][nx-2])/dx;
                        temppsixy[k][l] = 0.75 * (d3 + d4) - 0.25 * (d1 + d2);
                    }
                    
                }
                else{
                    dxy1 = (temppsiy[k][l+1] - temppsiy[k][l-1])/(2 * dx);
                    dxy2 = (temppsix[k+1][l] - temppsix[k-1][l])/(2 * dy);
                    temppsixy[k][l] = (dxy1 + dxy2)/2.0;
                }
    }
}

#endif
