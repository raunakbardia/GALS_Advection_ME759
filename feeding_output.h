// Feeding phi, psix, psiy and psixy values in their respective files

#ifndef _FEEDING_OUTPUT_H
#define _FEEDING_OUTPUT_H

        int check = t+1;
        if(check % printstep == 0){
            for(int i = 0; i < ny; i++){
                
                myfile << std::fixed << std::setprecision(10) <<  mphi[i][0];
                myfile1 << std::fixed << std::setprecision(10) << mpsix[i][0];
                myfile2 << std::fixed << std::setprecision(10) << mpsiy[i][0];
                myfile3 << std::fixed << std::setprecision(10) << mpsixy[i][0];
                
                u[i][0] = Velx(x[0],y[i],(t+1)*dt, T_period);
                uu << std::fixed << std::setprecision(10) << u[i][0];

                v[i][0] = Vely(x[0],y[i],(t+1)*dt, T_period);
                vv << std::fixed << std::setprecision(10) << v[i][0];
                
                for(int j = 1; j < nx; j++){
                    myfile << ",";      myfile << std::fixed << std::setprecision(10) << mphi[i][j];
                    myfile1 << ",";     myfile1 << std::fixed << std::setprecision(10) << mpsix[i][j];
                    myfile2 << ",";     myfile2 << std::fixed << std::setprecision(10) << mpsiy[i][j];
                    myfile3 << ",";     myfile3 << std::fixed << std::setprecision(10) << mpsixy[i][j];
                    
                    u[i][j] = Velx(x[j],y[i],(t+1)*dt, T_period);
                    uu << "," << std::fixed << std::setprecision(10) << u[i][j];

                    v[i][j] = Vely(x[j],y[i],(t+1)*dt, T_period);
                    vv << "," << std::fixed << std::setprecision(10) << v[i][j];
                }
                myfile << '\n';            myfile1 << '\n';                myfile2 << '\n';                myfile3 << '\n';             uu << '\n';                vv << '\n';
            }
            myfile << '\n';                myfile1 << '\n';                myfile2 << '\n';                myfile3 << '\n';             uu << '\n';                vv << '\n';
        }

#endif
