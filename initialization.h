// Initializing at t = 0

#ifndef _INITIALIZATION_H
#define _INITIALIZATION_H

for(int i = 0; i < ny; i++){
        
        // Initializing the first element of each row separately to avoid trailing commas in the file
        // Storing these values in the master variable as well as temporary variable matrices
        mphi[i][0] = initialize(x[0], y[i]);
        tempphi[i][0] = mphi[i][0];
        myfile << std::fixed << std::setprecision(10) << mphi[i][0];
        mpsix[i][0] = derivxinit(x[0], y[i]);
        temppsix[i][0] = mpsix[i][0];
        myfile1 << std::fixed << std::setprecision(10) << mpsix[i][0];
        mpsiy[i][0] = derivyinit(x[0], y[i]);
        temppsiy[i][0] = mpsiy[i][0];
        myfile2 << std::fixed << std::setprecision(10) << mpsiy[i][0];
        mpsixy[i][0] = derivxyinit(x[0], y[i]);
        temppsixy[i][0] = mpsixy[i][0];
        myfile3 << std::fixed << std::setprecision(10) << mpsixy[i][0];
        
        u[i][0] = Velx(x[0],y[i],0, T_period);
        uu << std::fixed << std::setprecision(10) << u[i][0];
        v[i][0] = Vely(x[0],y[i],0, T_period);
        vv << std::fixed << std::setprecision(10) << v[i][0];
        
        for(int j = 1; j < nx; j++){
            mphi[i][j] = initialize(x[j], y[i]);
            tempphi[i][j] = mphi[i][j];
            myfile << ",";
            myfile << std::fixed << std::setprecision(10) << mphi[i][j];
            mpsix[i][j] = derivxinit(x[j], y[i]);
            temppsix[i][j] = mpsix[i][j];
            myfile1 << ",";
            myfile1 << std::fixed << std::setprecision(10) << mpsix[i][j];
            mpsiy[i][j] = derivyinit(x[j], y[i]);
            temppsiy[i][j] = mpsiy[i][j];
            myfile2 << ",";
            myfile2 << std::fixed << std::setprecision(10) << mpsiy[i][j];
            mpsixy[i][j] = derivxyinit(x[j], y[i]);
            temppsixy[i][j] = mpsixy[i][j];
            myfile3 << ",";
            myfile3 << std::fixed << std::setprecision(10) << mpsixy[i][j];
            
            u[i][j] = Velx(x[j],y[i],0, T_period);
            uu << "," << std::fixed << std::setprecision(10) << u[i][j];
            v[i][j] = Vely(x[j],y[i],0, T_period);
            vv << "," << std::fixed << std::setprecision(10) << v[i][j];
        }
        myfile << '\n';
        myfile1 << '\n';
        myfile2 << '\n';
        myfile3 << '\n';
        uu << '\n';
        vv << '\n';
    }
#endif
