//
//  TimeSteppingMethods.h
//
//
//  Created by Raunak Bardia on 12/01/17.
//
//

namespace galsfunctions
{
    void advection_point(vectorarray& x, vectorarray& y, gridarray& xadv, gridarray &yadv, int t,
            double dt, double T_period, char backtrace_scheme[])
    {
        double c1 = (1/6.0);
        double c2 = (1/6.0);
        double c3 = (2/3.0);    //RK-3 constants
        
        int nx = x.size();
        int ny = y.size();
        xadv.resize(ny,vectorarray(nx,0.0));
        yadv.resize(ny,vectorarray(nx,0.0));
        
        for(int i = 0; i < ny; i++){  // loop for y - rows in the array
            for(int j = 0; j < nx; j++){  // loop for x - columns in the array
                
                double ux = Velx(x[j], y[i], (t + 1) * dt, T_period);
                double vy = Vely(x[j], y[i], (t + 1) * dt, T_period);
                
                // Advected Points - 1 Step
                xadv[i][j] = x[j] - ux * dt;
                yadv[i][j] = y[i] - vy * dt;
                
                if(strcmp("RK3",backtrace_scheme) == 0)
                {
                    double ux1 = Velx(xadv[i][j], yadv[i][j], t * dt, T_period);
                    double vy1 = Vely(xadv[i][j], yadv[i][j], t * dt, T_period);
                    
                    // 2 Step
                    xadv[i][j] = x[j] - ((dt)/2.0 * (0.5 * ux + 0.5 * ux1));
                    yadv[i][j] = y[i] - ((dt)/2.0 * (0.5 * vy + 0.5 * vy1));
                    
                    double ux2 = Velx(xadv[i][j], yadv[i][j],(t * dt + dt/2.0), T_period);
                    double vy2 = Vely(xadv[i][j], yadv[i][j],(t * dt + dt/2.0), T_period);
                    
                    // 3 Step
                    xadv[i][j] = x[j] - (dt * (c1 * ux + c2 * ux1 + c3 * ux2));
                    yadv[i][j] = y[i] - (dt * (c1 * vy + c2 * vy1 + c3 * vy2));
                }
            }
        }
    }
}
