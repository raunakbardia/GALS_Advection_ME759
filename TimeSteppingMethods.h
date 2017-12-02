//
//  TimeSteppingMethods.h
//
//
//  Created by Raunak Bardia on 12/01/17.
//
//

namespace galsfunctions
{
void advection_point(vectorarray& x, vectorarray& y, int tempindex_x, int tempindex_y, int t, double dt, double T_period, double& xadv, double& yadv, char backtrace_scheme[])
{ 
    	double c1 = (1/6.0);
    	double c2 = (1/6.0);
    	double c3 = (2/3.0);    //RK-3 constants
	double ux = Velx(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
        double vy = Vely(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                      
        // Advected Points - 1 Step
        xadv = x[tempindex_x] - ux * dt;
        yadv = y[tempindex_y] - vy * dt;

        if(strcmp("RK3",backtrace_scheme) == 0)
	{                 
        	double ux1 = Velx(xadv, yadv, t * dt, T_period);
                double vy1 = Vely(xadv, yadv, t * dt, T_period);
                            
                // 2 Step
                xadv = x[tempindex_x] - ((dt)/2.0 * (0.5 * ux + 0.5 * ux1));
                yadv = y[tempindex_y] - ((dt)/2.0 * (0.5 * vy + 0.5 * vy1));
                            
                double ux2 = Velx(xadv, yadv,(t * dt + dt/2.0), T_period);
                double vy2 = Vely(xadv, yadv,(t * dt + dt/2.0), T_period);
                            
                // 3 Step
                xadv = x[tempindex_x] - (dt * (c1 * ux + c2 * ux1 + c3 * ux2));
                yadv = y[tempindex_y] - (dt * (c1 * vy + c2 * vy1 + c3 * vy2));
	}
}
}                     
