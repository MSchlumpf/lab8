// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx,
	    double* const k1, double* const k2, double* const k3, double* const k4);
void btheta(double* const b, double theta);
void interstep(double* const b, double* const k1, double* const k2, double* const k3, double* const k4,
	       double* const y0, double* const yn, const double dx);

//--------------------
using namespace std;
//--------------------

int main(void)
{
  ofstream out("solution");
  const int dim = 2;
  double dx = 0.1,x=0;
  const double L = 100;
  
  
  for (double p0=0.1; p0 < 5; p0 += 0.1){
    x = 0;
    double y0[dim] = {p0, 0};
    double yn[dim];
    double k1[dim], k2[dim], k3[dim], k4[dim];
//   out << x << "\t" << y0[0] << "\t" << y0[1] << "\t" << endl;

	while(x<=L)
	{
		x += dx;
		RKstep(yn, y0, x, dx, k1, k2, k3, k4);
		if (y0[1]>0 && yn[1]<0)
		  break;
		
		for(int i=0; i<dim; i++) y0[i] = yn[i];
// 		out << x << "\t" << y0[0] << "\t" << y0[1] << endl;
	}
//   out << x << "\t" << yn[0] << "\t" << yn[1] << endl;
//   out.close();
	//starting interpolation
	double b[4], c;
	
	double x0 = x-dx;
	
	double xl = x0;
	double xr = x;
	double xm = (xl+xr)/2;
	double theta = (xm -x0) / dx;
      while (abs(xl-xr) > 1E-5){
	
	btheta(b, theta);
	interstep(b, k1, k2, k3, k4, y0, yn, dx);
	c = y0[1]*yn[1];
	
	if (c > 0){
	  xl = xm;
	}
	
	else
	  xr = xm; 
	
	xm = (xl+xr)/2;
	theta = (xm - x0) / dx;
	  
      }
    out << p0 << "\t" << x0 + theta*dx << endl;   
  }
  out.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx, double* const k1, double* const k2, double* const k3, double* const k4)
{
	const int dim = 2;

	for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1, x);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dx * k1[i];
	f(k2, x+0.5*dx);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dx * k2[i];
	f(k3, x+0.5*dx);

	for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * k3[i];
	f(k4,  x+dx);

	for(int i=0;i<dim; i++)
	  yn[i] = y0[i] + 1./6.*dx*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}
//-------------------
// Lorenz model
void f(double* const y0, const double x){
	
	double y[2] = { y0[0], y0[1] };

	y0[0] = y[1];
	y0[1] = -y[0]/sqrt(1+y[0]*y[0]);
}
//-------------------
void btheta(double* const b, double theta){
  b[0] = theta - 3*theta*theta/2 + 2* theta*theta*theta/3;
  b[1] = theta*theta - 2*theta*theta*theta/3;
  b[2] = b[1]; 
  b[3] = -theta*theta/2 + 2*theta*theta*theta/3;
}
//-------------------
void interstep(double* const b, double* const k1, double* const k2, double* const k3, double* const k4,
	       double* const y0, double* const yn, const double dx){
  const int dim = 2;
  for(int i=0;i<dim; i++)
    yn[i] = y0[i] + dx *(b[0]*k1[i]+b[1]*k2[i]+b[2]*k3[i]+b[3]*k4[i]);
}