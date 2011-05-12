#include <string>
#include <blitz/array.h>

using namespace blitz;

#include "GL.h"

void integr_1D_X_W(Array<double, 1>& X, Array<double, 1>& W, const double & a, const double & b, const int N_points, const string & METHOD)
{
  double h;
  X.resize(N_points);
  W.resize(N_points);

  if (METHOD == "TRAP")
    {
      h = (b-a)/(N_points-1); // trapezoidal rule
      for (int j=0 ; j<N_points ; j++) X(j) = a + j*h;
      W = h;
      W(0) /= 2.0;
      W(N_points-1) /= 2.0;
    }
  else if (METHOD == "PONCELET")
    {
      h = (b-a)/N_points; // mid-point method
      for (int j=0 ; j<N_points ; j++) X(j) = a + j*h + h/2.0;
      W = h;
    }
  else if (METHOD == "GAUSSL")
    {
      double Dx, center;
      Array<double, 1> XGL, WGL;
      if (N_points<=20) {
        Gauss_Legendre(XGL, WGL, N_points); 
        Dx = 0.5 * (b - a);
        center = 0.5 * (b + a);
        for (int j=0 ; j<N_points ; j++) X(j) = center + Dx * XGL(j); 
        W = abs(Dx) * WGL;
      }
    }

  else
    {
      cout << "integr_1D_X_W(): You have not entered a correct method. It should be TRAP, PONCELET or GAUSSL" << endl;
      exit(1);
    }
}
