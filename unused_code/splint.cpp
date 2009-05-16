#include <complex>
#include <blitz/array.h>

using namespace blitz;

#include "splint.h"

void spline_y2 (Array<complex<double>, 1>& y2, const Array<double, 1>& x0, const Array<complex<double>, 1>& y0) {
  // x0[0..n-1], y0[0..n-1]
  // computation of the second derivatives y2 of input vector y0
  // spline is supposed "natural": second derivative is zero at extremities of domain
  // this code is inspired from "Numerical Recipes in C"
  int j, n = x0.size();
  if (n<4) {
    cout << "spline_y2(): too small vector for spline calculation" << endl;
    y2 = 0.0;
  }
  else {
    double sig;
    complex<double> p;
    Array<complex<double>, 1> u(n-1); // u[0..n-2]

    y2(0) = 0.0; // second derivative at origin 
    u(0) = 0.0;
    for (j=1 ; j<=n-2 ; j++) {
      sig = (x0(j)-x0(j-1))/(x0(j+1)-x0(j-1));
      p = sig*y2(j-1)+2.0;
      y2(j) = (sig-1.0)/p;
      u(j) = (y0(j+1)-y0(j))/(x0(j+1)-x0(j)) - (y0(j)-y0(j-1))/(x0(j)-x0(j-1));
      u(j) = (6.0*u(j)/(x0(j+1)-x0(j-1))-sig*u(j-1))/p;
    }
    y2(n-1) = 0.0; // second derivative at end
    for (j=n-2 ; j>=0 ; j--)  y2(j) = y2(j)*y2(j+1)+u(j);
  }
}

void splint_1D_point (complex<double>& f, const double x, const Array<double, 1>& x0, const Array<complex<double>, 1>& f0, const Array<complex<double>, 1>& f2) {
  // OUTPUT f is the interpolated value at abscissa x
  // INPUT x is the desired abscissa
  // INPUT vector x0 contains the original abscissas
  // INPUT vector f0 is the original function, given at abscissas x0
  // INPUT f2 must be precomputed (for gaining time)
  int j, ind_x, ind_x_sup;
  double a, b, h;
  ind_x = find_index(x, x0);
  ind_x_sup = ind_x + 1;

  h = x0(ind_x_sup)-x0(ind_x);
  a = (x0(ind_x_sup)-x)/h;
  b = (x-x0(ind_x))/h;
  f = a*f0(ind_x) + b*f0(ind_x_sup) + ((a*a*a-a)*f2(ind_x) + (b*b*b-b)*f2(ind_x_sup))*(h*h)/6.0;
}

void splint_1D_vector (Array<complex<double>, 1>& f, const Array<double, 1>& x, const Array<double, 1>& x0, const Array<complex<double>, 1>& f0, const Array<complex<double>, 1>& f2) {
  // OUTPUT f is the interpolated vector at points x
  // INPUT vector x contains the desired abscissas
  // INPUT vector x0 contains the original abscissas
  // INPUT vector f0 is the original function, given at abscissas x0
  // INPUT f2 must be precomputed (for gaining time)
  int j, N_x = x.size();
  f.resize(N_x);
  for (j=0 ; j<N_x ; j++) splint_1D_point(f(j), x(j), x0, f0, f2);
}

void decimate_2D_splint (Array<double, 1>& x1_int, Array<double, 1>& x2_int, Array<complex<double>, 2>& Y_int, const Array<double, 1>& x1, const Array<double, 1>& x2, const Array<complex<double>, 2>& Y) {
  // OUTPUT:
  // x1_int, x2_int: the first and second dimension abscissas for the
  // interpolated array Y_int. These abscissas will be equispaced
  // Y_int: the decimated array
  // 
  // INPUT:
  // x1_int, x2_int: the first and second dimension abscissas for the
  // input array Y. These abscissas are supposed to be equispaced
  // Y: the array to be decimated
  Range all = Range::all();
  int j, N_x1 = x1.size(), N_x2 = x2.size(), N_x1_int, N_x2_int;
  Y_int.resize(2*N_x1 - 1, 2*N_x2 - 1);
  Y_int(Range(0, toEnd, 2), Range(0, toEnd, 2)) = Y; // these values
						     // do not have to
						     // be interpolated

  // we now resize the x1_int and x2_int arrays (interpolation abscissas)
  x1_int.resize(2*N_x1 - 1), x2_int.resize(2*N_x2 - 1);
  N_x1_int = x1_int.size(), N_x2_int = x2_int.size();
  double D_x1_int = (x1(N_x1-1) - x1(0))/(N_x1_int-1), D_x2_int = (x2(N_x2-1) - x2(0))/(N_x2_int-1);
  for (j=0 ; j<N_x1_int ; j++) x1_int(j) = x1(0) + j * D_x1_int;
  for (j=0 ; j<N_x2_int ; j++) x2_int(j) = x2(0) + j * D_x2_int;

  Array<complex<double>, 1> y2_tmp, y_tmp; // temporary arrays that will 
              				   // be used in the computations
  // we first interpolate following the 1st dimension (i.e. along the
  // columns of the array)
  y2_tmp.resize(N_x1); 
  y_tmp.resize(N_x1-1);
  for (j=0 ; j<N_x2 ; j++) {
    spline_y2(y2_tmp, x1, Y(all, j));
    splint_1D_vector(y_tmp, x1_int(Range(1, toEnd, 2)), x1, Y(all, j), y2_tmp);
    Y_int(Range(1, toEnd, 2), 2*j) = y_tmp;
  }
  // we then interpolate following x2 (along the lines). 
  // Since we now must use Y_int for this purpose, pay attention to the fact
  // that we now have N_x2_int = (2*N_x2 - 1) values following x2!!
  y2_tmp.resize(N_x2); 
  y_tmp.resize(N_x2-1);
  for (j=0 ; j<N_x1_int ; j++) {
    spline_y2(y2_tmp, x2, Y_int(j, Range(0, toEnd, 2)));
    splint_1D_vector(y_tmp, x2_int(Range(1, toEnd, 2)), x2, Y_int(j, Range(0, toEnd, 2)), y2_tmp);
    Y_int(j, Range(1, toEnd, 2)) = y_tmp;
  }
}


// NOT WORKING YET!!!
// void decimate_3D_splint (Array<double, 1>& x1_int, Array<double, 1>& x2_int, Array<double, 1>& x3_int, Array<complex<double>, 3>& M_int, const Array<double, 1>& x1, const Array<double, 1>& x2, const Array<double, 1>& x3, const Array<complex<double>, 3>& M, const int N_decim_x1, const int N_decim_x2, const int N_decim_x3) {
//   // this function decimates a complex 3D array using cubic splines
//   // x1, x2, x3 are the OUTPUT abscissas in the 1st, 2nd, 3rd dimensions
//   // respectively
//   // M will be the OUTPUT array
//   // INPUT N_decim_xi are the decimation factors: 1, 2, etc...
//   // INPUT x1_int, x2_int, x3_int 



//   firstIndex k;
//   Array<complex<double>, 1> f2, f0, f;
//   int m, n, p, N_x1 = x1.size(), N_x2 = x2.size(), N_x3 = x3.size();
//   int N_x1_int = (N_x1-1)*N_decim_x1 + N_x1, N_x2_int = (N_x2-1)*N_decim_x2 + N_x2, N_x3_int = (N_x3-1)*N_decim_x3 + N_x3;
//   x1_int.resize(N_x1_int);
//   x2_int.resize(N_x2_int);
//   x3_int.resize(N_x3_int);
//   double Delta_x1_int, Delta_x2_int, Delta_x3_int; 
//   Delta_x1_int = (N_x1>1) ? (x1 (1)-x1 (0))/(N_decim_x1+1) : 0.0;
//   Delta_x2_int = (N_x2>1) ? (x2 (1)-x2 (0))/(N_decim_x2+1) : 0.0;
//   Delta_x3_int = (N_x3>1) ? (x3 (1)-x3 (0))/(N_decim_x3+1) : 0.0;
//   x1_int = x1 (0) + k*Delta_x1_int;
//   x2_int = x2 (0) + k*Delta_x2_int;
//   x3_int = x3 (0) + k*Delta_x3_int;

//   M_int.resize(N_x1_int, N_x2_int, N_x3_int);
//   M_int = 0.0;
//   M_int (Range (0, toEnd, N_decim_x1+1), Range (0, toEnd, N_decim_x2+1), Range (0, toEnd, N_decim_x3+1)) = M;
//   if (N_decim_x3>0) {
//     for (m=0 ; m<N_x1 ; m++) {
//       for (n=0 ; n<N_x2 ; n++) { // we interpolate each (m, n, all) of M
//         f0.resize (N_x3);
//         f0 = M_int (m*(N_decim_x1+1), n*(N_decim_x2+1), Range (0, toEnd, N_decim_x3+1));
//         spline_y2 (f2, x3, f0);
//         f.resize (N_x3_int);
//         splint_1D_vector (x3, f0, f2, x3_int, f);
//         M_int (m*(N_decim_x1+1), n*(N_decim_x2+1), Range::all()) = f;
//       }
//     }
//   }
//   if (N_decim_x2>0) {
//     for (m=0 ; m<N_x1 ; m++) {
//       for (p=0 ; p<N_x3_int ; p++) { // we interpolate each (m*(N_decim_x1+1), all, p) of M_int
//         f0.resize (N_x2);
//         f0 = M_int (m*(N_decim_x1+1), Range(0,toEnd,N_decim_x2+1), p);
//         spline_y2 (f2, x2, f0);
//         f.resize (N_x2_int);
//         splint_1D_vector (x2, f0, f2, x2_int, f);
//         M_int (m*(N_decim_x1+1), Range::all(), p) = f;
//       }
//     }
//   }
//   if (N_decim_x1>0) {
//     for (n=0 ; n<N_x2_int ; n++) {
//       for (p=0 ; p<N_x3_int ; p++) { // we interpolate each (all, n, p) of M_int
//         f0.resize (N_x1);
//         f0 = M_int (Range(0,toEnd,N_decim_x1+1), n, p);
//         spline_y2 (f2, x1, f0);
//         f.resize (N_x1_int);
//         splint_1D_vector (x1, f0, f2, x1_int, f);
//         M_int (Range::all(), n, p) = f;
//       }
//     }
//   }
// }
