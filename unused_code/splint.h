#ifndef SPLINT_H
#define SPLINT_H
using namespace blitz;

template <typename T>
int find_index(const T x, const Array<T, 1>& x_i) {
  // a log(n) routine for finding the index of a x in an array x_i of size n
  // the array is supposed to be monotonically increasing
  int ind_inf, ind_sup, ind_mid;
  const int n = x_i.size();
  if ( (x < x_i(0)) || (x > x_i(n-1)) ) {
    cout << "find_index(): x outside range. x_i = " << x_i << ", x = " << x << endl;
    exit(1);
  }
  else {
    ind_inf = 0;
    ind_sup = n-1;
    while(ind_sup-ind_inf > 1) {
      ind_mid = (ind_sup+ind_inf)/2;
      if (x <= x_i(ind_mid)) ind_sup = ind_mid;
      else ind_inf = ind_mid;
    }
    return ind_inf;
  }
}

void spline_y2 (Array<complex<double>, 1>& y2, const Array<double, 1>& x0, const Array<complex<double>, 1>& y0);

void splint_1D_point (complex<double>& f, const double x, const Array<double, 1>& x0, const Array<complex<double>, 1>& f0, const Array<complex<double>, 1>& f2);

void splint_1D_vector (Array<complex<double>, 1>& f, const Array<double, 1>& x, const Array<double, 1>& x0, const Array<complex<double>, 1>& f0, const Array<complex<double>, 1>& f2);

void decimate_2D_splint (Array<double, 1>& x1_int, Array<double, 1>& x2_int, Array<complex<double>, 2>& Y_int, const Array<double, 1>& x1, const Array<double, 1>& x2, const Array<complex<double>, 2>& Y);

#endif
