#include <string>
#include <blitz/array.h>

using namespace std;

void integr_1D_X_W(blitz::Array<double, 1>& X, blitz::Array<double, 1>& W, const double & a, const double & b, const int N_points, const string & METHOD);
