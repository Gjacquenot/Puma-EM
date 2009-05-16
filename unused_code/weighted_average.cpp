/***************************************************************************
 * weighted_average.cpp  Acceleration of Sommerfeld integral tails
 *                       Based upon Michalski paper, AP-IEEE 1998 
 *
 * Copyright (C) 2000-2005 Idesbald van den Bosch <vandenbosch@emic.ucl.ac.be>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * Suggestions:          vandenbosch@emic.ucl.ac.be
 * Bugs:                 vandenbosch@emic.ucl.ac.be
 *
 * For more information, please see the ... Home Page:
 *
 ****************************************************************************/
#include <complex>
#include <blitz/array.h>

using namespace blitz;

complex<double> weighted_average (const int R, const Array<double, 1>& xi, const Array<complex<double>, 1>& u_i) {

  // notations are the same as in Michalski 98, 
  // "extrapolation methods for Sommerfeld integral tails" 
  int j, l, k;
  complex<double> S_accel;
  if ((R < 2) || (u_i (R) == 0.0))  S_accel = u_i(1); // avoid dividing by 0
  else {
    Array<complex<double>, 1> S_n (Range(1,R));
    for (j=1 ; j<=R ; j++) S_n(j) = sum(u_i(Range(1,j)));
    Array<double, 1> eta_0 (xi(Range(2, R))/xi(Range(1, R-1)));
    Array<complex<double>, 1> omega (u_i(Range(1, R-1))/u_i(Range(2, R)));
    Array<complex<double>, 1> eta_k (complex<double> (1.0, 0.0) * eta_0);
    Array<complex<double>, 2> S_n_k (Range(1, R), Range(1, R));
    S_n_k = 0.0;
    S_n_k(1, Range::all()) = S_n(Range(1, R));

    for (l=2 ; l<=R ; l++) {
      k = l-2;
      eta_k = -omega*pow(eta_0, k);
      S_n_k(l, Range(1, R-l+1)) = (S_n_k(l-1, Range(1, R-l+1)) + eta_k(Range(1, R-l+1)) * S_n_k(l-1, Range(2, R-l+2)))/(1.0 + eta_k(Range(1, R-l+1)));
    }

    S_accel = S_n_k(R,1);

  }
  return S_accel;
}


