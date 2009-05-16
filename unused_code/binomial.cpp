/***************************************************************************
 * AIM_S.cpp  computes the inverse of the Vandermonde matrix
 *            Based upon Bleszynski paper, Radio Science September-October 1996
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
#include <iostream>
#include <complex>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

using namespace blitz;

int factorial (const int N) {

  int i, result = 1;
  if (N<0) {
    cout << "Error in factorial: N must be positive" << endl;
    exit(1);
  }
  else {
    for (i = 1 ; i < N+1 ; i++) result *= i;
  }

  return result;

}

int factorial_N_k (const int N, const int k) {

  int i, result;

  if (k>N) {
    cout << "Error in factorial: k must be smaller than N" << endl;
    exit(1);
  }
  else if ((N<0) || (k<0)) {
    cout << "Error in factorial: N and/or k must be positive" << endl;
    exit(1);
  }

  result = 1;

  if ((N>1) && (k<N)) {
    for (i = N-k+1 ; i < N+1 ; i++) result *= i;
  }

  return result;
}

int binomial (const int N, const int k) {

  return factorial_N_k (N, k)/factorial(k);

}
