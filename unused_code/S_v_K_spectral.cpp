/***************************************************************************
 * S_v_K_spectral.cpp  Sommerfeld integration of spectral domain DGFs
 *                     Based upon Michalski paper, AP-IEEE 1997-1998 
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
#include <float.h>
#include <complex>
#include <blitz/array.h>

using namespace blitz;

const complex<double> I (0.0, 1.0);

#include "layers_constants.h"
#include "GK.h"
#include "weighted_average.h"

complex<double> int_ponc_K_spectral (complex<double> (*f)(const complex<double>, const double, const double, const double, const double, const int, const int, const layers_constants &), const double v, const complex<double> a, const complex<double> b, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC, const int nb_intervals) {

  int j;
  firstIndex i;
  complex<double> k_rho_deriv;
  Array<complex<double>, 1> k_rho(nb_intervals+1);
  complex<double> K (0.0, 0.0);
  k_rho = a + (i+0.0) * ((b-a)/(double) nb_intervals);
  k_rho_deriv = (b-a)/abs(b-a);
  for (j=0 ; j<nb_intervals ; j++) K += GK31 (f, v, k_rho (j), k_rho (j+1), rho, z, z_prime, m, n, LC);
  return 0.5/M_PI * K * k_rho_deriv;
}

complex<double> S_v_K_spectral (complex<double> (*f)(const complex<double>, const double, const double, const double, const double, const int, const int, const layers_constants &), const double v, const double rho, const double z, const double z_prime, const layers_constants & LC) {

  // we must find in which layer are z and z_prime
  int n = count(LC.z_i<z_prime), m = count(LC.z_i<z);

  int nb_zeros = 25, j=0;
  firstIndex i;
  double e, e_max, q, v0 = 0.0, first_zero, k_min;
  complex<double> K, int_ellipse_K;
  Array<double, 1> zeros_int(Range(1, nb_zeros));

  //-------------------------------------------------------------------
  // integration along the deformed  path. We must find the value of e,
  // which is the minimum value for the length of the real part of the 
  // deformed path (a rectangle). In general, e = k_0*(max(sqrt(eps_i))+1).
  //-------------------------------------------------------------------
  e = LC.k_0+max(real(LC.k_i));
  e_max = 500;
  if (e>e_max) e = e_max;

  //---------------------------------------------------------------------
  // for a Bessel function of order v, J_v(x), the zeros are given by the 
  // equation : x - v*pi/2 - pi_4 = (2*k+1)*pi/2, where k is a natural.
  // Hereafter we try to find the first zero after the deformed path,
  // in order to integrate up to there and after return to the real axis
  // in order to proceed the sommerfeld tail through extrapolation
  //---------------------------------------------------------------------
  if (rho>0) {
    q = M_PI/rho;
    k_min = ceil(e*rho/M_PI - 0.25 - v/2.0 - 0.5);
    first_zero = ( (2*k_min+1)*M_PI/2 + M_PI/4 + v*M_PI/2 )/rho;
    zeros_int = first_zero + (i-1)*q;
  }
  else {
    q = 1e2;
    zeros_int = e + (i-1)*q;
  }

  //---------------------------------------------------------
  // integration over the deformed path (a rectangle).
  // We have k_rho = a + (b-a)*x, where a and/or b may be 
  // complex, but x is real. 
  //---------------------------------------------------------
  int nb_interv_1 = 1, nb_interv_2 = 5;
  complex<double> a1 (0.0, 1e-8), b1 (0.0, 10.0);
  complex<double> a2 (0.0, 10.0), b2 (zeros_int(1), 10.0);
  complex<double> a3 (zeros_int(1), 10.0), b3 (zeros_int(1), 0.0);
  int_ellipse_K = int_ponc_K_spectral (f, v, a1, b1, rho, z, z_prime, m, n, LC, nb_interv_1) + int_ponc_K_spectral (f, v, a2, b2, rho, z, z_prime, m, n, LC, nb_interv_2) + int_ponc_K_spectral (f, v, a3, b3, rho, z, z_prime, m, n, LC, nb_interv_1);

  //---------------------------------------------------------
  // integration over the real axis. We calculate the terms
  // as long as the partial integrals are big enough compared 
  // to the deformed integral. The serie is after accelerated 
  // using the extrapolation method developed in michalski 98
  //---------------------------------------------------------
  j = 0;
  int cond = (abs(int_ellipse_K) == 0.0) ? 0 : 1, cond_1, cond_2;
  Array<complex<double>, 1> u_i(Range(1, nb_zeros-1));
  while (cond) {
    j++;
    u_i(j) = int_ponc_K_spectral (f, v, zeros_int(j), zeros_int(j+1), rho, z, z_prime, m, n, LC, nb_interv_1*3);
    cond_1 = abs(real(u_i(j))) > 1e-8*abs(real(int_ellipse_K));
    cond_2 = abs(imag(u_i(j))) > 1e-8*abs(imag(int_ellipse_K));
    cond = ( (cond_1 || cond_2) && (j < nb_zeros-1) );
  }

  // extrapolation part
  int R = j;
  K = (R>=1) ? int_ellipse_K + weighted_average (R, zeros_int, u_i) : int_ellipse_K;
  return K;
}

