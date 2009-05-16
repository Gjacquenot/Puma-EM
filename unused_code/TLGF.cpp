/***************************************************************************
 * TLGF.cpp  Transmission Line Green's Functions
 *           Based upon Michalski paper, AP-IEEE 1997 
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
#include <blitz/tinyvec-et.h>

using namespace blitz;

const complex<double> I (0.0, 1.0);

#include "layers_constants.h"

void t_G_left_right (Array<complex<double>, 1>& t, Array<complex<double>, 1>& G_left, Array<complex<double>, 1>& G_right, const Array<complex<double>, 1>& Z_i, const Array<complex<double>, 1>& k_z_i, const int m, const int n, const layers_constants & LC) {

  int N = LC.N, l;
  G_left (0) = 0.0;
  G_right (N-1) = 0.0;
  complex<double> G;

  // computation of t; t(0) = t(end) = 0.
  t = exp (-2.0 * I * k_z_i * LC.d_i);
  t (0) = 0.0;
  t (N-1) = 0.0;

  // computation of G_left; G_left(0) = 0.
  for (l=1 ; l<=n ; l++) {
    G = (Z_i (l-1) - Z_i (l))/(Z_i (l-1) + Z_i (l));
    G_left (l) = (G + G_left (l-1)*t (l-1))/(1.0 + G * G_left (l-1) * t (l-1));
  }

  // computation of G_right; G_right(N-1) = 0.
  for (l=N-2 ; l>=m ; l--) {
    G = (Z_i (l+1) - Z_i (l))/(Z_i (l+1) + Z_i (l));
    G_right (l) = (G + G_right (l+1) * t (l+1))/(1.0 + G * G_right (l+1) * t (l+1));
  }

}

TinyVector<complex<double>, 4> gamma_n (const int n, const double z, const double z_prime, const layers_constants & LC) {

  double z_plus_z_prime = z+z_prime, z_minus_z_prime = z-z_prime;
  TinyVector<complex<double>, 4> g;
  if (n==LC.N-1) g = 0.0, z_plus_z_prime - 2.0*LC.z_i(n-1), 0.0, 0.0;
  else if (n>0) g = 2.0 * LC.z_i (n) - z_plus_z_prime, z_plus_z_prime - 2.0 * LC.z_i (n-1), 2.0 * LC.d_i (n) + z_minus_z_prime, 2.0 * LC.d_i (n) - z_minus_z_prime;
  else g = 2.0*LC.z_i (n) - z_plus_z_prime, 0.0, 0.0, 0.0;
  return g;
}

complex<double> T_left (const int m, const int n, const Array<complex<double>, 1>& t, const Array<complex<double>, 1>& G_left, const Array<complex<double>, 1>& k_z_i, const layers_constants & LC) {

  int q;
  complex<double> T = 1.0;
  if ((m+1)<=(n-1)) {
    for (q=m+1 ; q<=n-1 ; q++) T *= (1.0+G_left(q))/(1.0+G_left(q)*t(q)) * exp(-I * k_z_i (q) * LC.d_i (q));
  }
  return T;
}

complex<double> V_i_mn (const Array<complex<double>, 1>& Z_i, const Array<complex<double>, 1>& k_z_i, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> V, D_n;
  TinyVector<complex<double>, 4> R_n;
  Array<complex<double>, 1> t (LC.N), G_left (LC.N), G_right (LC.N);

  t_G_left_right (t, G_left, G_right, Z_i, k_z_i, m, n, LC);

  D_n = 1.0 - G_left(n)*G_right(n)*t(n);
  if (m==n) {
    R_n = G_right(n), G_left(n), G_left(n)*G_right(n), G_left(n)*G_right(n);
    V = 0.5*Z_i(n)/D_n * sum (R_n * exp (-I * k_z_i(n) * gamma_n (n, z, z_prime, LC)));
  }
  else {
    R_n = G_right(n), G_left(n), G_left(n)*G_right(n), G_left(n)*G_right(n);
    V = 0.5 * Z_i (n) * ( exp(-I*k_z_i(n)*abs(LC.z_i(n-1)-z_prime)) +  1.0/D_n * sum (R_n * exp (-I * k_z_i (n) * gamma_n (n, LC.z_i (n-1), z_prime, LC))) ) * T_left (m, n, t, G_left, k_z_i, LC) * exp(-I*k_z_i(m)*(LC.z_i(m)-z));
    if (m>0) V *= (1.0+G_left(m)*exp(-2.0*I*k_z_i(m)*(z-LC.z_i(m-1))))/(1.0+G_left(m)*t(m));
  }
  return V;
}

complex<double> V_v_mn (const Array<complex<double>, 1>& Z_i, const Array<complex<double>, 1>& k_z_i, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> V, D_n;
  TinyVector<complex<double>, 4> R_n;
  Array<complex<double>, 1> t (LC.N), G_left (LC.N), G_right (LC.N);

  t_G_left_right (t, G_left, G_right, Z_i, k_z_i, m, n, LC);

  D_n = 1.0 - G_left(n)*G_right(n)*t(n);
  if (m==n) {
    R_n = G_right(n), -G_left(n), G_left(n)*G_right(n), -G_left(n)*G_right(n);
    V = 0.5/D_n * sum (R_n * exp (-I * k_z_i(n) * gamma_n (n, z, z_prime, LC)));      
  }
  else {
    R_n = G_right(n), -G_left(n), G_left(n)*G_right(n), -G_left(n)*G_right(n);
    V = 0.5 * ( -exp(-I*k_z_i(n)*abs(LC.z_i(n-1)-z_prime)) + 1.0/D_n * sum (R_n * exp (-I * k_z_i(n) * gamma_n (n, LC.z_i (n-1), z_prime, LC))) ) * T_left (m, n, t, G_left, k_z_i, LC) * exp(-I*k_z_i(m)*(LC.z_i(m)-z));
    if (m>0) V *= (1.0+G_left(m)*exp(-2.0*I*k_z_i(m)*(z-LC.z_i(m-1))))/(1.0+G_left(m)*t(m));
  }
  return V;
}

complex<double> I_v_mn (const Array<complex<double>, 1>& Z_i, const Array<complex<double>, 1>& k_z_i, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> I_v, D_n;
  TinyVector<complex<double>, 4> R_n;
  Array<complex<double>, 1> t (LC.N), G_left (LC.N), G_right (LC.N);

  t_G_left_right (t, G_left, G_right, Z_i, k_z_i, m, n, LC);

  D_n = 1.0 - G_left(n)*G_right(n)*t(n);
  if (m==n) {
    R_n = -G_right(n), -G_left(n), G_left(n)*G_right(n), G_left(n)*G_right(n);
    I_v = 0.5/(Z_i (n) * D_n) * sum (R_n * exp (-I * k_z_i (n) * gamma_n (n, z, z_prime, LC)));
  }
  else {
    R_n = G_right(n), -G_left(n), G_left(n)*G_right(n), -G_left(n)*G_right(n); // R_n of V_v
    I_v = 0.5 * ( -exp(-I*k_z_i(n)*abs(LC.z_i(n-1)-z_prime)) + 1.0/D_n * sum (R_n * exp (-I * k_z_i(n) * gamma_n (n, LC.z_i (n-1), z_prime, LC))) ) * T_left (m, n, t, G_left, k_z_i, LC) * exp(-I*k_z_i(m)*(LC.z_i(m)-z))/(-Z_i(m));
    if (m>0) I_v *= (1.0-G_left(m)*exp(-2.0*I*k_z_i(m)*(z-LC.z_i(m-1))))/(1.0+G_left(m)*t(m));
  }
  return I_v;
}

complex<double> I_i_mn (const Array<complex<double>, 1>& Z_i, const Array<complex<double>, 1>& k_z_i, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> I_i, D_n;
  TinyVector<complex<double>, 4> R_n;
  Array<complex<double>, 1> t (LC.N), G_left (LC.N), G_right (LC.N);

  t_G_left_right (t, G_left, G_right, Z_i, k_z_i, m, n, LC);

  D_n = 1.0 - G_left(n)*G_right(n)*t(n);
  if (m==n) {
    R_n = -G_right(n), G_left(n), G_left(n)*G_right(n), -G_left(n)*G_right(n);
    I_i = 0.5/D_n * sum (R_n * exp (-I * k_z_i (n) * gamma_n (n, z, z_prime, LC)));
  }
  else {
    R_n = G_right(n), G_left(n), G_left(n)*G_right(n), G_left(n)*G_right(n); // R_n of V_i
    I_i =  0.5 * Z_i (n) * ( exp(-I*k_z_i(n)*abs(LC.z_i(n-1)-z_prime)) +  1.0/D_n * sum (R_n * exp (-I * k_z_i (n) * gamma_n (n, LC.z_i (n-1), z_prime, LC))) ) * T_left (m, n, t, G_left, k_z_i, LC) * exp(-I*k_z_i(m)*(LC.z_i(m)-z))/(-Z_i(m));
    if (m>0) I_i *= (1.0-G_left(m)*exp(-2.0*I*k_z_i(m)*(z-LC.z_i(m-1))))/(1.0+G_left(m)*t(m));
  }
  return I_i;
}

complex<double> V_i (const Array<complex<double>, 1>& Z_i, const Array<complex<double>, 1>& k_z_i, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> result;
  if (m<=n) result = V_i_mn  (Z_i, k_z_i, z, z_prime, m, n, LC);
  else result = V_i_mn  (Z_i, k_z_i, z_prime, z, n, m, LC);
  return result;
}

complex<double> V_v (const Array<complex<double>, 1>& Z_i, const Array<complex<double>, 1>& k_z_i, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> result;
  if (m<=n) result = V_v_mn  (Z_i, k_z_i, z, z_prime, m, n, LC);
  else result = -I_i_mn  (Z_i, k_z_i, z_prime, z, n, m, LC);
  return result;
}

complex<double> I_v (const Array<complex<double>, 1>& Z_i, const Array<complex<double>, 1>& k_z_i, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> result;
  if (m<=n) result = I_v_mn  (Z_i, k_z_i, z, z_prime, m, n, LC);
  else result = I_v_mn  (Z_i, k_z_i, z_prime, z, n, m, LC);
  return result;
}

complex<double> I_i (const Array<complex<double>, 1>& Z_i, const Array<complex<double>, 1>& k_z_i, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> result;
  if (m<=n) result = I_i_mn  (Z_i, k_z_i, z, z_prime, m, n, LC);
  else result = -V_v_mn  (Z_i, k_z_i, z_prime, z, n, m, LC);
  return result;
}
