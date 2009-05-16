/***************************************************************************
 * K_G_J_spectral.cpp  Components of DGF in the spectral domain
 *                     Based upon Michalski paper, AP-IEEE 1997 
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

const complex<double> I (0.0, 1.0);

#include "besselj_complex.h"
#include "layers_constants.h"
#include "TLGF.h"

/* dyadic elements related to J */

void k_z_i_Z_i_h_Z_i_e (Array<complex<double>, 1>& k_z_i, Array<complex<double>, 1>& Z_i_h, Array<complex<double>, 1>& Z_i_e, const complex<double> k_rho, const layers_constants & LC) {

  Array<complex<double>, 1> k_z_i_tmp (LC.N);
  k_z_i_tmp = sqrt(LC.k_i*LC.k_i - k_rho*k_rho);
  k_z_i = real(k_z_i_tmp)-I*abs(imag(k_z_i_tmp));
  Z_i_h = LC.w*LC.mu_0*LC.mu_i/k_z_i;
  Z_i_e = k_z_i/(LC.w*LC.eps_0*LC.eps_i);
}

complex<double> KA_xx_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> V_i_h;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  V_i_h = V_i (Z_i_h, k_z_i, z, z_prime, m, n, LC);

  return -I/LC.w * V_i_h * besselj_complex(k_rho * rho, v, 1e-15) * k_rho;
}

complex<double> KA_zx_zy_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> I_i_h, I_i_e;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  I_i_h = I_i (Z_i_h, k_z_i, z, z_prime, m, n, LC);
  I_i_e = I_i (Z_i_e, k_z_i, z, z_prime, m, n, LC);

  return -LC.mu_0 * LC.mu_i(m) * (I_i_h - I_i_e) * besselj_complex(k_rho * rho, v, 1e-15);
}

complex<double> KA_xz_yz_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> V_v_h, V_v_e;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  V_v_e = V_v (Z_i_e, k_z_i, z, z_prime, m, n, LC);
  V_v_h = V_v (Z_i_h, k_z_i, z, z_prime, m, n, LC);

  return -LC.mu_0 * LC.mu_i(n) * (V_v_h - V_v_e) * besselj_complex(k_rho * rho, v, 1e-15);
}

complex<double> KA_zz_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> I_v_h, I_v_e;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  I_v_e = I_v (Z_i_e, k_z_i, z, z_prime, m, n, LC);
  I_v_h = I_v (Z_i_h, k_z_i, z, z_prime, m, n, LC);

  return -I * LC.mu_0 * ( LC.mu_i(m)/(LC.w*LC.eps_0*LC.eps_i(n)) * I_v_e + LC.mu_i(n) * (LC.w*LC.mu_0*LC.mu_i(m)*I_v_h - k_z_i(m)*k_z_i(m)/(LC.w*LC.eps_0*LC.eps_i(m)) * I_v_e) / (k_rho*k_rho) ) * besselj_complex(k_rho * rho, v, 1e-15) * k_rho;
}

complex<double> K_phi_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> V_i_h, V_i_e;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  V_i_e = V_i (Z_i_e, k_z_i, z, z_prime, m, n, LC);
  V_i_h = V_i (Z_i_h, k_z_i, z, z_prime, m, n, LC);

  return I * LC.w * (V_i_e - V_i_h) * besselj_complex(k_rho * rho, v, 1e-15)/k_rho;
}

complex<double> grad_K_phi_x_y_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> V_i_h, V_i_e;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  V_i_e = V_i (Z_i_e, k_z_i, z, z_prime, m, n, LC);
  V_i_h = V_i (Z_i_h, k_z_i, z, z_prime, m, n, LC);

  return I * LC.w * (V_i_h - V_i_e) * besselj_complex(k_rho * rho, v, 1e-15);
}

complex<double> grad_K_phi_z_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> I_i_h, I_i_e;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  I_i_h = I_i (Z_i_h, k_z_i, z, z_prime, m, n, LC);
  I_i_e = I_i (Z_i_e, k_z_i, z, z_prime, m, n, LC);

  return -LC.w * (LC.w*LC.mu_0*LC.mu_i(m) * I_i_h - k_z_i(m)*k_z_i(m)/(LC.w*LC.eps_0*LC.eps_i(m)) * I_i_e)/k_rho * besselj_complex(k_rho * rho, v, 1e-15);
}

complex<double> grad_prime_K_phi_x_y_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {
  
  return -grad_K_phi_x_y_spectral (k_rho, v, rho, z, z_prime, m, n, LC);
}

complex<double> grad_prime_K_phi_z_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> V_v_h, V_v_e;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  V_v_e = V_v (Z_i_e, k_z_i, z, z_prime, m, n, LC);
  V_v_h = V_v (Z_i_h, k_z_i, z, z_prime, m, n, LC);

  return LC.w * (LC.w*LC.mu_0*LC.mu_i(n) * V_v_h - k_z_i(n)*k_z_i(n)/(LC.w*LC.eps_0*LC.eps_i(n)) * V_v_e)/k_rho * besselj_complex(k_rho * rho, v, 1e-15);
}

/* G_HJ */

complex<double> G_HJ_xx_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> I_i_e, I_i_h;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  I_i_h = I_i (Z_i_h, k_z_i, z, z_prime, m, n, LC);
  I_i_e = I_i (Z_i_e, k_z_i, z, z_prime, m, n, LC);

  return (I_i_e - I_i_h) * k_rho * besselj_complex(k_rho * rho, v, 1e-15);
}

complex<double> G_HJ_xy_yx_1_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> I_i_e, I_i_h;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  I_i_h = I_i (Z_i_h, k_z_i, z, z_prime, m, n, LC);
  I_i_e = I_i (Z_i_e, k_z_i, z, z_prime, m, n, LC);

  return (I_i_e + I_i_h) * k_rho * besselj_complex(k_rho * rho, v, 1e-15);
}

complex<double> G_HJ_xz_yz_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  // j/(LC.w*LC.eps_0*LC.eps_i(n)) * I_v_e * k_rho*k_rho * besselj_complex(k_rho * rho, v, 1e-15)
  complex<double> I_v_e;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  I_v_e = I_v (Z_i_e, k_z_i, z, z_prime, m, n, LC);

  return I/(LC.w*LC.eps_0*LC.eps_i(n)) * I_v_e * k_rho*k_rho * besselj_complex(k_rho * rho, v, 1e-15);
}

complex<double> G_HJ_zx_zy_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  // j/(LC.w*LC.mu_0*LC.mu_i(m)) * V_i_h * k_rho*k_rho * besselj_complex(k_rho * rho, v, 1e-15)
  complex<double> V_i_h;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  V_i_h = V_i (Z_i_h, k_z_i, z, z_prime, m, n, LC);

  return I/(LC.w*LC.mu_0*LC.mu_i(m)) * V_i_h * k_rho*k_rho * besselj_complex(k_rho * rho, v, 1e-15);
}

/* G_EJ */

complex<double> G_EJ_xx_1_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> V_i_e, V_i_h;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  V_i_e = V_i (Z_i_e, k_z_i, z, z_prime, m, n, LC);
  V_i_h = V_i (Z_i_h, k_z_i, z, z_prime, m, n, LC);

  return (V_i_e + V_i_h) * k_rho * besselj_complex(k_rho * rho, v, 1e-15);
}

complex<double> G_EJ_xx_2_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> V_i_e, V_i_h;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  V_i_e = V_i (Z_i_e, k_z_i, z, z_prime, m, n, LC);
  V_i_h = V_i (Z_i_h, k_z_i, z, z_prime, m, n, LC);

  return (V_i_e - V_i_h) * k_rho * besselj_complex(k_rho * rho, v, 1e-15);
}

complex<double> G_EJ_xz_yz_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> V_v_e;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  V_v_e = V_v (Z_i_e, k_z_i, z, z_prime, m, n, LC);

  return -I/(LC.w*LC.eps_0*LC.eps_i(n)) * V_v_e * k_rho*k_rho * besselj_complex(k_rho * rho, v, 1e-15);
}

complex<double> G_EJ_zx_zy_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> I_i_e;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  I_i_e = I_i (Z_i_e, k_z_i, z, z_prime, m, n, LC);

  return -I/(LC.w*LC.eps_0*LC.eps_i(m)) * I_i_e * k_rho*k_rho * besselj_complex(k_rho * rho, v, 1e-15);
}

complex<double> G_EJ_zz_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> I_v_e;
  Array<complex<double>, 1> k_z_i (LC.N), Z_i_h (LC.N), Z_i_e (LC.N);
  k_z_i_Z_i_h_Z_i_e (k_z_i, Z_i_h, Z_i_e, k_rho, LC);

  I_v_e = I_v (Z_i_e, k_z_i, z, z_prime, m, n, LC);

  return -1.0/(LC.w*LC.w*LC.eps_0*LC.eps_0*LC.eps_i(m)*LC.eps_i(n)) * I_v_e * k_rho*k_rho*k_rho * besselj_complex(k_rho * rho, v, 1e-15);
}
