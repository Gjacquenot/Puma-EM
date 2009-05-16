/***************************************************************************
 * K_G_J.cpp  Components of DGF in the spatial domain 
 *            = Sommerfeld integration of spectral domain DGF 
 *            Based upon Michalski paper, AP-IEEE 1997 
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

using namespace blitz;

const complex<double> I (0.0, 1.0);

#include "layers_constants.h"
#include "K_G_J_spectral.h"
#include "S_v_K_spectral.h"

complex<double> KA_xx (const double rho, const double z, const double z_prime, const layers_constants & LC) {
  return S_v_K_spectral (&KA_xx_spectral, 0.0, rho, z, z_prime, LC);
}

complex<double> KA_zx_zy (const double rho, const double z, const double z_prime, const layers_constants & LC) {
  return S_v_K_spectral (&KA_zx_zy_spectral, 1.0, rho, z, z_prime, LC);
}

complex<double> KA_xz_yz (const double rho, const double z, const double z_prime, const layers_constants & LC) {
  return S_v_K_spectral (&KA_xz_yz_spectral, 1.0, rho, z, z_prime, LC);
}

complex<double> KA_zz (const double rho, const double z, const double z_prime, const layers_constants & LC) {
  return S_v_K_spectral (&KA_zz_spectral, 0.0, rho, z, z_prime, LC);
}

complex<double> K_phi (const double rho, const double z, const double z_prime, const layers_constants & LC) {
  return S_v_K_spectral (&K_phi_spectral, 0.0, rho, z, z_prime, LC);
}

complex<double> grad_K_phi_xy (const double rho, const double z, const double z_prime, const layers_constants & LC) {
  return S_v_K_spectral (&grad_K_phi_x_y_spectral, 1.0, rho, z, z_prime, LC);
}

complex<double> grad_K_phi_z (const double rho, const double z, const double z_prime, const layers_constants & LC) {
  return S_v_K_spectral (&grad_K_phi_z_spectral, 0.0, rho, z, z_prime, LC);
}

complex<double> grad_prime_K_phi_xy (const double rho, const double z, const double z_prime, const layers_constants & LC) {
  return S_v_K_spectral (&grad_prime_K_phi_x_y_spectral, 1.0, rho, z, z_prime, LC);
}

complex<double> grad_prime_K_phi_z (const double rho, const double z, const double z_prime, const layers_constants & LC) {
  return S_v_K_spectral (&grad_prime_K_phi_z_spectral, 0.0, rho, z, z_prime, LC);
}

complex<double> G_HJ_xx (const double rho, const double z, const double z_prime, const layers_constants & LC) {
  return S_v_K_spectral (&G_HJ_xx_spectral, 2.0, rho, z, z_prime, LC); // S_2 {I_i_e - I_i_h}
}

complex<double> G_HJ_xy_1 (const double rho, const double z, const double z_prime, const layers_constants & LC) {
  return S_v_K_spectral (&G_HJ_xy_yx_1_spectral, 0.0, rho, z, z_prime, LC); // S_0 {I_i_h + I_i_e}
}

complex<double> G_HJ_xz_yz (const double rho, const double z, const double z_prime, const layers_constants & LC) {
  return S_v_K_spectral (&G_HJ_xz_yz_spectral, 1.0, rho, z, z_prime, LC); // j/(w*LC.eps_0*LC.eps(n)) * S_1 {k_rho * I_v_e}
}

complex<double> G_HJ_zx_zy (const double rho, const double z, const double z_prime, const layers_constants & LC) {
  return S_v_K_spectral (&G_HJ_zx_zy_spectral, 1.0, rho, z, z_prime, LC); // j/(w*LC.mu_0*LC.mu(m)) * S_1 {k_rho * V_i_h}
}

Array<complex<double>, 1> G_EJ (const double rho, const double z, const double z_prime, const layers_constants & LC) {

  Array<complex<double>, 1> G_EJ (5);
  G_EJ (0) = S_v_K_spectral (&G_EJ_xx_1_spectral, 0.0, rho, z, z_prime, LC); // S_0 {V_i_e + V_i_h}
  G_EJ (1) = S_v_K_spectral (&G_EJ_xx_2_spectral, 2.0, rho, z, z_prime, LC); // S_2 {V_i_e - V_i_h}
  G_EJ (2) = S_v_K_spectral (&G_EJ_xz_yz_spectral, 1.0, rho, z, z_prime, LC); // -j/(w*LC.eps_0*LC.eps(n)) * S_1 {k_rho * V_v_e}
  G_EJ (3) = S_v_K_spectral (&G_EJ_zx_zy_spectral, 1.0, rho, z, z_prime, LC); // -j/(w*LC.eps_0*LC.eps(m)) * S_1 {k_rho * I_i_e}
  G_EJ (4) = S_v_K_spectral (&G_EJ_zz_spectral, 0.0, rho, z, z_prime, LC); // -1/(w^2*LC.eps_0^2*LC.eps(n)*LC.eps(m)) * S_0 {k_rho^2 * I_v_e}
  return G_EJ;
}
