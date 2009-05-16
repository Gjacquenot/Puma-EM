/***************************************************************************
 * interp_K_G_grid.cpp  routine for DGF tables interpolation
 *                      the grids are regularly spaced
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

#include "K_G_grid.h"

int index_search (const Array<double, 1> & x0, const double x, int & ind_x, int & ind_x_sup, double & mu_x) {
  // only valid for equispaced position vectors
  int N_x = x0.size();
  if (x>x0(N_x-1)) {
    cout << "out of interpolation range, too big: x > x0 (end): " << x << " > " << x0(N_x-1) << endl;
    exit (1);
  }
  if (x<x0(0)) {
    cout << "out of interpolation range, too small: x < x0 (0): " << x << " < " << x0(0) << endl;
    exit (1);
  }
  if (N_x>1) {
    double Delta_x = x0(1) - x0(0);
    ind_x = (int) floor((x-x0(0))/Delta_x);
    if (ind_x+1==N_x) ind_x = ind_x-1;
    ind_x_sup = ind_x + 1;
    mu_x = (x-x0(ind_x))/Delta_x;
  }
  else {
    ind_x = 0;
    ind_x_sup = 0;
    mu_x = 0.0;
  }
  return 0;
}

complex<double> interp_K_tab (const Array<complex<double>, 3>& K_tab, const int ind_z, const int ind_z_sup, const int ind_z_prime, const int ind_z_prime_sup, const int ind_rho, const int ind_rho_sup, const double mu_z, const double mu_z_prime, const double mu_rho) {

  double eta;

  complex<double> K_000, K_100, K_010, K_110, K_001, K_101, K_011, K_111;
  complex<double> K_00z, K_01z, K_10z, K_11z, K_0yz, K_1yz, K_xyz;
  
  K_000 = K_tab (ind_z, ind_z_prime, ind_rho);
  K_001 = K_tab (ind_z, ind_z_prime, ind_rho_sup);
  K_010 = K_tab (ind_z, ind_z_prime_sup, ind_rho);
  K_011 = K_tab (ind_z, ind_z_prime_sup, ind_rho_sup);
  K_100 = K_tab (ind_z_sup, ind_z_prime, ind_rho);
  K_101 = K_tab (ind_z_sup, ind_z_prime, ind_rho_sup);
  K_110 = K_tab (ind_z_sup, ind_z_prime_sup, ind_rho);
  K_111 = K_tab (ind_z_sup, ind_z_prime_sup, ind_rho_sup);

  // z = x, z_prime = y, rho = z

  eta = 1-mu_rho;
  K_00z = K_000*eta + K_001*mu_rho;
  K_01z = K_010*eta + K_011*mu_rho;
  K_10z = K_100*eta + K_101*mu_rho;
  K_11z = K_110*eta + K_111*mu_rho;

  eta = 1.0-mu_z_prime;
  K_0yz = K_00z*eta + K_01z*mu_z_prime;
  K_1yz = K_10z*eta + K_11z*mu_z_prime;
  
  eta = 1.0-mu_z;
  K_xyz = K_0yz*eta + K_1yz*mu_z;

  return K_xyz;
}

void interp_KA_mm (Array<complex<double>, 2>& KA, const KA_grid & KA_tab, const TinyVector<double,3>& r, const TinyVector<double,3>& r_prime) {

  // this is a 3D interpolation scheme; we have to treat the case
  // where the array is a degenerated 2D-array

  double rho = sqrt( pow2 (r(0)-r_prime(0)) + pow2 (r(1)-r_prime(1)) );
  double z = r(2), z_prime = r_prime(2), mu_z, mu_z_prime, mu_rho;

  int ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup;
  index_search (KA_tab.z_values, z, ind_z, ind_z_sup, mu_z);
  index_search (KA_tab.z_prime_values, z_prime, ind_z_prime, ind_z_prime_sup, mu_z_prime);
  index_search (KA_tab.rho_values, rho, ind_rho, ind_rho_sup, mu_rho);

  KA (1, 0) = 0.0;
  KA (0, 1) = 0.0;
  KA (0, 0) = interp_K_tab (KA_tab.xx, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
  KA (1, 1) = KA (0, 0);
  if (rho>0.0) {
    double cos_phi = (r(0)-r_prime(0))/rho, sin_phi = (r(1)-r_prime(1))/rho;
    complex<double> KA_zx_zy, KA_xz_yz;
    KA_zx_zy = interp_K_tab (KA_tab.zx_zy, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
    KA_xz_yz = interp_K_tab (KA_tab.xz_yz, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
    KA (0,2) = cos_phi * KA_xz_yz;
    KA (1,2) = sin_phi * KA_xz_yz;
    KA (2,0) = cos_phi * KA_zx_zy;
    KA (2,1) = sin_phi * KA_zx_zy;
  }
  else {
    KA (0,2) = 0.0;
    KA (1,2) = 0.0;
    KA (2,0) = 0.0;
    KA (2,1) = 0.0;    
  }
  KA (2,2) = interp_K_tab (KA_tab.zz, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
}

void interp_Kphi_mm (complex<double>& K_phi, const Kphi_grid & Kphi_tab, const TinyVector<double,3>& r, const TinyVector<double,3>& r_prime) {

  // this is a 3D interpolation scheme; we have to treat the case
  // where the array is a degenerated 2D-array

  double rho = sqrt( pow2 (r(0)-r_prime(0)) + pow2 (r(1)-r_prime(1)) );
  double z = r(2), z_prime = r_prime(2), mu_z, mu_z_prime, mu_rho;

  int ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup;
  index_search (Kphi_tab.z_values, z, ind_z, ind_z_sup, mu_z);
  index_search (Kphi_tab.z_prime_values, z_prime, ind_z_prime, ind_z_prime_sup, mu_z_prime);
  index_search (Kphi_tab.rho_values, rho, ind_rho, ind_rho_sup, mu_rho);

  K_phi = interp_K_tab (Kphi_tab.K, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
}

void interp_grad_K_phi (TinyVector<complex<double>, 3>& grad_K_phi, const Kphi_grid & Kphi_tab, const TinyVector<double,3>& r, const TinyVector<double,3>& r_prime) {

  // this is a 3D interpolation scheme; we have to treat the case
  // where the array is a degenerated 2D-array

  double rho = sqrt( pow2 (r(0)-r_prime(0)) + pow2 (r(1)-r_prime(1)) );
  double z = r(2), z_prime = r_prime(2), mu_z, mu_z_prime, mu_rho;

  int ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup;
  index_search (Kphi_tab.z_values, z, ind_z, ind_z_sup, mu_z);
  index_search (Kphi_tab.z_prime_values, z_prime, ind_z_prime, ind_z_prime_sup, mu_z_prime);
  index_search (Kphi_tab.rho_values, rho, ind_rho, ind_rho_sup, mu_rho);
  if (rho>0.0) {
    double cos_phi = (r(0)-r_prime(0))/rho, sin_phi = (r(1)-r_prime(1))/rho;
    complex<double> grad_K_phi_xy;
    grad_K_phi_xy = interp_K_tab (Kphi_tab.grad_x, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
    grad_K_phi (0) = cos_phi * grad_K_phi_xy;
    grad_K_phi (1) = sin_phi * grad_K_phi_xy;
  }
  else {
    grad_K_phi (0) = 0.0;
    grad_K_phi (1) = 0.0;
  }
  grad_K_phi (2) = interp_K_tab (Kphi_tab.grad_z, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
}

void interp_grad_prime_K_phi (TinyVector<complex<double>, 3>& grad_prime_K_phi, const Kphi_grid & Kphi_tab, const TinyVector<double,3>& r, const TinyVector<double,3>& r_prime) {

  // this is a 3D interpolation scheme; we have to treat the case
  // where the array is a degenerated 2D-array

  double rho = sqrt( pow2 (r(0)-r_prime(0)) + pow2 (r(1)-r_prime(1)) );
  double z = r(2), z_prime = r_prime(2), mu_z, mu_z_prime, mu_rho;

  int ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup;
  index_search (Kphi_tab.z_values, z, ind_z, ind_z_sup, mu_z);
  index_search (Kphi_tab.z_prime_values, z_prime, ind_z_prime, ind_z_prime_sup, mu_z_prime);
  index_search (Kphi_tab.rho_values, rho, ind_rho, ind_rho_sup, mu_rho);
  if (rho>0.0) {
    double cos_phi = (r(0)-r_prime(0))/rho, sin_phi = (r(1)-r_prime(1))/rho;
    complex<double> grad_prime_K_phi_xy;
    grad_prime_K_phi_xy = interp_K_tab (Kphi_tab.grad_prime_x, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
    grad_prime_K_phi (0) = cos_phi * grad_prime_K_phi_xy;
    grad_prime_K_phi (1) = sin_phi * grad_prime_K_phi_xy;
  }
  else {
    grad_prime_K_phi (0) = 0.0;
    grad_prime_K_phi (1) = 0.0;
  }
  grad_prime_K_phi (2) = interp_K_tab (Kphi_tab.grad_prime_z, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
}

void interp_G_HJ_mm (Array<complex<double>, 2>& G_HJ, const G_HJ_grid & G_HJ_tab, const TinyVector<double,3>& r, const TinyVector<double,3>& r_prime) {

  // this is a 3D interpolation scheme; we have to treat the case
  // where the array is a degenerated 2D-array

  double rho = sqrt( pow2 (r(0)-r_prime(0)) + pow2 (r(1)-r_prime(1)) );
  double z = r(2), z_prime = r_prime(2), mu_z, mu_z_prime, mu_rho;

  int ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup;
  index_search (G_HJ_tab.z_values, z, ind_z, ind_z_sup, mu_z);
  index_search (G_HJ_tab.z_prime_values, z_prime, ind_z_prime, ind_z_prime_sup, mu_z_prime);
  index_search (G_HJ_tab.rho_values, rho, ind_rho, ind_rho_sup, mu_rho);

  complex<double> G_HJ_xy_1 = interp_K_tab (G_HJ_tab.xy_1, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
  G_HJ (0, 1) = 0.5 * G_HJ_xy_1;
  G_HJ (1, 0) = -G_HJ (0, 1); // = -0.5 * G_HJ_xy_1
  if (rho>0.0) {
    double cos_phi = (r(0)-r_prime(0))/rho, sin_phi = (r(1)-r_prime(1))/rho, cos_2_phi = cos_phi*cos_phi - sin_phi*sin_phi, sin_2_phi = 2.0*sin_phi*cos_phi;
    complex<double> G_HJ_xx_tmp, G_HJ_zx_zy, G_HJ_xz_yz;
    G_HJ_xx_tmp = interp_K_tab (G_HJ_tab.xx, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
    G_HJ (0, 0) = -0.5 * sin_2_phi * G_HJ_xx_tmp;
    G_HJ (1, 1) = -G_HJ (0, 0);

    G_HJ (0, 1) += 0.5 * cos_2_phi * G_HJ_xx_tmp;
    G_HJ (1, 0) += 0.5 * cos_2_phi * G_HJ_xx_tmp;

    G_HJ_zx_zy = interp_K_tab (G_HJ_tab.zx_zy, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
    G_HJ_xz_yz = interp_K_tab (G_HJ_tab.xz_yz, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
    G_HJ (0, 2) = sin_phi * G_HJ_xz_yz;
    G_HJ (1, 2) = -cos_phi * G_HJ_xz_yz;
    G_HJ (2, 0) = -sin_phi * G_HJ_zx_zy;
    G_HJ (2, 1) = cos_phi * G_HJ_zx_zy;
  }
  else {
    G_HJ (0, 0) = 0.0;
    G_HJ (1, 1) = 0.0;

    G_HJ (0, 2) = 0.0;
    G_HJ (1, 2) = 0.0;
    G_HJ (2, 0) = 0.0;
    G_HJ (2, 1) = 0.0;    
  }
  G_HJ (2,2) = 0.0;
}

void interp_G_EJ_mm (Array<complex<double>, 2>& G_EJ, const G_EJ_grid & G_EJ_tab, const TinyVector<double,3>& r, const TinyVector<double,3>& r_prime) {

  // this is a 3D interpolation scheme; we have to treat the case
  // where the array is a degenerated 2D-array

  double rho = sqrt( pow2 (r(0)-r_prime(0)) + pow2 (r(1)-r_prime(1)) );
  double z = r(2), z_prime = r_prime(2), mu_z, mu_z_prime, mu_rho;

  int ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup;
  index_search (G_EJ_tab.z_values, z, ind_z, ind_z_sup, mu_z);
  index_search (G_EJ_tab.z_prime_values, z_prime, ind_z_prime, ind_z_prime_sup, mu_z_prime);
  index_search (G_EJ_tab.rho_values, rho, ind_rho, ind_rho_sup, mu_rho);

  G_EJ (0, 0) = -0.5 * interp_K_tab (G_EJ_tab.xx_1, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
  G_EJ (1, 1) = G_EJ (0, 0);
  if (rho>0.0) {
    double cos_phi = (r(0)-r_prime(0))/rho, sin_phi = (r(1)-r_prime(1))/rho, cos_2_phi = cos_phi*cos_phi - sin_phi*sin_phi, sin_2_phi = 2.0*sin_phi*cos_phi;

    complex<double> G_EJ_xx_2 = interp_K_tab (G_EJ_tab.xx_2, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
    G_EJ (0, 0) += 0.5 * cos_2_phi *  G_EJ_xx_2;
    G_EJ (1, 1) -= 0.5 * cos_2_phi *  G_EJ_xx_2;

    G_EJ (1, 0) = 0.5 * sin_2_phi * G_EJ_xx_2;
    G_EJ (0, 1) = G_EJ (1, 0);

    complex<double> G_EJ_zx_zy, G_EJ_xz_yz;
    G_EJ_zx_zy = interp_K_tab (G_EJ_tab.zx_zy, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
    G_EJ_xz_yz = interp_K_tab (G_EJ_tab.xz_yz, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
    G_EJ (0, 2) = cos_phi * G_EJ_xz_yz;
    G_EJ (1, 2) = sin_phi * G_EJ_xz_yz;
    G_EJ (2, 0) = cos_phi * G_EJ_zx_zy;
    G_EJ (2, 1) = sin_phi * G_EJ_zx_zy;
  }
  else {
    G_EJ (1, 0) = 0.0;
    G_EJ (0, 1) = 0.0;

    G_EJ (0, 2) = 0.0;
    G_EJ (1, 2) = 0.0;
    G_EJ (2, 0) = 0.0;
    G_EJ (2, 1) = 0.0;    
  }
  G_EJ (2,2) = interp_K_tab (G_EJ_tab.zz, ind_z, ind_z_sup, ind_z_prime, ind_z_prime_sup, ind_rho, ind_rho_sup, mu_z, mu_z_prime, mu_rho);
}
