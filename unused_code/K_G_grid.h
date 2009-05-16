/***************************************************************************
 * K_G_grid.h  structures containing the grids for the DGFs
 *         
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

/* dyadics related to J */
typedef struct {
  Array<complex<double>, 3> xx;
  Array<complex<double>, 3> xz_yz;
  Array<complex<double>, 3> zx_zy;
  Array<complex<double>, 3> zz;
  Array<double, 1> rho_values;
  Array<double, 1> z_values;
  Array<double, 1> z_prime_values;
} KA_grid;

typedef struct {
  Array<complex<double>, 3> K;
  Array<complex<double>, 3> grad_x;
  Array<complex<double>, 3> grad_z;
  Array<complex<double>, 3> grad_prime_x;
  Array<complex<double>, 3> grad_prime_z;
  Array<double, 1> rho_values;
  Array<double, 1> z_values;
  Array<double, 1> z_prime_values;
} Kphi_grid;

typedef struct {
  Array<complex<double>, 3> xx; // -1/2 * S_2 {I_i_e - I_i_h}
  Array<complex<double>, 3> xy_1; // 1/2 * S_0 {I_i_h + I_i_e}
  Array<complex<double>, 3> xz_yz; // j/(w*eps_0*LC.eps(n)) * S_1 {k_rho * I_v_e}
  Array<complex<double>, 3> zx_zy; // j/(w*mu_0*LC.mu(m)) * S_1 {k_rho * V_i_h}
  Array<double, 1> rho_values;
  Array<double, 1> z_values;
  Array<double, 1> z_prime_values;
} G_HJ_grid;

typedef struct {
  Array<complex<double>, 3> xx_1; // -1/2 * S_0 {V_i_e + V_i_h}
  Array<complex<double>, 3> xx_2; // 1/2 * S_2 {V_i_e - V_i_h}
  Array<complex<double>, 3> xz_yz; // -j/(w*eps_0*LC.eps(n)) * S_1 {k_rho * V_v_e}
  Array<complex<double>, 3> zx_zy; // -j/(w*eps_0*LC.eps(m)) * S_1 {k_rho * I_i_e}
  Array<complex<double>, 3> zz; // -1/(w^2*eps_0^2*LC.eps(n)*LC.eps(m)) * S_0 {k_rho^2 * I_v_e}
  Array<double, 1> rho_values;
  Array<double, 1> z_values;
  Array<double, 1> z_prime_values;
} G_EJ_grid;

/* dyadics related to M */
typedef struct {
  Array<complex<double>, 3> xx;
  Array<complex<double>, 3> xz_yz;
  Array<complex<double>, 3> zx_zy;
  Array<complex<double>, 3> zz;
  Array<double, 1> rho_values;
  Array<double, 1> z_values;
  Array<double, 1> z_prime_values;
} KF_grid;

typedef struct {
  Array<complex<double>, 3> K;
  Array<complex<double>, 3> grad_x;
  Array<complex<double>, 3> grad_z;
  Array<complex<double>, 3> grad_prime_x;
  Array<complex<double>, 3> grad_prime_z;
  Array<double, 1> rho_values;
  Array<double, 1> z_values;
  Array<double, 1> z_prime_values;
} Kpsi_grid;

typedef struct {
  Array<complex<double>, 3> xx; 
  Array<complex<double>, 3> xy_1; 
  Array<complex<double>, 3> xz_yz;
  Array<complex<double>, 3> zx_zy;
  Array<double, 1> rho_values;
  Array<double, 1> z_values;
  Array<double, 1> z_prime_values;
} G_EM_grid;

typedef struct {
  Array<complex<double>, 3> xx_1;
  Array<complex<double>, 3> xx_2;
  Array<complex<double>, 3> xz_yz;
  Array<complex<double>, 3> zx_zy;
  Array<complex<double>, 3> zz;
  Array<double, 1> rho_values;
  Array<double, 1> z_values;
  Array<double, 1> z_prime_values;
} G_HM_grid;
