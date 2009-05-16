/***************************************************************************
 * K_G_J_spectral.h  Components of DGF in the spectral domain
 *                   Based upon Michalski paper, AP-IEEE 1997 
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
using namespace blitz;

complex<double> KA_xx_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> KA_zx_zy_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> KA_xz_yz_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> KA_zz_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> K_phi_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> grad_K_phi_x_y_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> grad_K_phi_z_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> grad_prime_K_phi_x_y_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> grad_prime_K_phi_z_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> G_HJ_xx_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> G_HJ_xy_yx_1_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> G_HJ_xz_yz_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> G_HJ_zx_zy_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> G_EJ_xx_1_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> G_EJ_xx_2_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> G_EJ_xz_yz_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> G_EJ_zx_zy_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> G_EJ_zz_spectral (const complex<double> k_rho, const double v, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);
