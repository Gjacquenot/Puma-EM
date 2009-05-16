/***************************************************************************
 * K_G_J.h  Components of DGF in the spatial domain 
 *          = Sommerfeld integration of spectral domain DGF 
 *          Based upon Michalski paper, AP-IEEE 1997 
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

complex<double> KA_xx (const double rho, const double z, const double z_prime, const layers_constants & LC);

complex<double> KA_zx_zy (const double rho, const double z, const double z_prime, const layers_constants & LC);

complex<double> KA_xz_yz (const double rho, const double z, const double z_prime, const layers_constants & LC);

complex<double> KA_zz (const double rho, const double z, const double z_prime, const layers_constants & LC);

complex<double> K_phi (const double rho, const double z, const double z_prime, const layers_constants & LC);

complex<double> grad_K_phi_xy (const double rho, const double z, const double z_prime, const layers_constants & LC);

complex<double> grad_K_phi_z (const double rho, const double z, const double z_prime, const layers_constants & LC);

complex<double> grad_prime_K_phi_xy (const double rho, const double z, const double z_prime, const layers_constants & LC);

complex<double> grad_prime_K_phi_z (const double rho, const double z, const double z_prime, const layers_constants & LC);

complex<double> G_HJ_xx (const double rho, const double z, const double z_prime, const layers_constants & LC);

complex<double> G_HJ_xy_1 (const double rho, const double z, const double z_prime, const layers_constants & LC);

complex<double> G_HJ_xz_yz (const double rho, const double z, const double z_prime, const layers_constants & LC);

complex<double> G_HJ_zx_zy (const double rho, const double z, const double z_prime, const layers_constants & LC);

//Array<complex<double>, 1> G_HJ (const double rho, const double z, const double z_prime, const layers_constants & LC);

Array<complex<double>, 1> G_EJ (const double rho, const double z, const double z_prime, const layers_constants & LC);
