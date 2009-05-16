/***************************************************************************
 * generate_K_G_J_grid.h  function that generates 3D tables for the DGFs
 *                        those tables will be interpolated in the MoM
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
void generate_grid_points (Array<double, 1>& rho_values, Array<double, 1>& z_values, Array<double, 1>& z_prime_values, const double rho_min, const double rho_max, const double z_min, const double z_max, const double z_prime_min, const double z_prime_max, const double lambda_layer);

void generate_KA (const int USE_SYMMETRY, KA_grid & KA_tab, const Array<double, 1>& z_values, const Array<double, 1>& z_prime_values, const Array<double, 1>& rho_values, const int N_decim_z, const int N_decim_z_prime, const int N_decim_rho, const layers_constants & LC);

void generate_Kphi (const int USE_SYMMETRY, Kphi_grid & Kphi_tab, const Array<double, 1>& z_values, const Array<double, 1>& z_prime_values, const Array<double, 1>& rho_values, const int N_decim_z, const int N_decim_z_prime, const int N_decim_rho, const layers_constants & LC);

void generate_grad_Kphi (const int USE_SYMMETRY, Kphi_grid & Kphi_tab, const Array<double, 1>& z_values, const Array<double, 1>& z_prime_values, const Array<double, 1>& rho_values, const int N_decim_z, const int N_decim_z_prime, const int N_decim_rho, const layers_constants & LC);

void generate_grad_prime_Kphi (const int USE_SYMMETRY, Kphi_grid & Kphi_tab, const Array<double, 1>& z_values, const Array<double, 1>& z_prime_values, const Array<double, 1>& rho_values, const int N_decim_z, const int N_decim_z_prime, const int N_decim_rho, const layers_constants & LC);

void generate_G_HJ (const int USE_SYMMETRY, G_HJ_grid & G_HJ_tab, const Array<double, 1>& z_values, const Array<double, 1>& z_prime_values, const Array<double, 1>& rho_values, const int N_decim_z, const int N_decim_z_prime, const int N_decim_rho, const layers_constants & LC);

void generate_G_EJ (const int USE_SYMMETRY, G_EJ_grid & G_EJ_tab, const Array<double, 1>& z_values, const Array<double, 1>& z_prime_values, const Array<double, 1>& rho_values, const int N_decim_z, const int N_decim_z_prime, const int N_decim_rho, const layers_constants & LC);
