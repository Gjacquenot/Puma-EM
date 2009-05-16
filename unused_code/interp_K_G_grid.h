/***************************************************************************
 * interp_K_G_grid.h  routine for DGF tables interpolation
 *                    the grids are regularly spaced
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
void interp_KA_mm (Array<complex<double>, 2>& KA, const KA_grid & KA_tab, const TinyVector<double,3>& r, const TinyVector<double,3>& r_prime);

void interp_Kphi_mm (complex<double>& Kphi, const Kphi_grid & Kphi_tab, const TinyVector<double,3>& r, const TinyVector<double,3>& r_prime);

void interp_grad_K_phi (TinyVector<complex<double>, 3>& grad_K_phi, const Kphi_grid & Kphi_tab, const TinyVector<double,3>& r, const TinyVector<double,3>& r_prime);

void interp_grad_prime_K_phi (TinyVector<complex<double>, 3>& grad_prime_K_phi, const Kphi_grid & Kphi_tab, const TinyVector<double,3>& r, const TinyVector<double,3>& r_prime);

void interp_G_HJ_mm (Array<complex<double>, 2>& G_HJ, const G_HJ_grid & G_HJ_tab, const TinyVector<double,3>& r, const TinyVector<double,3>& r_prime);

void interp_G_EJ_mm (Array<complex<double>, 2>& G_EJ, const G_EJ_grid & G_EJ_tab, const TinyVector<double,3>& r, const TinyVector<double,3>& r_prime);
