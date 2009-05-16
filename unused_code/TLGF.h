/***************************************************************************
 * TLGF.h  Transmission Line Green's Functions
 *         Based upon Michalski paper, AP-IEEE 1997 
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

complex<double> V_i (const Array<complex<double>, 1>& Z_i, const Array<complex<double>, 1>& k_z_i, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> V_v (const Array<complex<double>, 1>& Z_i, const Array<complex<double>, 1>& k_z_i, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> I_v (const Array<complex<double>, 1>& Z_i, const Array<complex<double>, 1>& k_z_i, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

complex<double> I_i (const Array<complex<double>, 1>& Z_i, const Array<complex<double>, 1>& k_z_i, const double z, const double z_prime, const int m, const int n, const layers_constants & LC);

