/***************************************************************************
 * S_v_K_spectral.h    Sommerfeld integration of spectral domain DGFs
 *                     Based upon Michalski paper, AP-IEEE 1997-1998 
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
complex<double> S_v_K_spectral (complex<double> (*f)(const complex<double>, const double, const double, const double, const double, const int, const int, const layers_constants &), const double v, const double rho, const double z, const double z_prime, const layers_constants & LC);
