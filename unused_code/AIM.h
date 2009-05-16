/***************************************************************************
 * AIM_S.cpp  computes the inverse of the Vandermonde matrix
 *            Based upon Bleszynski paper, Radio Science September-October 1996
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

void AIM_Q_numeric (Array<double, 3>& Phi, Array<double, 4>& Phi_r, Array<double, 4>& Phi_n_hat_X_r, const TinyVector<double, 3>& R, const Array<double, 2>&  vertexes_coord, const Array<int, 2>& triangles_vertexes, const int triangle_index, const int M, const int N_points);

void AIM_S (Array<double, 2>& S_x, const Array<double, 1>& x);
