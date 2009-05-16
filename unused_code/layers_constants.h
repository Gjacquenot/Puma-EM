/***************************************************************************
 * layers_constants.h  structure containing multilayer medium properties
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
using namespace blitz;

typedef struct {
  int N;
  double w;
  double k_0;
  double eps_0;
  double mu_0;
  Array<double, 1> z_i;
  Array<double, 1> d_i;
  Array<complex<double>, 1> mu_i;
  Array<complex<double>, 1> eps_i;
  Array<complex<double>, 1> k_i;
} layers_constants;
