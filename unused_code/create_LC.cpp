/***************************************************************************
 * create_LC.cpp  function that initializes object "layers_constants.h"
 *                needs a file "layers_param.ini"
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
#include <fstream>
#include <complex>
#include <blitz/array.h>

using namespace blitz;

const complex<double> I (0.0, 1.0);

#include "layers_constants.h"

void create_LC (layers_constants & LC, const double eps_0, const double mu_0) {

  int j;
  Array<double, 2> layers_param;

  // reading layers_param.ini
  {
    char * filename = "layers_param.ini";
    ifstream ifs (filename);
    if (! ifs.is_open()) { 
      cout << "create_LC : Error opening file : " << filename << endl; 
      exit (1); 
    }
    else ifs >> layers_param;
    ifs.close();
  }

  LC.eps_0 = eps_0;
  LC.mu_0 = mu_0;
  LC.N = (int) layers_param(1,0);
  LC.w = 2*M_PI*layers_param(0,0);
  LC.k_0 = LC.w*sqrt(eps_0*mu_0);
  if (LC.N>1) {
    LC.z_i.resize(LC.N-1);
    LC.d_i.resize(LC.N);
    LC.mu_i.resize(LC.N);
    LC.eps_i.resize(LC.N);
    LC.k_i.resize(LC.N);
    for (j=0 ; j<LC.N ; j++) {
      if (j<LC.N-1) LC.z_i(j) = layers_param(j+2,4);
      LC.eps_i(j) = layers_param(j+2,0)-I*layers_param(j+2,1);
      LC.mu_i(j) = layers_param(j+2,2)-I*layers_param(j+2,3);
      LC.k_i(j) = LC.k_0*sqrt(LC.eps_i(j)*LC.mu_i(j));
    }
    LC.d_i (Range (1, LC.N-2))= LC.z_i (Range(1, LC.N-2)) - LC.z_i (Range(0, LC.N-3));
    LC.d_i (0) = 0.0;
    LC.d_i (LC.N-1) = 0.0;
  }
  else {
    LC.z_i.resize(1);
    LC.d_i.resize(LC.N);
    LC.mu_i.resize(LC.N);
    LC.eps_i.resize(LC.N);
    LC.k_i.resize(LC.N);

    j = 0;
    LC.z_i(j) = 0.0;
    LC.eps_i(j) = layers_param(j+2,0)-I*layers_param(j+2,1);
    LC.mu_i(j) = layers_param(j+2,2)-I*layers_param(j+2,3);
    LC.k_i(j) = LC.k_0*sqrt(LC.eps_i(j)*LC.mu_i(j));
  }
}
