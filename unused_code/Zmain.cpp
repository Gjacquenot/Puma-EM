/***************************************************************************
 * Zmain.cpp  main function that creates the DGF grids and the MoM matrices
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
#include <iostream>
#include <fstream>
#include <complex>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

using namespace blitz;

const complex<double> I (0.0, 1.0);
const double mu_0 = M_PI*4.0e-7;
const double c = 2.99792458e8;
const double eps_0 = 1.0/(c*c*mu_0);

#include "layers_constants.h"
#include "create_LC.h"
#include "K_G_grid.h"
#include "mesh.h"
#include "create_mesh.h"
#include "generate_K_G_J_grid.h"
#include "interp_K_G_grid.h"
#include "Z_EJ_Z_HJ.h"

int is_multilayer (layers_constants & LC) {
  int i, M = 0;
  for (i=0 ; i<LC.N ; i++) M += ( (real(LC.eps_i(i))==real(LC.eps_i(0))) && (imag(LC.eps_i(i))==imag(LC.eps_i(0))) );
  return (M!=LC.N);
}

void Zmain (Array<complex<double>, 2>& Z_TE_J, Array<complex<double>, 2>& Z_NE_J, Array<complex<double>, 2>& Z_TH_J, Array<complex<double>, 2>& Z_NH_J, Array<complex<double>, 2>& Z_TE_M, Array<complex<double>, 2>& Z_NE_M, Array<complex<double>, 2>& Z_TH_M, Array<complex<double>, 2>& Z_NH_M) {

  int i, j, m, n, p;

  layers_constants LC;
  create_LC (LC, eps_0, mu_0);
  cout << "z_i = " << LC.z_i << endl;
  cout << "eps_i = " << LC.eps_i << endl;
  cout << "mu_i = " << LC.mu_i << endl;
  cout << "k_i = " << LC.k_i << endl;

  mesh MESH_TARGET;
  char * nodes_target = "nodes_target.txt", * triangles_target = "triangles_target.txt", * edges_target = "edges_target.txt";
  create_mesh (MESH_TARGET, nodes_target, triangles_target, edges_target);

  int E = MESH_TARGET.edges.rows()/2, MULTILAYER_MEDIUM = (LC.N>1 && is_multilayer(LC)), DIELECTRIC_TARGET = 1;
  cout << "MULTILAYER_MEDIUM = " << MULTILAYER_MEDIUM << endl;
  Z_TE_J.resize (E, E);
  Z_NE_J.resize (E, E);
  Z_TH_J.resize (E, E);
  Z_NH_J.resize (E, E);
  Z_TE_M.resize (E, E);
  Z_NE_M.resize (E, E);
  Z_TH_M.resize (E, E);
  Z_NH_M.resize (E, E);
  //Z_EJ_Z_HJ_Z_EM_Z_HM_FS (Z_TE_J, Z_NE_J, Z_TH_J, Z_NH_J, Z_TE_M, Z_NE_M, Z_TH_M, Z_NH_M, MESH_TARGET, LC);
  
  if (MULTILAYER_MEDIUM) {
    cout << "Multilayer part" << endl;
    cout << "Generation of multilayered media green functions for J currents" << endl;

    // Green function domain definition
    double eps_limit = 1.0e-8;
    double z_min = MESH_TARGET.z_min-eps_limit, z_max = MESH_TARGET.z_max+eps_limit, z_prime_min = MESH_TARGET.z_min-eps_limit, z_prime_max = MESH_TARGET.z_max+eps_limit;
    double rho_max = MESH_TARGET.rho_max+eps_limit, rho_min = 0.0;
    cout << "z_min = " << z_min << ", z_max = " << z_max << ", rho_max = " << rho_max << endl;

    // wavelength in layer containing the target
    int ind_layer = count(LC.z_i<z_min);
    if (ind_layer!=count(LC.z_i<z_max)) cout << "BIG PROBLEM : OBJECT CROSSES INTERFACES!!" << endl;
    double lambda_layer = 2 * M_PI/real(LC.k_i (ind_layer));

    // construction of the points on which the DGF will be evaluated
    Array<double, 1> rho_values, z_values, z_prime_values;
    generate_grid_points (rho_values, z_values, z_prime_values, rho_min, rho_max, z_min, z_max, z_prime_min, z_prime_max, lambda_layer);

    // grid construction
    int USE_SYMMETRY = 1;
    int N_decim_z = (z_values.size()>3) ? 1 : 0, N_decim_z_prime = (z_prime_values.size()>3) ? 1 : 0, N_decim_rho = (rho_values.size()>3) ? 1 : 0;
    KA_grid KA_tab_object_object;
    generate_KA (USE_SYMMETRY, KA_tab_object_object, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC);

    Kphi_grid Kphi_tab_object_object;
    generate_Kphi (USE_SYMMETRY, Kphi_tab_object_object, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC);
    generate_grad_Kphi (USE_SYMMETRY, Kphi_tab_object_object, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC);

    G_HJ_grid G_HJ_tab_object_object;
    generate_G_HJ (USE_SYMMETRY, G_HJ_tab_object_object, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC);

    {
      Array<complex<double>, 2> Z_TE_J_ML (E, E), Z_NE_J_ML (E, E), Z_TH_J_ML (E, E), Z_NH_J_ML (E, E);
      Z_EJ_Z_HJ_ML (Z_TE_J_ML, Z_NE_J_ML, Z_TH_J_ML, Z_NH_J_ML, MESH_TARGET, KA_tab_object_object, Kphi_tab_object_object, G_HJ_tab_object_object, LC);
      Z_TE_J += Z_TE_J_ML;
      Z_NE_J += Z_NE_J_ML;
      Z_TH_J += Z_TH_J_ML;
      Z_NH_J += Z_NH_J_ML;
    }

    // dual green functions for M currents
    if (DIELECTRIC_TARGET == 1) {
      cout << "Generation of multilayered media dual green functions for M currents" << endl;
      layers_constants LC_dual;
      create_LC (LC_dual, mu_0, eps_0);
      LC_dual.eps_i = LC.mu_i;
      LC_dual.mu_i = LC.eps_i;

      KA_grid KA_tab_object_object_dual;
      generate_KA (USE_SYMMETRY, KA_tab_object_object_dual, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC_dual);

      Kphi_grid Kphi_tab_object_object_dual;
      generate_Kphi (USE_SYMMETRY, Kphi_tab_object_object_dual, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC_dual);
      generate_grad_Kphi (USE_SYMMETRY, Kphi_tab_object_object_dual, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC_dual);

      G_HJ_grid G_HJ_tab_object_object_dual;
      generate_G_HJ (USE_SYMMETRY, G_HJ_tab_object_object_dual, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC_dual);

      {
        Array<complex<double>, 2> Z_TE_M_ML (E, E), Z_NE_M_ML (E, E), Z_TH_M_ML (E, E), Z_NH_M_ML (E, E);
        Z_EJ_Z_HJ_ML (Z_TH_M_ML, Z_NH_M_ML, Z_TE_M_ML, Z_NE_M_ML, MESH_TARGET, KA_tab_object_object_dual, Kphi_tab_object_object_dual, G_HJ_tab_object_object_dual, LC_dual);
        Z_TE_M -= Z_TE_M_ML;
        Z_NE_M -= Z_NE_M_ML;
        Z_TH_M += Z_TH_M_ML;
        Z_NH_M += Z_NH_M_ML;
      }

    }

  } // end (MULTILAYER_MEDIUM == 1)

}
