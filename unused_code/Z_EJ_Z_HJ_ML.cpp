/***************************************************************************
 * Z_EJ_Z_HJ_ML.cpp  function that computes the ML medium MoM matrix
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
#include <iomanip>
#include <complex>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

using namespace blitz;

const complex<double> I (0.0, 1.0);

#include "layers_constants.h"
#include "K_G_grid.h"
#include "mesh.h"
#include "triangle_struct.h"
#include "vector_functions.h"
#include "triangle_int.h"

void Z_EJ_Z_HJ_ML (Array<complex<double>,2>& Z_TE_J, Array<complex<double>,2>& Z_NE_J, Array<complex<double>,2>& Z_TH_J, Array<complex<double>,2>& Z_NH_J, const mesh & MESH, const KA_grid & KA_tab, const Kphi_grid & Kphi_tab, const G_HJ_grid & G_HJ_tab, const layers_constants & LC) {

  // def of k, mu_i, eps_i, obs layer
  int ind_layer = count(LC.z_i<KA_tab.z_values(0));
  int N_triangles = MESH.triangles.rows(), E = MESH.edges.rows()/2;
  int i, j, m, n, p, q, r, s, cpt, N_points_o, N_points_s, IS_TOUCH, IS_NEAR;
  int sign_edge_p, sign_edge_q;
  double l_p, l_q, C_pq, R_os;

  // declaration of all the geometrical vectors and triangles
  TinyVector<double, 3> r_oc, r_sc, r_q, r_p, n_hat_X_r_p;
  triangle T_src, T_obs;

  // declaration of the scalars, vectors and arrays needed in the integrations
  // Z_EJ
  complex<double> ITo_r_dot_ITs_KA_dot_rprime, ITo_n_hat_X_r_dot_ITs_KA_dot_rprime, p_dot_ITo_ITs_KA_dot_rprime, n_hat_X_r_p_dot_ITo_ITs_KA_dot_rprime, ITo_ITs_K_phi, ITo_n_hat_X_r_dot_ITs_grad_K_phi, n_hat_X_r_p_dot_ITo_ITs_grad_K_phi, a_pq, phi_pq;
  TinyVector<complex<double>, 3> ITo_r_dot_ITs_KA, ITo_ITs_KA_dot_rprime, ITo_n_hat_X_r_dot_ITs_KA, p_dot_ITo_ITs_KA, n_hat_X_r_p_dot_ITo_ITs_KA, ITo_ITs_grad_K_phi;
  Array<complex<double>, 2> ITo_ITs_KA (3, 3);

  // Z_HJ
  complex<double> ITo_r_dot_ITs_G_HJ_dot_rprime, ITo_n_hat_X_r_dot_ITs_G_HJ_dot_rprime, r_p_dot_ITo_ITs_G_HJ_dot_rprime, n_hat_X_r_p_dot_ITo_ITs_G_HJ_dot_rprime;
  TinyVector<complex<double>, 3> ITo_r_dot_ITs_G_HJ, ITo_ITs_G_HJ_dot_rprime, ITo_n_hat_X_r_dot_ITs_G_HJ, r_p_dot_ITo_ITs_G_HJ, n_hat_X_r_p_dot_ITo_ITs_G_HJ;
  Array<complex<double>, 2> ITo_ITs_G_HJ (3, 3);

  Z_TE_J = 0.0;
  Z_NE_J = 0.0;
  Z_TH_J = 0.0;
  Z_NH_J = 0.0;

  Array<triangle, 1> T_array (N_triangles);
  for (r=0 ; r<N_triangles ; r++) construct_triangle (T_array (r), MESH, r);

  cout << "computation of the multilayer Z matrix" << endl;
  for (r=0 ; r<N_triangles ; r++) { // loop on the observation triangles
    cout << "\r" << r*100/N_triangles << " \%";
    flush(cout);

    for (s=0 ; s<N_triangles ; s++) { // loop on the source triangles
      R_os = sqrt (dot (T_array (r).r_grav-T_array (s).r_grav, T_array (r).r_grav-T_array (s).r_grav));
      IS_TOUCH = (R_os - T_array (r).R_max - T_array (s).R_max <= 0.0);
      IS_NEAR = (R_os - 1.5*T_array (r).R_max - 1.5*T_array (s).R_max <= 0.0); 
      N_points_o = (N_points_s = 3);
      if (r==s || IS_TOUCH) {
        N_points_o = (N_points_s = 9);
      }
      else if (IS_NEAR) {
        N_points_o = (N_points_s = 6);
      }

      ITo_ITs_KA_ML (ITo_ITs_KA, ITo_r_dot_ITs_KA, ITo_ITs_KA_dot_rprime, ITo_r_dot_ITs_KA_dot_rprime, ITo_n_hat_X_r_dot_ITs_KA, ITo_n_hat_X_r_dot_ITs_KA_dot_rprime, KA_tab, T_array (r), T_array (s), LC.k_i(ind_layer), N_points_o, N_points_s);

      ITo_ITs_Kphi_ML (ITo_ITs_K_phi, ITo_ITs_grad_K_phi, ITo_n_hat_X_r_dot_ITs_grad_K_phi, Kphi_tab, T_array (r), T_array (s), LC.k_i(ind_layer), N_points_o, N_points_s);

      ITo_ITs_G_HJ_ML (ITo_ITs_G_HJ, ITo_r_dot_ITs_G_HJ, ITo_ITs_G_HJ_dot_rprime, ITo_r_dot_ITs_G_HJ_dot_rprime, ITo_n_hat_X_r_dot_ITs_G_HJ, ITo_n_hat_X_r_dot_ITs_G_HJ_dot_rprime, G_HJ_tab, T_array (r), T_array (s), LC.k_i(ind_layer), N_points_o, N_points_s);

      int E_and_J = (T_array (r).E_H_J_M (0) && T_array (s).E_H_J_M (2)), E_and_M = (T_array (r).E_H_J_M (0) && T_array (s).E_H_J_M (3));
      int H_and_J = (T_array (r).E_H_J_M (1) && T_array (s).E_H_J_M (2)), H_and_M = (T_array (r).E_H_J_M (1) && T_array (s).E_H_J_M (3));
      m = 0;
      while (MESH.ind_edges_of_tr_i(r, m) != -1) {
	j = MESH.ind_edges_of_tr_i(r, m);
        m++;
	p = MESH.edges(j, 1);
        r_node (r_p, MESH.nodes, MESH.edges(j, 3));
	l_p = MESH.edge_length(j);
        sign_edge_p = MESH.edges(j, 2);
        n_hat_X_r_p = cross (T_array (r).n_hat, r_p);

        // temporary elements for Z_TE_J
        p_dot_ITo_ITs_KA_dot_rprime = dot(r_p, ITo_ITs_KA_dot_rprime);
        for (n=0 ; n<3 ; n++) p_dot_ITo_ITs_KA (n) = r_p (0) * ITo_ITs_KA (0, n) + r_p (1) * ITo_ITs_KA (1, n) + r_p (2) * ITo_ITs_KA (2, n);

        // temporary elements for Z_NE_J
        n_hat_X_r_p_dot_ITo_ITs_KA_dot_rprime = dot (n_hat_X_r_p, ITo_ITs_KA_dot_rprime);
        for (n=0 ; n<3 ; n++) n_hat_X_r_p_dot_ITo_ITs_KA (n) = n_hat_X_r_p (0) * ITo_ITs_KA (0, n) + n_hat_X_r_p (1) * ITo_ITs_KA (1, n) + n_hat_X_r_p (2) * ITo_ITs_KA (2, n);
        n_hat_X_r_p_dot_ITo_ITs_grad_K_phi = dot (n_hat_X_r_p, ITo_ITs_grad_K_phi);

        // temporary elements for Z_TH_J
        r_p_dot_ITo_ITs_G_HJ_dot_rprime = dot (r_p, ITo_ITs_G_HJ_dot_rprime);
        for (n=0 ; n<3 ; n++) r_p_dot_ITo_ITs_G_HJ (n) = r_p (0) * ITo_ITs_G_HJ (0, n) + r_p (1) * ITo_ITs_G_HJ (1, n) + r_p (2) * ITo_ITs_G_HJ (2, n);

        // temporary elements for Z_NH_J
        n_hat_X_r_p_dot_ITo_ITs_G_HJ_dot_rprime = dot (n_hat_X_r_p, ITo_ITs_G_HJ_dot_rprime);
        for (n=0 ; n<3 ; n++) n_hat_X_r_p_dot_ITo_ITs_G_HJ (n) = n_hat_X_r_p (0) * ITo_ITs_G_HJ (0, n) + n_hat_X_r_p (1) * ITo_ITs_G_HJ (1, n) + n_hat_X_r_p (2) * ITo_ITs_G_HJ (2, n);

        n = 0;
        while (MESH.ind_edges_of_tr_i(s, n) != -1) {
          i = MESH.ind_edges_of_tr_i(s, n);
	  n++;
	  q = MESH.edges(i, 1);
          r_node (r_q, MESH.nodes, MESH.edges(i, 3));
          l_q = MESH.edge_length(i);
          sign_edge_q = MESH.edges(i, 2);
          C_pq = sign_edge_p*sign_edge_q*l_p*l_q/(4.0*T_array (r).A*T_array (s).A);

          // <f_m ; EFIE> : Z_TE_J computation
          a_pq = C_pq * ( ITo_r_dot_ITs_KA_dot_rprime - dot(ITo_r_dot_ITs_KA, r_q) - p_dot_ITo_ITs_KA_dot_rprime + dot(p_dot_ITo_ITs_KA, r_q) );
          phi_pq = 4.0 * C_pq * ITo_ITs_K_phi;
          if (E_and_J || H_and_M) Z_TE_J (p-1, q-1) += (I*LC.w * a_pq - I/LC.w * phi_pq);

          // <n x f_m ; EFIE> : Z_NE_J computation
          a_pq = C_pq * ( ITo_n_hat_X_r_dot_ITs_KA_dot_rprime - dot(ITo_n_hat_X_r_dot_ITs_KA, r_q) - n_hat_X_r_p_dot_ITo_ITs_KA_dot_rprime + dot (n_hat_X_r_p_dot_ITo_ITs_KA, r_q) );
          phi_pq = 2.0 * C_pq * (ITo_n_hat_X_r_dot_ITs_grad_K_phi - n_hat_X_r_p_dot_ITo_ITs_grad_K_phi);
          if (E_and_J || H_and_M) Z_NE_J (p-1, q-1) += (I*LC.w * a_pq + I/LC.w * phi_pq);

          // <f_m ; MFIE> : Z_TH_J computation
          if (H_and_J || E_and_M) Z_TH_J (p-1, q-1) -= C_pq * ( ITo_r_dot_ITs_G_HJ_dot_rprime - dot (r_q, ITo_r_dot_ITs_G_HJ) - r_p_dot_ITo_ITs_G_HJ_dot_rprime + dot (r_p_dot_ITo_ITs_G_HJ, r_q) );

          // <n x f_m ; MFIE> : Z_NH_J computation
	  if (H_and_J || E_and_M) Z_NH_J (p-1, q-1) -= C_pq * ( ITo_n_hat_X_r_dot_ITs_G_HJ_dot_rprime - dot (ITo_n_hat_X_r_dot_ITs_G_HJ, r_q) - n_hat_X_r_p_dot_ITo_ITs_G_HJ_dot_rprime + dot (n_hat_X_r_p_dot_ITo_ITs_G_HJ, r_q) );

	} // for i

      } // for j

    } // for s

  } // for r

  printf("\n");
}
