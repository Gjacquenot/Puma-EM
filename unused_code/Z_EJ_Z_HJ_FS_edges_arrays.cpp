/***************************************************************************
 * Z_EJ_Z_HJ_FS.cpp  function that computes the free space MoM matrix
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
const double mu_0 = M_PI*4.0e-7;
const double c = 2.99792458e8;
const double eps_0 = 1.0/(c*c*mu_0);

#include "K_G_grid.h"
#include "triangle_int.h"

void Z_EJ_Z_HJ_FS_edges_arrays (Array<complex<double>, 2>& Z_TE_J, Array<complex<double>, 2>& Z_NE_J, Array<complex<double>, 2>& Z_TH_J, Array<complex<double>, 2>& Z_NH_J, const Array<double, 2>& vertexes_coord, const Array<int, 2>& triangles_vertexes, const Array<int, 2>& triangles_edges_indexes, const Array<int, 2>& edges_vertexes, const Array<int, 1>& edges_kinds, const Array<int, 1>& edges_numbers, const Array<int, 1>& edges_signs, const Array<double, 1>& edges_lengths, const double w, const complex<double>& eps_r, const complex<double>& mu_r) {

  // def of k, mu_i, eps_i, obs layer
  int N_triangles = triangles_vertexes.rows();
  int p, q, r, s, N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R, IS_TOUCH, IS_NEAR;
  int sign_edge_p, sign_edge_q, index_edge_p, index_edge_q, number_edge_p, number_edge_q;
  double l_p, l_q, C_pq, R_os;
  complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu);

  // declaration of all the geometrical vectors
  TinyVector<double, 3> r_oc, r_sc, r_q, r_p, n_hat_X_r_p;

  // declaration of the scalars and vectors needed in the integrations
  double IT_r_square, n_hat_X_r_p_dot_IT_r;
  TinyVector<double, 3> IT_r, IT_n_hat_X_r;

  complex<double> ITo_ITs_G, ITo_r_dot_ITs_G_rprime, p_dot_ITo_ITs_G_rprime, ITo_n_hat_X_r_dot_ITs_G_rprime, n_hat_X_r_p_dot_ITo_ITs_G_rprime, IDTo_m_hat_dot_n_hat_X_r_ITs_G, n_hat_X_r_p_dot_IDTo_m_hat_ITs_G, IDTo_l_hat_dot_r_ITs_G;
  complex<double> a_pq, phi_pq, z_pq;
  TinyVector<complex<double>, 3> ITo_r_ITs_G, ITo_ITs_G_rprime, ITo_n_hat_X_r_ITs_G, IDTo_m_hat_ITs_G, IDTo_l_hat_ITs_G;

  complex<double> ITo_n_hat_X_r_dot_r_X_ITs_grad_G, n_hat_X_r_p_dot_ITo_r_X_ITs_grad_G, n_hat_dot_ITo_r_X_ITs_grad_G;
  TinyVector<complex<double>, 3> ITo_ITs_grad_G, ITo_r_X_ITs_grad_G, ITo_n_hat_X_r_X_ITs_grad_G, r_p_X_ITo_ITs_grad_G, n_hat_X_r_p_X_ITo_ITs_grad_G;

  Array<triangle, 1> T_array (N_triangles);
  for (r=0 ; r<N_triangles ; r++) construct_triangle (T_array (r), vertexes_coord, triangles_vertexes, r);
  
  cout << "computation of the free-space Z matrix: " << N_triangles << " triangles" << endl;
  for (r=0 ; r<N_triangles ; r++) { // loop on the observation triangles
    cout << "\r" << r*100/N_triangles << " \%";
    flush(cout);

    IT_fm_fn (IT_r_square, IT_r, T_array (r)); // serves for <f_m ; f_n> and <f_m ; n x f_n>
    IT_n_hat_X_r = cross (T_array (r).n_hat, IT_r); // serves for <f_m ; n x f_n>

    for (s=0 ; s<N_triangles ; s++) { // loop on the source triangles
      R_os = sqrt (dot (T_array(r).r_grav-T_array(s).r_grav, T_array(r).r_grav-T_array(s).r_grav));
      IS_TOUCH = (R_os - T_array(r).R_max - T_array(s).R_max <= 0.0);
      IS_NEAR = (R_os - 1.5*T_array(r).R_max - 1.5*T_array(s).R_max <= 0.0); 

      EXTRACT_1_R = (EXTRACT_R = 0); N_points_o = (N_points_s = 3);
      if (r==s || IS_TOUCH) {
        EXTRACT_1_R = (EXTRACT_R = 1); N_points_o = (N_points_s = 13);
      }
      else if (IS_NEAR) {
        EXTRACT_1_R = (EXTRACT_R = 1); N_points_o = (N_points_s = 9);
      }

      ITo_ITs_free (ITo_ITs_G, ITo_r_ITs_G, ITo_ITs_G_rprime, ITo_r_dot_ITs_G_rprime, ITo_n_hat_X_r_ITs_G, ITo_n_hat_X_r_dot_ITs_G_rprime, ITo_ITs_grad_G, ITo_r_X_ITs_grad_G, ITo_n_hat_X_r_dot_r_X_ITs_grad_G, ITo_n_hat_X_r_X_ITs_grad_G, T_array (r), T_array (s), k, N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R);

      n_hat_dot_ITo_r_X_ITs_grad_G = dot (T_array (r).n_hat, ITo_r_X_ITs_grad_G);

      if (IS_NEAR) IDTo_ITs_free (IDTo_m_hat_dot_n_hat_X_r_ITs_G, IDTo_m_hat_ITs_G, IDTo_l_hat_dot_r_ITs_G, IDTo_l_hat_ITs_G, T_array (r), T_array (s), k, 15, N_points_s, EXTRACT_1_R, EXTRACT_R);

      for (p=0 ; p<3 ; p++) { // loop on the test basis
	index_edge_p = triangles_edges_indexes(r, p);
	if (edges_kinds(index_edge_p) > 1) { // if the edge is not border
	  number_edge_p = edges_numbers(index_edge_p); // this has to be modified!!
	  r_node(r_p, vertexes_coord, edges_vertexes(index_edge_p, 2));
	  l_p = edges_lengths(index_edge_p);
	  sign_edge_p = edges_signs(index_edge_p);
	  n_hat_X_r_p = cross(T_array(r).n_hat, r_p);

	  // temporary elements for Z_TE_J
	  p_dot_ITo_ITs_G_rprime = dot (r_p, ITo_ITs_G_rprime);

	  // temporary elements for Z_NE_J
	  n_hat_X_r_p_dot_ITo_ITs_G_rprime = dot(n_hat_X_r_p, ITo_ITs_G_rprime);
	  n_hat_X_r_p_dot_IDTo_m_hat_ITs_G = dot(n_hat_X_r_p, IDTo_m_hat_ITs_G);

	  // temporary elements for Z_TH_J
	  n_hat_X_r_p_dot_IT_r = dot(n_hat_X_r_p, IT_r);
	  r_p_X_ITo_ITs_grad_G = cross_real_complex(r_p, ITo_ITs_grad_G);

	  // temporary elements for Z_NH_J
	  n_hat_X_r_p_X_ITo_ITs_grad_G = cross_real_complex (n_hat_X_r_p, ITo_ITs_grad_G);
	  n_hat_X_r_p_dot_ITo_r_X_ITs_grad_G = dot (n_hat_X_r_p, ITo_r_X_ITs_grad_G);

	  for (q=0 ; q<3 ; q++) { // loop on the source basis
	    index_edge_q = triangles_edges_indexes(s, q);
	    if (edges_kinds(index_edge_q) > 1) { // if the edge is not border
	      number_edge_q = edges_numbers(index_edge_q); // this has to be modified!!
	      r_node(r_q, vertexes_coord, edges_vertexes(index_edge_q, 2));
	      l_q = edges_lengths(index_edge_q);
	      sign_edge_q = edges_signs(index_edge_q);
	      C_pq = sign_edge_p * sign_edge_q * l_p * l_q/(4.0 * T_array(r).A * T_array(s).A);

	      // <f_p ; EFIE> : Z_TE_J computation
	      a_pq = C_pq * mu * (ITo_r_dot_ITs_G_rprime - dot(ITo_r_ITs_G, r_q) - p_dot_ITo_ITs_G_rprime + dot(r_p, r_q) * ITo_ITs_G);
	      phi_pq = 4.0 * C_pq/eps * ITo_ITs_G;
	      z_pq = (I* w * a_pq - I/w * phi_pq);
	      Z_TE_J (number_edge_p, number_edge_q) += z_pq;

	      // <n x f_p ; EFIE> : Z_NE_J computation
	      a_pq = C_pq * mu * (ITo_n_hat_X_r_dot_ITs_G_rprime - dot(ITo_n_hat_X_r_ITs_G, r_q) - n_hat_X_r_p_dot_ITo_ITs_G_rprime + dot (n_hat_X_r_p, r_q) * ITo_ITs_G);
	      if (IS_NEAR) phi_pq = 2.0 * C_pq/eps * (-IDTo_l_hat_dot_r_ITs_G + dot(r_p, IDTo_l_hat_ITs_G));
	      else phi_pq = 2.0 * C_pq/eps * ( n_hat_dot_ITo_r_X_ITs_grad_G - dot(n_hat_X_r_p, ITo_ITs_grad_G));
	      z_pq = (I * w * a_pq + I/w * phi_pq);
	      Z_NE_J (number_edge_p, number_edge_q) += z_pq;

	      // <f_p ; MFIE> : Z_TH_J computation
	      if (r!=s) z_pq = C_pq * dot (r_p-r_q, ITo_r_X_ITs_grad_G - r_p_X_ITo_ITs_grad_G);
	      else z_pq = 0.0;
	      Z_TH_J (number_edge_p, number_edge_q) -= z_pq;

	      // <n x f_p ; MFIE> : Z_NH_J computation
	      if (r!=s) z_pq = C_pq * ( ITo_n_hat_X_r_dot_r_X_ITs_grad_G + dot (r_q, ITo_n_hat_X_r_X_ITs_grad_G - n_hat_X_r_p_X_ITo_ITs_grad_G) - n_hat_X_r_p_dot_ITo_r_X_ITs_grad_G  );
	      else z_pq = 0.0;
	      Z_NH_J (number_edge_p, number_edge_q) += z_pq;
	    } // end if
	  } // end for(q...)

	} // end if
      } // end for(p...)

    } // for s

  } // for r
  printf("\n");
}
