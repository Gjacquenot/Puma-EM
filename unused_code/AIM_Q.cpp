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
#include <iostream>
#include <complex>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

using namespace blitz;

#include "triangle_struct.h"
#include "vector_functions.h"
#include "GK_triangle.h"
#include "binomial.h"

void AIM_Q_numeric (Array<double, 3>& Phi, Array<double, 4>& Phi_r, Array<double, 4>& Phi_n_hat_X_r, const TinyVector<double, 3>& R, const Array<double, 2>&  vertexes_coord, const Array<int, 2>& triangles_vertexes, const int triangle_index, const int M, const int N_points) {

  // R is the point around which the multipole expansion happens
  int q1, q2, q3, i, j;
  double term, phi, sum_weigths;
  TinyVector<double, 3> r, r_R, n_hat_X_r;
  Array<double, 1> phi_r(3), n_hat_X_phi_r(3);

  // choice of the weights and integration position points
  const double *xi, *eta, *weigths;
  IT_points (xi, eta, weigths, sum_weigths, N_points);

  // construction of the triangle
  triangle T;
  construct_triangle(T, vertexes_coord, triangles_vertexes, triangle_index);
  double norm_factor = T.A/sum_weigths; // 2.0 or not? that's the question...

  // we now compute Q(q1, q2, q3)
  for (q1=0 ; q1<M ; q1++) {
    for (q2=0 ; q2<M ; q2++) {
      for (q3=0 ; q3<M ; q3++) {

        phi = 0.0;
        phi_r = 0.0;
        for (j=0 ; j<N_points ; j++) {

          r = T.r_nodes(0) * xi[j] + T.r_nodes(1) * eta[j] + T.r_nodes(2) * (1.0-xi[j]-eta[j]);
          n_hat_X_r = cross(T.n_hat, r);
          r_R = r - R;
          term = pow(r_R(0), q1) * pow(r_R(1), q2) * pow(r_R(2), q3) * weigths[j];
          phi += term;

          // Phi_r: we must go from "TinyVector" to "Array" type
          for (i=0 ; i<3 ; i++) {
            phi_r(i) += r(i) * term;
            n_hat_X_phi_r(i) += n_hat_X_r(i) * term;
          } // i

        } // j
        Phi(q1, q2, q3) = phi * norm_factor;
        Phi_r(q1, q2, q3, Range::all()) = phi_r * norm_factor;
        Phi_n_hat_X_r(q1, q2, q3, Range::all()) = n_hat_X_phi_r * norm_factor;
      } // q3
    } // q2
  } //q1

}

// void AIM_Q_analytic (Array<double, 3>& Phi, Array<double, 4>& Phi_r, Array<double, 4>& Phi_n_hat_X_r, const TinyVector<double, 3>& R, const Array<double, 2>&  vertexes_coord, const Array<int, 2>& triangles_vertexes, const int triangle_index, const int M, const int N_points) {
// 
//   int q1, q2, q3, p1, p2, p3, k1, k2, k3, m;
//   double result_temp;
//   TinyVector<double, 3> r0_r2, r1_r2, r2_R;
// 
//   triangle T;
//   construct_triangle(T, vertexes_coord, triangles_vertexes, triangle_index);
//   r0_r2 = T.r_nodes(0) - T.r_nodes(2);
//   r1_r2 = T.r_nodes(1) - T.r_nodes(2);
//   r2_R = T.r_nodes(2) - R;
// 
//   for (q1=0 ; q1<M ; q1++) {
//     result_temp = 0;
//     for (q2=0 ; q2<M ; q2++) {
//       for (q3=0 ; q3<M ; q3++) {
// 
//         for (p1=0 ; p1<q1+1 ; p1++) {
//           for (p2=0 ; p2<q2+1 ; p2++) {
//             for (p3=0 ; p3<q3+1 ; p3++) {
// 
//               for (k1=0 ; k1<p1+1 ; k1++) {
//                 for (k2=0 ; k2<p2+1 ; k2++) {
//                   for (k3=0 ; k3<p3+1 ; k3++) {
// 
//                      for (m=0 ; m<k1+k2+k3+2 ; m++) result_temp += binomial (k1+k2+k3+1, m) * pow(-1, k1+k2+k3+1-m) * 1.0/(p1+p2+p3-m);
//                      result_temp /= (k1+k2+k3+1);
//                      result_temp *= binomial (p3, k3) * pow(r0_r2(0), k1) * pow(r1_r2(0), p1-k1)
// 
//   binomial (q1, p1) * binomial (q2, p2) * binomial (q3, p3);
// 
// 
// 
// }
