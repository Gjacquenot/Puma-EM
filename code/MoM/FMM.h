#ifndef FMM_H
#define FMM_H

#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

using namespace blitz;

#include "EMConstants.h"
#include "integr_1D_X_W.h"
#include "GK_triangle.h"
#include "./amos/zbesh/zbesh_interface.h"

inline void rNode(TinyVector<double,3>& r, const Array<double,2>& vertexes_coord, const int node_index) {
  for (int i=0 ; i<3 ; i++) r(i) = vertexes_coord(node_index, i);
}

inline void rNode(TinyVector<float,3>& r, const Array<double,2>& vertexes_coord, const int node_index) {
  for (int i=0 ; i<3 ; i++) r(i) = static_cast<float>(vertexes_coord(node_index, i));
}

template <typename T>
void P_Legendre (Array<T, 1>& P,
                 const T z)
/**
 * This function computes the Legendre polynomials, given by the following recursive formula:
 *
 *  \f[ P \left( i+1 \right) = \big( \left(2 i + 1\right) z P\left(i\right) - i P\left(i-1\right) \big)/\left(i+1\right) \f]
 *
 * with
 *
 * \f[ P(0) = 1, P(1) = z \f]
 */
{
  int L = P.size() - 1;
  P(0) = 1;
  if (L >= 1) {
    P(1) = z;
    for (int i=1 ; i<L ; i++) P(i+1) = ( (2.0*i+1.0) * z * P(i) - i * P(i-1) )/(i+1.0);
  }
}


/***********************************************************/
/******************** alpha computation ********************/
/***********************************************************/


/***********************************************************/
/********************** B computation **********************/
/***********************************************************/

// void IT_g_IT_g_r (TinyVector<complex<double>, 3>& IT_g_r, /**< OUTPUT: \f$ \int_T \mathbf{r} \; e^{-jk \left(\widehat{\mathbf{k}} \cdot \mathbf{r} \right)} \f$ */
//                   complex<double>& IT_g, /**< OUTPUT: \f$ \int_T e^{-jk \left(\widehat{\mathbf{k}} \cdot \mathbf{r} \right)} \f$ */
//                   const double exp_arg_sign, /**< INPUT: sign to be taken in argument of exponential */
//                   const TinyVector<double, 3>& r0, /**< INPUT: position vector of node 0 */
//                   const TinyVector<double, 3>& r1, /**< INPUT: position vector of node 1 */
//                   const TinyVector<double, 3>& r2, /**< INPUT: position vector of node 2 */
//                   const double triangle_area, /**< INPUT: area of triangle */
//                   const TinyVector<double, 3>& k_hat, /**< INPUT: unit vector \f$ \widehat{\mathbf{k}} \f$ of plane wave */
//                   const complex<double>& k, /**< INPUT: the wavenumber */
//                   const int N_points) /**< INPUT: number of Gaussian points for triangle integration */
// {
//   double sum_weigths, norm_factor;
//   complex<double> I_k_r, IT_g_j;
//   TinyVector<double, 3> r;
//   const double *xi, *eta, *weigths;
//   IT_points (xi, eta, weigths, sum_weigths, N_points);
// 
//   IT_g = 0.0;
//   IT_g_r = 0.0;
// 
//   for (int j=0 ; j<N_points ; j++) {
//     r = r0 * xi[j] + r1 * eta[j] + r2 * (1-xi[j]-eta[j]);
//     I_k_r = I * k * dot(k_hat, r);
//     IT_g_j = exp(exp_arg_sign * I_k_r) * weigths[j];
// 
//     IT_g += IT_g_j;
//     IT_g_r += IT_g_j * r;
//   }
//   norm_factor = triangle_area/sum_weigths;
//   IT_g *= norm_factor;
//   IT_g_r *= norm_factor;
// }
// 
// 
// inline void cart2sph (TinyVector<complex<double>, 2>& A_sph, 
//                       const TinyVector<complex<double>, 3>& A_cart, 
//                       const double cos_theta, 
//                       const double phi) 
// {
//   double sin_theta = sqrt(1.0-pow2(cos_theta));
//   A_sph = A_cart(0) * cos_theta*cos(phi) + A_cart(1) * cos_theta*sin(phi) - A_cart(2) * sin_theta,
//           -A_cart(0) * sin(phi) + A_cart(1) * cos(phi);
// }
// 
// template <typename T>
// void B_EJ_B_HJ (Array<complex<T>, 4>& B, /**< OUTPUT: 4D Array \f$ B( RWG, coord, \theta,\phi) \f$ */
//                 const Array<complex<double>, 1>& CFIE, /**< combination factors for integral equation */
//                 const Array<int, 1>& indexes_triangles, 
//                 const Array<int, 1>& edges_numbers_local_edges_numbers, 
//                 const Array<double, 2>& vertexes_coord, 
//                 const Array<int, 2>& triangles_vertexes, 
//                 const Array<int, 2>& triangles_edges_kinds, 
//                 const Array<int, 2>& triangles_edges_numbers, 
//                 const Array<int, 2>& triangles_edges_signs, 
//                 const Array<double, 2>& triangles_edges_lengths, 
//                 const Array<int, 2>& triangles_edges_opp_vertexes, 
//                 const Array<double, 2>& triangles_normals, 
//                 const Array<double, 2>& triangles_areas, 
//                 const Array<double, 2>& edges_numbers_cubes_centroids, 
//                 const int IS_SRC, /**< 1 if we compute the radiation function, 0 if receiving function */
//                 const complex<double>& k, /**< INPUT: the wavenumber */
//                 const double w, /**< INPUT: pulsation */
//                 const complex<double>& mu_r, /**< INPUT: relative magnetic permeability */
//                 const double theta, /**< INPUT: \f$ \theta \f$ */
//                 const double phi, /**< INPUT: \f$ \phi \f$ */
//                 const int index_theta, /**< INPUT: index of \f$ \theta \f$ */
//                 const int index_phi, /**< INPUT: index of \f$ \phi \f$ */
//                 const int N_points) /**< INPUT: the number of points used in Gaussian integration */
// {
//   int N_triangles = indexes_triangles.size(), index_T, number_edge_p, local_number_edge_p;
//   double exp_arg_sign = (IS_SRC==1) ? 1.0 : -1.0, sin_theta = sin(theta), cos_theta = cos(theta);
//   Range all = Range::all();
// 
//   // geometrical entities
//   TinyVector<double, 3> r0, r1, r2, r_p, r_cube_centroid, n_hat, n_hat_X_r_p;
//   TinyVector<double, 3> k_hat, theta_hat, phi_hat;
//   k_hat = sin_theta*cos(phi), sin_theta*sin(phi), cos_theta;
//   theta_hat = cos_theta*cos(phi), cos_theta*sin(phi), -sin_theta;
//   phi_hat = -sin(phi), cos(phi), 0.0;
// 
//   int sign_edge_p;
//   double l_p, C_p;
//   complex<double> IT_g;
//   TinyVector<complex<double>, 3> IT_g_r, b_p, n_hat_X_IT_g_r, n_hat_X_b_p;
//   TinyVector<complex<double>, 3> I_kk_dot_b_p_cartesian, I_kk_dot_n_hat_X_b_p_cartesian, k_hat_X_b_p_cartesian, k_hat_X_n_hat_X_b_p_cartesian;
//   TinyVector<complex<double>, 2> I_kk_dot_b_p_spherical, I_kk_dot_n_hat_X_b_p_spherical, k_hat_X_b_p_spherical, k_hat_X_n_hat_X_b_p_spherical, resultTmp;
// 
//   for (int r=0 ; r<N_triangles ; r++) { // loop on the triangles
//     index_T = indexes_triangles(r);
//     rNode(r0, vertexes_coord, triangles_vertexes(index_T, 0));
//     rNode(r1, vertexes_coord, triangles_vertexes(index_T, 1));
//     rNode(r2, vertexes_coord, triangles_vertexes(index_T, 2));
//     rNode(n_hat, triangles_normals, index_T);
//     IT_g_IT_g_r (IT_g_r, IT_g, exp_arg_sign, r0, r1, r2, triangles_areas(index_T, 0), k_hat, k, N_points);
//     n_hat_X_IT_g_r = n_hat(1)*IT_g_r(2)-n_hat(2)*IT_g_r(1),
//                      n_hat(2)*IT_g_r(0)-n_hat(0)*IT_g_r(2),
//                      n_hat(0)*IT_g_r(1)-n_hat(1)*IT_g_r(0);
// 
//       for (int p=0 ; p<3 ; p++) { // loop on the basis
// 
// 	if (triangles_edges_kinds(index_T, p) > 1) { // if the edge is not border
// 	  number_edge_p = triangles_edges_numbers(index_T, p); 
// 	  local_number_edge_p = edges_numbers_local_edges_numbers(number_edge_p);
// 	  if (local_number_edge_p > -1) { // if the edge has received a local number (which enables its participation)
// 
//             rNode(r_cube_centroid, edges_numbers_cubes_centroids, number_edge_p);
//             rNode(r_p, vertexes_coord, triangles_edges_opp_vertexes(index_T, p));
//             n_hat_X_r_p = cross(n_hat, r_p);
// 
// 	    l_p = triangles_edges_lengths(index_T, p);
// 	    sign_edge_p = triangles_edges_signs(index_T, p);
// 	    C_p = sign_edge_p * l_p / (2.0 * triangles_areas(index_T, 0));
//             complex<double> b_p_tmp(exp(-exp_arg_sign * I * k * dot(k_hat, r_cube_centroid)) * C_p);
//             b_p = b_p_tmp * (IT_g_r - IT_g * r_p);
//             I_kk_dot_b_p_spherical = dot(theta_hat, b_p), dot(phi_hat, b_p);
// 
//             if (IS_SRC==1) {
//               B(local_number_edge_p, 0, index_theta, index_phi) += static_cast< complex<T> > (I_kk_dot_b_p_spherical(0));
//               B(local_number_edge_p, 1, index_theta, index_phi) += static_cast< complex<T> > (I_kk_dot_b_p_spherical(1));
//             }
//             else if (IS_SRC==0) { // we perform the hereafter computations
// 			          // only if we are computing B_test arrays
//               n_hat_X_b_p = b_p_tmp * (n_hat_X_IT_g_r - IT_g * n_hat_X_r_p);
//               I_kk_dot_n_hat_X_b_p_spherical = dot(theta_hat, n_hat_X_b_p), dot(phi_hat, n_hat_X_b_p);
//               k_hat_X_b_p_cartesian =  k_hat(1)*b_p(2)-k_hat(2)*b_p(1),
//                                        k_hat(2)*b_p(0)-k_hat(0)*b_p(2),
//                                        k_hat(0)*b_p(1)-k_hat(1)*b_p(0);
//               cart2sph (k_hat_X_b_p_spherical, k_hat_X_b_p_cartesian, cos_theta, phi);
// 
//               k_hat_X_n_hat_X_b_p_cartesian = dot(k_hat, b_p) * n_hat - dot(k_hat, n_hat) * b_p;
//               cart2sph (k_hat_X_n_hat_X_b_p_spherical, k_hat_X_n_hat_X_b_p_cartesian, cos_theta, phi);
//               /// CFIE : [tEJ, nEJ, tHJ, nHJ]
//               resultTmp = (CFIE(0) * -I * w * mu_0 * mu_r) * I_kk_dot_b_p_spherical + (CFIE(1) * -I * w * mu_0 * mu_r) * I_kk_dot_n_hat_X_b_p_spherical + (CFIE(2) * I * k) * k_hat_X_b_p_spherical + (CFIE(3) * I * k) * k_hat_X_n_hat_X_b_p_spherical;
//               B(local_number_edge_p, 0, index_theta, index_phi) += static_cast< complex<T> > (resultTmp(0));
//               B(local_number_edge_p, 1, index_theta, index_phi) += static_cast< complex<T> > (resultTmp(1));
//             } // end if (IS_SRC==0)
// 	  } // end if (local_number_edge_p > -1)
// 	} // end if (triangles_edges_kinds(index_T, p) > 1)
//       } // end for (p...)
// 
//   } // end for (r...)
// }

template <typename T>
void IT_theta_IT_phi_B_EJ_B_HJ (Array<complex<T>, 4>& B, 
                                const Array<complex<double>, 1>& CFIE, 
                                const Array<int, 1>& indexes_triangles, 
                                const Array<int, 1>& edges_numbers_local_edges_numbers, 
                                const Array<double, 2>& vertexes_coord, 
                                const Array<int, 2>& triangles_vertexes, 
                                const Array<int, 2>& triangles_edges_kinds, 
                                const Array<int, 2>& triangles_edges_numbers, 
                                const Array<int, 2>& triangles_edges_signs, 
                                const Array<double, 2>& triangles_edges_lengths, 
                                const Array<int, 2>& triangles_edges_opp_vertexes, 
                                const Array<double, 2>& triangles_normals, 
                                const Array<double, 2>& triangles_areas, 
                                const Array<double, 2>& edges_numbers_cubes_centroids, 
                                const int IS_SRC, 
                                const complex<double>& k,
                                const double w,
                                const complex<double>& mu_r,
                                const Array<T, 1>& Xtheta, //XcosTheta,
                                const Array<T, 1>& Xphi,
                                const int N_points) /**< INPUT: the number of points used in Gaussian integration */
{
  B = 0.0; // initialization
  for (int index_theta=0 ; index_theta<Xtheta.size() ; index_theta++) {
    cout << "\r" << index_theta*100/Xtheta.size() << " \% of B_EJ_B_HJ completed ";
    flush(cout);
    for (int index_phi=0 ; index_phi<Xphi.size() ; index_phi++) {
      double theta = static_cast<double>(Xtheta(index_theta)), phi = static_cast<double>(Xphi(index_phi));
      B_EJ_B_HJ (B, CFIE, indexes_triangles, edges_numbers_local_edges_numbers, vertexes_coord, triangles_vertexes, triangles_edges_kinds, triangles_edges_numbers, triangles_edges_signs, triangles_edges_lengths, triangles_edges_opp_vertexes, triangles_normals, triangles_areas, edges_numbers_cubes_centroids, IS_SRC, k, w, mu_r, theta, phi, index_theta, index_phi, N_points);
    }
  }
  cout << endl;
}

/***********************************************************/
/*********************** S computation *********************/
/***********************************************************/
template <typename T>
void S_EH_J_computation (Array<complex<T>, 4>& S_EH_J, 
                         const Array<complex<T>, 4>& B_EH_src, 
                         const Array<complex<T>, 1>& I_PQ, 
                         const Array<int, 2>& cubes_edges_numbers, 
                         const Array<int, 1>& cubes_Ei) 
{
  Range all = Range::all();
  int Ej, C = cubes_Ei.size();
  S_EH_J = 0.0;
  for (int j=0 ; j<C ; j++) {
    Ej = cubes_Ei(j);
    Array<int, 1> edges_numbers_src(Ej);
    edges_numbers_src = cubes_edges_numbers(j, Range(0, Ej-1));
    for (int i=0 ; i<Ej ; i++) S_EH_J(j, all, all, all) += B_EH_src(edges_numbers_src(i), all, all, all) * I_PQ(edges_numbers_src(i));
  }
}


/***********************************************************/
/******************** matvec computation *******************/
/***********************************************************/

template <typename T>
void matvec_Z_PQ_near (Array<complex<T>, 1>& ZI_PQ,
                       const Array<complex<T>, 1>& Z_PQ_near,
                       const Array<complex<T>, 1>& I_PQ,
                       const Array<int, 2>& pq_array)
/**
 * matrix-vector multiplication for a sparse matrix data structure such as the coordinate scheme
 * (i, j, a_ij)
 */
{
  int i, p, q;
  int N_near = Z_PQ_near.size();
  for (i=0 ; i<N_near ; i++) {
    p = pq_array(i, 0);
    q = pq_array(i, 1);
    ZI_PQ(p) += Z_PQ_near(i) * I_PQ(q);
  }
}

template <typename T>
void matvec_Z_PQ_near (Array<complex<T>, 1>& ZI_PQ,
                       const Array<complex<T>, 1>& Z_PQ_near,
                       const Array<complex<T>, 1>& I_PQ,
                       const Array<int, 1>& q_array,
                       const Array<int, 2>& rowIndexToColumnIndexes,
                       const Array<int, 1>& RWG_numbers)
/**
 * matrix-vector multiplication for a sparse matrix data structure
 * such as the compressed row storage scheme
 */
{
  const int N_RWG = rowIndexToColumnIndexes.extent(0), N_near = Z_PQ_near.size();
  for (int i=0 ; i<N_RWG ; i++) {
    const int RWG_number = RWG_numbers(i);
    const int startIndex = rowIndexToColumnIndexes(i,0), stopIndex = rowIndexToColumnIndexes(i,1);
    for (int j=startIndex ; j<=stopIndex ; j++) {
      ZI_PQ(RWG_number) += Z_PQ_near(j) * I_PQ(q_array(j));
    }
  }
}

template <typename T>
void matvec_Z_PQ_far_Ci_Cj (Array<complex<T>, 1>& ZI_PQ, 
                            const Array<complex<T>, 4>& B_CFIE_test, 
                            const Array<complex<T>, 3>& S_EH_J_cube_j,
                            const Array<int, 1>& edges_numbers_test)
{
  Range all = Range::all();
  int N_edges_test = edges_numbers_test.size();
  for (int j=0 ; j<N_edges_test ; j++) ZI_PQ(edges_numbers_test(j)) += sum(B_CFIE_test(edges_numbers_test(j), all, all, all) * S_EH_J_cube_j);
}

template <typename T>
void matvec_Z_PQ_far(Array<complex<T>, 1>& ZI_PQ_far, 
                     const Array<complex<T>, 1>& I_PQ, 
                     const Array<complex<T>, 4>& B_CFIE_test, 
                     const Array<complex<T>, 4>& B_tEJ_src, 
                     const Array<T, 2>& Weights2D,
                     const Array<complex<T>, 5>& alpha, 
                     const Array<int, 2>& cubes_edges_numbers, 
                     const Array<int, 1>& cubes_Ei, 
                     const Array<double, 2>& cubes_centroids, 
                     const Array<double, 2>& edges_numbers_cubes_centroids, 
                     const double a, 
                     const int L) 
{
  Range all = Range::all();
  int N_edges;
  const int C = cubes_Ei.size(), N_coord = B_tEJ_src.extent(1), N_theta = B_tEJ_src.extent(2), N_phi = B_tEJ_src.extent(3);
  const int N_Cx = (alpha.extent(0)+1)/2, N_Cy = (alpha.extent(1)+1)/2, N_Cz = (alpha.extent(2)+1)/2;
  double norm_r_ij;
  Array<double, 1> r_ij(3);
  Array<int, 1> mnp(3), offset_mnp(3); // cartesian 3D indexes
  offset_mnp = N_Cx-1, N_Cy-1, N_Cz-1;
  Array<complex<T>, 3> sum_S_EH_J(N_coord, N_theta, N_phi);
  Array<complex<T>, 4> S_EH_J(C, N_coord, N_theta, N_phi);
  S_EH_J_computation(S_EH_J, B_tEJ_src, I_PQ, cubes_edges_numbers, cubes_Ei);

  for (int i=0 ; i<C ; i++) 
  {
    sum_S_EH_J = 0.0;
    N_edges = cubes_Ei(i);
    Array<int, 1> edges_numbers_test( cubes_edges_numbers(i, Range(0, N_edges-1)) );
    for (int j=0 ; j<C ; j++) 
    {
      r_ij = cubes_centroids(i, all) - cubes_centroids(j, all);
      // never EVER change the 2 following lines!!
      mnp = static_cast<int>(round(r_ij(0)/a)), static_cast<int>(round(r_ij(1)/a)), static_cast<int>(round(r_ij(2)/a));
      mnp += offset_mnp;
      norm_r_ij = sqrt(sum(r_ij * r_ij));
      if ( norm_r_ij > sqrt(3.0)*a + 1.0e-5*a ) { // if cubes are not neighbors
        for (int l=0 ; l<N_coord ; l++) sum_S_EH_J(l, all, all) += alpha(mnp(0), mnp(1), mnp(2), all, all) * S_EH_J(j, l, all, all);
      }
    }
    for (int l=0 ; l<N_coord ; l++) sum_S_EH_J(l, all, all) *= Weights2D;
    matvec_Z_PQ_far_Ci_Cj(ZI_PQ_far, B_CFIE_test, sum_S_EH_J, edges_numbers_test);
  }
}

#endif
