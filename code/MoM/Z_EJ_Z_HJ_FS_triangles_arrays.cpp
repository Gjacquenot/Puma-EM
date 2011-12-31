#include <iostream>
#include <iomanip>
#include <complex>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <vector>
#include <algorithm>

using namespace blitz;

#include "EMConstants.h"
#include "triangle_int.h"
#include "mesh.h"

inline blitz::TinyVector<std::complex<double>, 3> cross_real_complex (const blitz::TinyVector<double, 3>& x, const blitz::TinyVector<std::complex<double>, 3>& y) {
  blitz::TinyVector<std::complex<double>, 3> z;
  z = x(1) * y(2) - y(1) * x(2),
      y(0) * x(2) - x(0) * y(2),
      x(0) * y(1) - y(0) * x(1);
  return z;
}

void Z_CFIE_J_computation (blitz::Array<std::complex<double>, 2>& Z_CFIE_J,
                           blitz::Array<std::complex<double>, 2>& Z_CFIE_M,
                           const blitz::Array<std::complex<double>, 1>& CFIE,
                           const double signSurfObs,
                           const double signSurfSrc,
                           const blitz::Array<int, 1>& numbers_RWG_test,
                           const blitz::Array<int, 1>& numbers_RWG_src,
                           const blitz::Array<int, 1>& testRWGNumber_CFIE_OK,
                           const blitz::Array<int, 1>& srcRWGNumber_CURRENT_M_OK,
                           const blitz::Array<int, 2>& RWGNumber_signedTriangles,
                           const blitz::Array<int, 2>& RWGNumber_nodes,
                           const blitz::Array<double, 2>& nodesCoord,
                           const double w,
                           const std::complex<double>& eps_r,
                           const std::complex<double>& mu_r,
                           const int TDS_APPROX, // do we compute for a surface impedance?
                           const std::complex<double>& Z_s, // surface impedance
                           const int FULL_PRECISION)
{
  blitz::Range all = blitz::Range::all();
  const int N_RWG_src = numbers_RWG_src.size(), N_RWG_test = numbers_RWG_test.size();
  // half RWGs construction
  std::vector<RWG> src_RWGs, test_RWGs;
  src_RWGs.reserve(N_RWG_src);
  test_RWGs.reserve(N_RWG_test);
  for (int i=0 ; i<N_RWG_src ; ++i) {
    const int RWGnumber = numbers_RWG_src(i);
    blitz::TinyVector<int, 2> triangle_numbers, triangle_signs;
    triangle_numbers(0) = abs(RWGNumber_signedTriangles(RWGnumber, 0));
    triangle_numbers(1) = abs(RWGNumber_signedTriangles(RWGnumber, 1));
    triangle_signs(0) = 1;
    triangle_signs(1) = -1;
    // the nodes of the RWG
    blitz::Array<double, 1> r0(3), r1(3), r2(3), r3(3);
    r0 = nodesCoord(RWGNumber_nodes(RWGnumber, 0), all);
    r1 = nodesCoord(RWGNumber_nodes(RWGnumber, 1), all);
    r2 = nodesCoord(RWGNumber_nodes(RWGnumber, 2), all);
    r3 = nodesCoord(RWGNumber_nodes(RWGnumber, 3), all);
    src_RWGs.push_back(RWG(RWGnumber, triangle_numbers, triangle_signs, r0, r1, r2, r3));
  }
  for (int i=0 ; i<N_RWG_test ; ++i) {
    const int RWGnumber = numbers_RWG_test(i);
    blitz::TinyVector<int, 2> triangle_numbers, triangle_signs;
    triangle_numbers(0) = abs(RWGNumber_signedTriangles(RWGnumber, 0));
    triangle_numbers(1) = abs(RWGNumber_signedTriangles(RWGnumber, 1));
    triangle_signs(0) = 1;
    triangle_signs(1) = -1;
    // the nodes of the RWG
    blitz::Array<double, 1> r0(3), r1(3), r2(3), r3(3);
    r0 = nodesCoord(RWGNumber_nodes(RWGnumber, 0), all);
    r1 = nodesCoord(RWGNumber_nodes(RWGnumber, 1), all);
    r2 = nodesCoord(RWGNumber_nodes(RWGnumber, 2), all);
    r3 = nodesCoord(RWGNumber_nodes(RWGnumber, 3), all);
    test_RWGs.push_back(RWG(RWGnumber, triangle_numbers, triangle_signs, r0, r1, r2, r3));
  }
  // triangles
  std::vector< Dictionary<int, int> > srcTriangleToRWG, testTriangleToRWG;
  srcTriangleToRWG.reserve(N_RWG_src*2), testTriangleToRWG.reserve(N_RWG_test*2);
  for (int i=0 ; i<src_RWGs.size() ; ++i) {
    srcTriangleToRWG.push_back(Dictionary<int, int> (src_RWGs[i].triangleNumbers(0), src_RWGs[i].number));
    srcTriangleToRWG.push_back(Dictionary<int, int> (src_RWGs[i].triangleNumbers(1), src_RWGs[i].number));
  }
  for (int i=0 ; i<test_RWGs.size() ; ++i) {
    testTriangleToRWG.push_back(Dictionary<int, int> (test_RWGs[i].triangleNumbers(0), test_RWGs[i].number));
    testTriangleToRWG.push_back(Dictionary<int, int> (test_RWGs[i].triangleNumbers(1), test_RWGs[i].number));
  }
  sort(srcTriangleToRWG.begin(), srcTriangleToRWG.end());
  sort(testTriangleToRWG.begin(), testTriangleToRWG.end());
  std::vector<Triangle> triangles_src, triangles_test;
  constructVectorTriangles(triangles_src, src_RWGs, srcTriangleToRWG);
  constructVectorTriangles(triangles_test, test_RWGs, testTriangleToRWG);

  // Z_CFIE computation
  // def of k, mu_i, eps_i
  const bool tds_approx = ((TDS_APPROX != 0) && (Z_s != 0.0));
  const bool tE_tmp = (CFIE(0)!=0.0), nE_tmp = ((CFIE(1)!=0.0) && (!tds_approx)), tH_tmp = ((CFIE(2)!=0.0) && (!tds_approx)), nH_tmp = ((CFIE(3)!=0.0) && (!tds_approx));
  bool IS_TOUCH, IS_NEAR;
  int N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R;
  complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu);

  // declaration of all the geometrical vectors
  blitz::TinyVector<double, 3> r_oc, r_sc;

  // declaration of the scalars and vectors needed in the integrations
  double IT_r_square, n_hat_X_r_p_dot_IT_r;
  blitz::TinyVector<double, 3> IT_r, IT_n_hat_X_r;

  std::complex<double> ITo_ITs_G, ITo_r_dot_ITs_G_rprime, p_dot_ITo_ITs_G_rprime, ITo_n_hat_X_r_dot_ITs_G_rprime, n_hat_X_r_p_dot_ITo_ITs_G_rprime, IDTo_l_hat_dot_r_ITs_G;
  std::complex<double> a_pq, phi_pq, z_pq;
  blitz::TinyVector<std::complex<double>, 3> ITo_r_ITs_G, ITo_ITs_G_rprime, ITo_n_hat_X_r_ITs_G, IDTo_l_hat_ITs_G;

  std::complex<double> ITo_n_hat_X_r_dot_r_X_ITs_grad_G, n_hat_X_r_p_dot_ITo_r_X_ITs_grad_G, n_hat_dot_ITo_r_X_ITs_grad_G;
  blitz::TinyVector<std::complex<double>, 3> ITo_ITs_grad_G, ITo_r_X_ITs_grad_G, ITo_n_hat_X_r_X_ITs_grad_G, r_p_X_ITo_ITs_grad_G, n_hat_X_r_p_X_ITo_ITs_grad_G;

  Z_CFIE_J = 0.0;
  for (int r=0 ; r<triangles_test.size() ; ++r) { // loop on the observation RWGs
    // computation of the triangle-to-triangle terms
    blitz::TinyVector<double, 3> n_hat(triangles_test[r].n_hat);
    IT_fm_fn (IT_r_square, IT_r, triangles_test[r]); // serves for <f_m ; f_n> and <f_m ; n x f_n>
    IT_n_hat_X_r = cross (n_hat, IT_r); // serves for <f_m ; n x f_n>
    // the RWGs concerned by the test triangle
    std::vector<int> RWGsIndexes_test(triangles_test[r].RWGIndexes);
    std::vector<int> triangleTest_indexesInRWGs(triangles_test[r].indexesInRWGs);
    std::vector<double> triangleTest_signsInRWGs(triangles_test[r].signInRWG);
    // we now start the loop on the src triangles
    for (int s=0 ; s<triangles_src.size() ; ++s) { // loop on the source RWGs
      // the RWGs concerned by the source triangle
      std::vector<int> RWGsIndexes_src(triangles_src[s].RWGIndexes);
      std::vector<int> triangleSrc_indexesInRWGs(triangles_src[s].indexesInRWGs);
      std::vector<double> triangleSrc_signsInRWGs(triangles_src[s].signInRWG);

      double R_os = sqrt (dot (triangles_test[r].r_grav - triangles_src[s].r_grav, triangles_test[r].r_grav - triangles_src[s].r_grav));
      IS_TOUCH = ((R_os - triangles_test[r].R_max - triangles_src[s].R_max) <= 0.0);
      IS_NEAR = ((R_os - 1.5 * triangles_test[r].R_max - 1.5 * triangles_src[s].R_max) <= 0.0);
      const bool IS_SAME_TR = (triangles_test[r].number == triangles_src[s].number);

      if (FULL_PRECISION!=0) {
        EXTRACT_1_R = (EXTRACT_R = 0); N_points_o = (N_points_s = 6);
        if ( IS_SAME_TR || IS_TOUCH) {
          EXTRACT_1_R = (EXTRACT_R = 1); N_points_o = (N_points_s = 13);
        }
        else if (IS_NEAR) {
          EXTRACT_1_R = (EXTRACT_R = 1); N_points_o = (N_points_s = 9);
        }
      }
      else {
        EXTRACT_1_R = (EXTRACT_R = 0); N_points_o = (N_points_s = 3);
        if (IS_SAME_TR) {
          EXTRACT_1_R = (EXTRACT_R = 1); N_points_o = (N_points_s = 9);
        }
        else if (IS_TOUCH) {
          EXTRACT_1_R = (EXTRACT_R = 1); N_points_o = (N_points_s = 6);
        }
        else if (IS_NEAR) {
          EXTRACT_1_R = (EXTRACT_R = 0); N_points_o = (N_points_s = 6);
        }
      }

      ITo_ITs_free (ITo_ITs_G, ITo_r_ITs_G, ITo_ITs_G_rprime, ITo_r_dot_ITs_G_rprime, ITo_n_hat_X_r_ITs_G, ITo_n_hat_X_r_dot_ITs_G_rprime, ITo_ITs_grad_G, ITo_r_X_ITs_grad_G, ITo_n_hat_X_r_dot_r_X_ITs_grad_G, ITo_n_hat_X_r_X_ITs_grad_G, triangles_test[r], triangles_src[s], k, N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R);

      n_hat_dot_ITo_r_X_ITs_grad_G = dot (n_hat, ITo_r_X_ITs_grad_G);

      if ((IS_TOUCH) && (nE_tmp || nH_tmp)) IDTo_ITs_free(IDTo_l_hat_dot_r_ITs_G, IDTo_l_hat_ITs_G, triangles_test[r], triangles_src[s], k, 3, N_points_s, EXTRACT_1_R, EXTRACT_R);

      for (int p=0 ; p<RWGsIndexes_test.size() ; ++p) {
        const int index_p = RWGsIndexes_test[p];
        const int local_number_edge_p = test_RWGs[index_p].number;
        const double l_p = test_RWGs[index_p].length;
        const double sign_edge_p = triangleTest_signsInRWGs[p];
        const double C_p = sign_edge_p * l_p * 0.25/triangles_test[r].A;
        blitz::TinyVector<double, 3> r_p;
        if (triangleTest_indexesInRWGs[p]==0) r_p = test_RWGs[index_p].vertexesCoord(0);
        else if (triangleTest_indexesInRWGs[p]==1) r_p = test_RWGs[index_p].vertexesCoord(3);
        const blitz::TinyVector<double, 3> n_hat_X_r_p(cross(n_hat, r_p));
        const int IS_CFIE = testRWGNumber_CFIE_OK(local_number_edge_p);
        const bool tEJ = tE_tmp, nEJ = (nE_tmp * IS_CFIE) , tHJ = (tH_tmp * IS_CFIE ), nHJ = (nH_tmp * IS_CFIE);

        // temporary elements for Z_tE_J
        p_dot_ITo_ITs_G_rprime = dot (r_p, ITo_ITs_G_rprime);

        // temporary elements for Z_nE_J
        n_hat_X_r_p_dot_ITo_ITs_G_rprime = dot(n_hat_X_r_p, ITo_ITs_G_rprime);

        // temporary elements for Z_tH_J
        n_hat_X_r_p_dot_IT_r = dot(n_hat_X_r_p, IT_r);
        r_p_X_ITo_ITs_grad_G = cross_real_complex(r_p, ITo_ITs_grad_G);

        // temporary elements for Z_nH_J
        n_hat_X_r_p_X_ITo_ITs_grad_G = cross_real_complex (n_hat_X_r_p, ITo_ITs_grad_G);
        n_hat_X_r_p_dot_ITo_r_X_ITs_grad_G = dot (n_hat_X_r_p, ITo_r_X_ITs_grad_G);

        for (int q=0 ; q<RWGsIndexes_src.size() ; ++q) {
          const int index_q = RWGsIndexes_src[q];
          const int local_number_edge_q = src_RWGs[index_q].number;
          const double l_q = src_RWGs[index_q].length;
          const double sign_edge_q = triangleSrc_signsInRWGs[q];
          const double C_pq = sign_edge_q * l_q * C_p/triangles_src[s].A;
          blitz::TinyVector<double, 3> r_q;
          if (triangleSrc_indexesInRWGs[q]==0) r_q = src_RWGs[index_q].vertexesCoord(0);
          else if (triangleSrc_indexesInRWGs[q]==1) r_q = src_RWGs[index_q].vertexesCoord(3);
          const bool M_CURRENT_OK = (srcRWGNumber_CURRENT_M_OK(local_number_edge_q)==1);
          const bool tHM = (tH_tmp && M_CURRENT_OK), nHM = (nH_tmp && M_CURRENT_OK), tEM = (tE_tmp && M_CURRENT_OK), nEM = (nE_tmp && M_CURRENT_OK);

          // <f_p ; EFIE> : Z_tE_J computation. Z_tH_M = eps/mu * Z_tE_J
          if (tEJ || tHM) {
            a_pq = C_pq * (ITo_r_dot_ITs_G_rprime - dot(ITo_r_ITs_G, r_q) - p_dot_ITo_ITs_G_rprime + dot(r_p, r_q) * ITo_ITs_G);
            phi_pq = 4.0 * C_pq * ITo_ITs_G;
            z_pq = -I * w * mu * (a_pq - 1.0/(k*k) * phi_pq);
            if (tds_approx && IS_SAME_TR) z_pq += 0.5 * C_pq * Z_s * ( IT_r_square - dot((r_p + r_q), IT_r) + dot(r_p, r_q) * triangles_test[r].A );
            if (tEJ) Z_CFIE_J (local_number_edge_p, local_number_edge_q) += signSurfObs * signSurfSrc * CFIE(0) * z_pq;
            if (tHM) Z_CFIE_M (local_number_edge_p, local_number_edge_q) += signSurfObs * signSurfSrc * CFIE(2) * eps/mu * z_pq;
          }
          // <n x f_p ; EFIE> : Z_nE_J computation. Z_nH_M = eps/mu * Z_nE_J
          if (nEJ || nHM) {
            a_pq = C_pq * (ITo_n_hat_X_r_dot_ITs_G_rprime - dot(ITo_n_hat_X_r_ITs_G, r_q) - n_hat_X_r_p_dot_ITo_ITs_G_rprime + dot (n_hat_X_r_p, r_q) * ITo_ITs_G);
            if (IS_TOUCH) phi_pq = 2.0 * C_pq * (-IDTo_l_hat_dot_r_ITs_G + dot(r_p, IDTo_l_hat_ITs_G));
            else phi_pq = 2.0 * C_pq * ( n_hat_dot_ITo_r_X_ITs_grad_G - dot(n_hat_X_r_p, ITo_ITs_grad_G));
            z_pq = -I * w * mu * (a_pq + 1.0/(k*k) * phi_pq);
            if (nEJ) Z_CFIE_J (local_number_edge_p, local_number_edge_q) += signSurfObs * signSurfSrc * CFIE(1) * z_pq;
            if (nHM) Z_CFIE_M (local_number_edge_p, local_number_edge_q) += signSurfObs * signSurfSrc * CFIE(3) * eps/mu * z_pq;
          }
          // <f_p ; MFIE> : Z_tH_J computation. Z_tE_M = -Z_tH_J
          if (tHJ || tEM) {
            if (!IS_SAME_TR) z_pq = signSurfObs * signSurfSrc * C_pq * dot (r_p-r_q, ITo_r_X_ITs_grad_G - r_p_X_ITo_ITs_grad_G);
            // else we have: z_pq = 0.5 * <f_p ; n x f_q>
            else z_pq = signSurfObs * 0.5 * C_pq * ( dot(n_hat_X_r_p, IT_r) + dot(r_q, IT_n_hat_X_r) - dot(r_q, n_hat_X_r_p) * triangles_test[r].A );
            if (tHJ) Z_CFIE_J (local_number_edge_p, local_number_edge_q) += CFIE(2) * z_pq;
            if (tEM) Z_CFIE_M (local_number_edge_p, local_number_edge_q) -= CFIE(0) * z_pq;
          }
          // <n x f_p ; MFIE> : Z_nH_J computation. Z_nE_M = -Z_nH_J
          if (nHJ || nEM) {
            if (!IS_SAME_TR) z_pq = -signSurfObs * signSurfSrc * C_pq * ( ITo_n_hat_X_r_dot_r_X_ITs_grad_G + dot (r_q, ITo_n_hat_X_r_X_ITs_grad_G - n_hat_X_r_p_X_ITo_ITs_grad_G) - n_hat_X_r_p_dot_ITo_r_X_ITs_grad_G  );
            // else we have: z_pq = 0.5 * <n x f_p ; n x f_q> = 0.5 * <f_p ; f_q>
            else z_pq = signSurfObs * 0.5 * C_pq * ( IT_r_square - dot((r_p + r_q), IT_r) + dot(r_p, r_q) * triangles_test[r].A );
            if (nHJ) Z_CFIE_J (local_number_edge_p, local_number_edge_q) += CFIE(3) * z_pq;
            if (nEM) Z_CFIE_M (local_number_edge_p, local_number_edge_q) -= CFIE(1) * z_pq;
          }

        } // for q
      } // for p
      
    } // for s
  } // for r
}












