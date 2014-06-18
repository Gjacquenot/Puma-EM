#include <iostream>
#include <iomanip>
#include <complex>
#include <blitz/array.h>
#include <vector>
#include <algorithm>

using namespace std;

#include "EMConstants.h"
#include "triangle_int.h"
#include "dictionary.h"

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
  const int N_RWG_src = numbers_RWG_src.size(), N_RWG_test = numbers_RWG_test.size();
  // half RWGs construction
  std::vector<RWG> src_RWGs, test_RWGs;
  src_RWGs.reserve(N_RWG_src);
  test_RWGs.reserve(N_RWG_test);
  for (int i=0 ; i<N_RWG_src ; ++i) {
    const int RWGnumber = numbers_RWG_src(i);
    int triangle_numbers[2], triangle_signs[2];
    triangle_numbers[0] = abs(RWGNumber_signedTriangles(RWGnumber, 0));
    triangle_numbers[1] = abs(RWGNumber_signedTriangles(RWGnumber, 1));
    triangle_signs[0] = 1;
    triangle_signs[1] = -1;
    // the nodes of the RWG
    const int n0 = RWGNumber_nodes(RWGnumber, 0);
    const int n1 = RWGNumber_nodes(RWGnumber, 1);
    const int n2 = RWGNumber_nodes(RWGnumber, 2);
    const int n3 = RWGNumber_nodes(RWGnumber, 3);
    const double r0[3] = { nodesCoord(n0, 0), nodesCoord(n0, 1), nodesCoord(n0, 2)};
    const double r1[3] = { nodesCoord(n1, 0), nodesCoord(n1, 1), nodesCoord(n1, 2)};
    const double r2[3] = { nodesCoord(n2, 0), nodesCoord(n2, 1), nodesCoord(n2, 2)};
    const double r3[3] = { nodesCoord(n3, 0), nodesCoord(n3, 1), nodesCoord(n3, 2)};
    src_RWGs.push_back(RWG(RWGnumber, triangle_numbers, triangle_signs, r0, r1, r2, r3));
  }
  for (int i=0 ; i<N_RWG_test ; ++i) {
    const int RWGnumber = numbers_RWG_test(i);
    int triangle_numbers[2], triangle_signs[2];
    triangle_numbers[0] = abs(RWGNumber_signedTriangles(RWGnumber, 0));
    triangle_numbers[1] = abs(RWGNumber_signedTriangles(RWGnumber, 1));
    triangle_signs[0] = 1;
    triangle_signs[1] = -1;
    // the nodes of the RWG
    const int n0 = RWGNumber_nodes(RWGnumber, 0);
    const int n1 = RWGNumber_nodes(RWGnumber, 1);
    const int n2 = RWGNumber_nodes(RWGnumber, 2);
    const int n3 = RWGNumber_nodes(RWGnumber, 3);
    const double r0[3] = { nodesCoord(n0, 0), nodesCoord(n0, 1), nodesCoord(n0, 2)};
    const double r1[3] = { nodesCoord(n1, 0), nodesCoord(n1, 1), nodesCoord(n1, 2)};
    const double r2[3] = { nodesCoord(n2, 0), nodesCoord(n2, 1), nodesCoord(n2, 2)};
    const double r3[3] = { nodesCoord(n3, 0), nodesCoord(n3, 1), nodesCoord(n3, 2)};
    test_RWGs.push_back(RWG(RWGnumber, triangle_numbers, triangle_signs, r0, r1, r2, r3));
  }
  // triangles
  std::vector< Dictionary<int, int> > srcTriangleToRWG, testTriangleToRWG;
  srcTriangleToRWG.reserve(N_RWG_src*2), testTriangleToRWG.reserve(N_RWG_test*2);
  for (unsigned int i=0 ; i<src_RWGs.size() ; ++i) {
    srcTriangleToRWG.push_back(Dictionary<int, int> (src_RWGs[i].triangleNumbers[0], src_RWGs[i].number));
    srcTriangleToRWG.push_back(Dictionary<int, int> (src_RWGs[i].triangleNumbers[1], src_RWGs[i].number));
  }
  for (unsigned int i=0 ; i<test_RWGs.size() ; ++i) {
    testTriangleToRWG.push_back(Dictionary<int, int> (test_RWGs[i].triangleNumbers[0], test_RWGs[i].number));
    testTriangleToRWG.push_back(Dictionary<int, int> (test_RWGs[i].triangleNumbers[1], test_RWGs[i].number));
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
  int N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R;
  const std::complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu), k_square = k*k;
  const std::complex<double> z_pq_factor = 1.0/(I*w*eps);

  Z_CFIE_J = 0.0; Z_CFIE_M = 0.0;
  for (unsigned int r=0 ; r<triangles_test.size() ; ++r) { // loop on the observation RWGs
    // computation of the triangle-to-triangle terms
    double *n_hat;
    n_hat = triangles_test[r].n_hat;
    double IT_r_square; // n_hat_X_r_p_dot_IT_r;
    double IT_r[3], IT_n_hat_X_r[3];
    IT_fm_fn (IT_r_square, IT_r, triangles_test[r]); // serves for <f_m ; f_n> and <f_m ; n x f_n>
    cross3D(IT_n_hat_X_r, n_hat, IT_r); // serves for <f_m ; n x f_n>
    // the RWGs concerned by the test triangle
    std::vector<int> RWGsIndexes_test(triangles_test[r].RWGIndexes);
    std::vector<int> triangleTest_indexesInRWGs(triangles_test[r].indexesInRWGs);
    std::vector<double> triangleTest_signsInRWGs(triangles_test[r].signInRWG);
    // we now start the loop on the src triangles
    for (unsigned int s=0 ; s<triangles_src.size() ; ++s) { // loop on the source RWGs
      // the RWGs concerned by the source triangle
      std::vector<int> RWGsIndexes_src(triangles_src[s].RWGIndexes);
      std::vector<int> triangleSrc_indexesInRWGs(triangles_src[s].indexesInRWGs);
      std::vector<double> triangleSrc_signsInRWGs(triangles_src[s].signInRWG);

      double r_grav_obs_r_grav_src[3] = {triangles_test[r].r_grav[0] - triangles_src[s].r_grav[0], triangles_test[r].r_grav[1] - triangles_src[s].r_grav[1], triangles_test[r].r_grav[2] - triangles_src[s].r_grav[2]};
      double R_os = sqrt(dot3D(r_grav_obs_r_grav_src, r_grav_obs_r_grav_src));
      const bool IS_TOUCH = ((R_os - triangles_test[r].R_max - triangles_src[s].R_max) <= 0.0);
      const bool IS_NEAR = ((R_os - 1.5 * triangles_test[r].R_max - 1.5 * triangles_src[s].R_max) <= 0.0);
      const bool IS_SAME_TR = (triangles_test[r].number == triangles_src[s].number);

      if (FULL_PRECISION!=0) {
        EXTRACT_1_R = (EXTRACT_R = 0); N_points_o = (N_points_s = 6);
        if (IS_SAME_TR) {
          EXTRACT_1_R = (EXTRACT_R = 1); N_points_o = (N_points_s = 13);
        }
        else if (IS_TOUCH) {
          EXTRACT_1_R = (EXTRACT_R = 1); N_points_o = (N_points_s = 9);
        }
        else if (IS_NEAR) {
          EXTRACT_1_R = (EXTRACT_R = 1); N_points_o = (N_points_s = 6);
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

      // declaration of the scalars and vectors needed in the integrations
      std::complex<double> ITo_ITs_G, ITo_r_dot_ITs_G_rprime, ITo_n_hat_X_r_dot_ITs_G_rprime, IDTo_l_hat_dot_r_ITs_G;
      std::complex<double> ITo_n_hat_X_r_dot_r_X_ITs_grad_G;
      std::complex<double> ITo_r_ITs_G[3], ITo_ITs_G_rprime[3], IDTo_l_hat_ITs_G[3];
      std::complex<double> ITo_ITs_grad_G[3], ITo_r_X_ITs_grad_G[3], ITo_n_hat_X_r_X_ITs_grad_G[3];

      ITo_ITs_free (ITo_ITs_G, ITo_r_ITs_G, ITo_ITs_G_rprime, ITo_r_dot_ITs_G_rprime, ITo_n_hat_X_r_dot_ITs_G_rprime, ITo_ITs_grad_G, ITo_r_X_ITs_grad_G, ITo_n_hat_X_r_dot_r_X_ITs_grad_G, ITo_n_hat_X_r_X_ITs_grad_G, triangles_test[r], triangles_src[s], k, N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R);

      const std::complex<double> ITo_n_hat_X_r_ITs_G[3] = {n_hat[1]*ITo_r_ITs_G[2] - n_hat[2]*ITo_r_ITs_G[1],
                                                           n_hat[2]*ITo_r_ITs_G[0] - n_hat[0]*ITo_r_ITs_G[2],
                                                           n_hat[0]*ITo_r_ITs_G[1] - n_hat[1]*ITo_r_ITs_G[0]};
      const std::complex<double> n_hat_dot_ITo_r_X_ITs_grad_G(n_hat[0]*ITo_r_X_ITs_grad_G[0] + n_hat[1]*ITo_r_X_ITs_grad_G[1] + n_hat[2]*ITo_r_X_ITs_grad_G[2]);

      if ((IS_TOUCH) && (nE_tmp || nH_tmp)) IDTo_ITs_free(IDTo_l_hat_dot_r_ITs_G, IDTo_l_hat_ITs_G, triangles_test[r], triangles_src[s], k, 3, N_points_s, EXTRACT_1_R, EXTRACT_R);

      for (unsigned int p=0 ; p<RWGsIndexes_test.size() ; ++p) {
        const int index_p = RWGsIndexes_test[p];
        const int local_number_edge_p = test_RWGs[index_p].number;
        const double l_p = test_RWGs[index_p].length;
        const double sign_edge_p = triangleTest_signsInRWGs[p];
        const double C_p = sign_edge_p * l_p * 0.5/triangles_test[r].A;
        double *r_p, n_hat_X_r_p[3];
        if (triangleTest_indexesInRWGs[p]==0) r_p = test_RWGs[index_p].vertexesCoord_0;
        else r_p = test_RWGs[index_p].vertexesCoord_3;
        cross3D(n_hat_X_r_p, n_hat, r_p);
        const int IS_CFIE = testRWGNumber_CFIE_OK(local_number_edge_p);
        const bool tEJ = tE_tmp, nEJ = nE_tmp, tHJ = (tH_tmp * IS_CFIE ), nHJ = (nH_tmp * IS_CFIE);

        // temporary elements for Z_tE_J
        const std::complex<double> p_dot_ITo_ITs_G_rprime(r_p[0]*ITo_ITs_G_rprime[0] + r_p[1]*ITo_ITs_G_rprime[1] + r_p[2]*ITo_ITs_G_rprime[2]);

        // temporary elements for Z_nE_J
        const std::complex<double> n_hat_X_r_p_dot_ITo_ITs_G_rprime(n_hat_X_r_p[0]*ITo_ITs_G_rprime[0] + n_hat_X_r_p[1]*ITo_ITs_G_rprime[1] + n_hat_X_r_p[2]*ITo_ITs_G_rprime[2]);
        const std::complex<double> n_hat_X_r_p_dot_ITo_ITs_grad_G(n_hat_X_r_p[0]*ITo_ITs_grad_G[0] + n_hat_X_r_p[1]*ITo_ITs_grad_G[1] + n_hat_X_r_p[2]*ITo_ITs_grad_G[2]);

        // temporary elements for Z_tH_J
        const std::complex<double> r_p_X_ITo_ITs_grad_G[3] = {r_p[1]*ITo_ITs_grad_G[2] - r_p[2]*ITo_ITs_grad_G[1],
                                                              r_p[2]*ITo_ITs_grad_G[0] - r_p[0]*ITo_ITs_grad_G[2],
                                                              r_p[0]*ITo_ITs_grad_G[1] - r_p[1]*ITo_ITs_grad_G[0]};
        const std::complex<double> ITo_r_rp_X_ITs_grad_G[3] = {ITo_r_X_ITs_grad_G[0]-r_p_X_ITo_ITs_grad_G[0],
                                                               ITo_r_X_ITs_grad_G[1]-r_p_X_ITo_ITs_grad_G[1],
                                                               ITo_r_X_ITs_grad_G[2]-r_p_X_ITo_ITs_grad_G[2]};
        const double n_hat_X_rp_dot_IT_r(n_hat_X_r_p[0]*IT_r[0] + n_hat_X_r_p[1]*IT_r[1] + n_hat_X_r_p[2]*IT_r[2]);

        // temporary elements for Z_nH_J
        const std::complex<double> n_hat_X_r_p_X_ITo_ITs_grad_G[3] = {n_hat_X_r_p[1]*ITo_ITs_grad_G[2] - n_hat_X_r_p[2]*ITo_ITs_grad_G[1],
                                                                      n_hat_X_r_p[2]*ITo_ITs_grad_G[0] - n_hat_X_r_p[0]*ITo_ITs_grad_G[2],
                                                                      n_hat_X_r_p[0]*ITo_ITs_grad_G[1] - n_hat_X_r_p[1]*ITo_ITs_grad_G[0]};
        const std::complex<double> ITo_n_hat_X_r_rp_X_ITs_grad_G[3] = {ITo_n_hat_X_r_X_ITs_grad_G[0] - n_hat_X_r_p_X_ITo_ITs_grad_G[0], 
                                                                       ITo_n_hat_X_r_X_ITs_grad_G[1] - n_hat_X_r_p_X_ITo_ITs_grad_G[1],
                                                                       ITo_n_hat_X_r_X_ITs_grad_G[2] - n_hat_X_r_p_X_ITo_ITs_grad_G[2]};
        const std::complex<double> n_hat_X_r_p_dot_ITo_r_X_ITs_grad_G(n_hat_X_r_p[0]*ITo_r_X_ITs_grad_G[0] + n_hat_X_r_p[1]*ITo_r_X_ITs_grad_G[1] + n_hat_X_r_p[2]*ITo_r_X_ITs_grad_G[2]);

        for (unsigned int q=0 ; q<RWGsIndexes_src.size() ; ++q) {
          const int index_q = RWGsIndexes_src[q];
          const int local_number_edge_q = src_RWGs[index_q].number;
          const double l_q = src_RWGs[index_q].length;
          const double sign_edge_q = triangleSrc_signsInRWGs[q];
          const double C_pq = C_p * sign_edge_q * l_q * 0.5/triangles_src[s].A;
          double *r_q;
          if (triangleSrc_indexesInRWGs[q]==0) r_q = src_RWGs[index_q].vertexesCoord_0;
          else r_q = src_RWGs[index_q].vertexesCoord_3;
          const bool M_CURRENT_OK = (srcRWGNumber_CURRENT_M_OK(local_number_edge_q)==1);
          const bool tHM = (tH_tmp && M_CURRENT_OK), nHM = (nH_tmp && M_CURRENT_OK), tEM = (tE_tmp && M_CURRENT_OK), nEM = (nE_tmp && M_CURRENT_OK);
          const double rp_dot_rq(r_p[0]*r_q[0] + r_p[1]*r_q[1] + r_p[2]*r_q[2]);
          const double rp_plus_rq_dot_ITr((r_p[0]+r_q[0])*IT_r[0] + (r_p[1]+r_q[1])*IT_r[1] + (r_p[2]+r_q[2])*IT_r[2]);
          std::complex<double> z_pq, D_mn_1, D_mn_2;

          // <f_p ; EFIE> : Z_tE_J computation. Z_tH_M = eps/mu * Z_tE_J
          if (tEJ || tHM) {
            D_mn_1 = (-4.0 * C_pq) * ITo_ITs_G;
            D_mn_2 = C_pq * (ITo_r_dot_ITs_G_rprime - (ITo_r_ITs_G[0]*r_q[0] + ITo_r_ITs_G[1]*r_q[1] + ITo_r_ITs_G[2]*r_q[2]) - p_dot_ITo_ITs_G_rprime + rp_dot_rq * ITo_ITs_G);
            z_pq = z_pq_factor * (D_mn_1 + k_square * D_mn_2); // z_pq_factor = 1.0/(I*w*eps)
            if (tds_approx && IS_SAME_TR) z_pq += 0.5 * C_pq * Z_s * ( IT_r_square - rp_plus_rq_dot_ITr + rp_dot_rq * triangles_test[r].A );
            if (tEJ) Z_CFIE_J (local_number_edge_p, local_number_edge_q) += (signSurfObs * signSurfSrc) * CFIE(0) * z_pq;
            if (tHM) Z_CFIE_M (local_number_edge_p, local_number_edge_q) += (signSurfObs * signSurfSrc) * CFIE(2) * eps/mu * z_pq;
          }
          // <n x f_p ; EFIE> : Z_nE_J computation. Z_nH_M = eps/mu * Z_nE_J
          if (nEJ || nHM) {
            if (IS_TOUCH) D_mn_1 = 2.0 * C_pq * (-IDTo_l_hat_dot_r_ITs_G + (r_p[0]*IDTo_l_hat_ITs_G[0] + r_p[1]*IDTo_l_hat_ITs_G[1] + r_p[2]*IDTo_l_hat_ITs_G[2]));
            else D_mn_1 = 2.0 * C_pq * ( n_hat_dot_ITo_r_X_ITs_grad_G - n_hat_X_r_p_dot_ITo_ITs_grad_G);
            D_mn_2 = C_pq * (ITo_n_hat_X_r_dot_ITs_G_rprime - (ITo_n_hat_X_r_ITs_G[0]*r_q[0] + ITo_n_hat_X_r_ITs_G[1]*r_q[1] + ITo_n_hat_X_r_ITs_G[2]*r_q[2]) - n_hat_X_r_p_dot_ITo_ITs_G_rprime + (n_hat_X_r_p[0]*r_q[0] + n_hat_X_r_p[1]*r_q[1] + n_hat_X_r_p[2]*r_q[2]) * ITo_ITs_G);
            z_pq = z_pq_factor * ( D_mn_1 + k_square*D_mn_2); // z_pq_factor = 1.0/(I*w*eps)
            if (nEJ) Z_CFIE_J (local_number_edge_p, local_number_edge_q) += (signSurfObs * signSurfSrc) * CFIE(1) * z_pq;
            if (nHM) Z_CFIE_M (local_number_edge_p, local_number_edge_q) += (signSurfObs * signSurfSrc) * CFIE(3) * eps/mu * z_pq;
          }
          // <f_p ; MFIE> : Z_tH_J computation. Z_tE_M = -Z_tH_J
          if (tHJ || tEM) {
            if (!IS_SAME_TR) z_pq = (signSurfObs * signSurfSrc * C_pq) * ((r_p[0]-r_q[0])*ITo_r_rp_X_ITs_grad_G[0] + (r_p[1]-r_q[1])*ITo_r_rp_X_ITs_grad_G[1] + (r_p[2]-r_q[2])*ITo_r_rp_X_ITs_grad_G[2]);
            // else we have: z_pq = 0.5 * <f_p ; n x f_q>
            else z_pq = (signSurfObs * 0.5 * C_pq) * ( n_hat_X_rp_dot_IT_r + (r_q[0]*IT_n_hat_X_r[0] + r_q[1]*IT_n_hat_X_r[1] + r_q[2]*IT_n_hat_X_r[2]) - (r_q[0]*n_hat_X_r_p[0] + r_q[1]*n_hat_X_r_p[1] + r_q[2]*n_hat_X_r_p[2]) * triangles_test[r].A );
            if (tHJ) Z_CFIE_J (local_number_edge_p, local_number_edge_q) += CFIE(2) * z_pq;
            if (tEM) Z_CFIE_M (local_number_edge_p, local_number_edge_q) -= CFIE(0) * z_pq;
          }
          // <n x f_p ; MFIE> : Z_nH_J computation. Z_nE_M = -Z_nH_J
          if (nHJ || nEM) {
            if (!IS_SAME_TR) z_pq = (-signSurfObs * signSurfSrc * C_pq) * ( ITo_n_hat_X_r_dot_r_X_ITs_grad_G + (r_q[0]*ITo_n_hat_X_r_rp_X_ITs_grad_G[0] + r_q[1]*ITo_n_hat_X_r_rp_X_ITs_grad_G[1] + r_q[2]*ITo_n_hat_X_r_rp_X_ITs_grad_G[2]) - n_hat_X_r_p_dot_ITo_r_X_ITs_grad_G  );
            // else we have: z_pq = 0.5 * <n x f_p ; n x f_q> = 0.5 * <f_p ; f_q>
            else z_pq = (signSurfObs * 0.5 * C_pq) * ( IT_r_square - rp_plus_rq_dot_ITr + rp_dot_rq * triangles_test[r].A );
            if (nHJ) Z_CFIE_J (local_number_edge_p, local_number_edge_q) += CFIE(3) * z_pq;
            if (nEM) Z_CFIE_M (local_number_edge_p, local_number_edge_q) -= CFIE(1) * z_pq;
          }

        } // for q
      } // for p
      
    } // for s
  } // for r
}


void Z_EH_J_computation (blitz::Array<std::complex<double>, 2>& Z_tE_J,
                         blitz::Array<std::complex<double>, 2>& Z_nE_J,
                         blitz::Array<std::complex<double>, 2>& Z_tH_J,
                         blitz::Array<std::complex<double>, 2>& Z_nH_J,
                         const blitz::Array<double, 1>& TENETHNH,
                         const double signSurfObs,
                         const double signSurfSrc,
                         const blitz::Array<int, 1>& numbers_RWG_test,
                         const blitz::Array<int, 1>& numbers_RWG_src,
                         const blitz::Array<int, 1>& testRWGNumber_CFIE_OK,
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
  const int N_RWG_src = numbers_RWG_src.size(), N_RWG_test = numbers_RWG_test.size();
  // half RWGs construction
  std::vector<RWG> src_RWGs, test_RWGs;
  src_RWGs.reserve(N_RWG_src);
  test_RWGs.reserve(N_RWG_test);
  for (int i=0 ; i<N_RWG_src ; ++i) {
    const int RWGnumber = numbers_RWG_src(i);
    int triangle_numbers[2], triangle_signs[2];
    triangle_numbers[0] = abs(RWGNumber_signedTriangles(RWGnumber, 0));
    triangle_numbers[1] = abs(RWGNumber_signedTriangles(RWGnumber, 1));
    triangle_signs[0] = 1;
    triangle_signs[1] = -1;
    // the nodes of the RWG
    const int n0 = RWGNumber_nodes(RWGnumber, 0);
    const int n1 = RWGNumber_nodes(RWGnumber, 1);
    const int n2 = RWGNumber_nodes(RWGnumber, 2);
    const int n3 = RWGNumber_nodes(RWGnumber, 3);
    const double r0[3] = { nodesCoord(n0, 0), nodesCoord(n0, 1), nodesCoord(n0, 2)};
    const double r1[3] = { nodesCoord(n1, 0), nodesCoord(n1, 1), nodesCoord(n1, 2)};
    const double r2[3] = { nodesCoord(n2, 0), nodesCoord(n2, 1), nodesCoord(n2, 2)};
    const double r3[3] = { nodesCoord(n3, 0), nodesCoord(n3, 1), nodesCoord(n3, 2)};
    src_RWGs.push_back(RWG(RWGnumber, triangle_numbers, triangle_signs, r0, r1, r2, r3));
  }
  for (int i=0 ; i<N_RWG_test ; ++i) {
    const int RWGnumber = numbers_RWG_test(i);
    int triangle_numbers[2], triangle_signs[2];
    triangle_numbers[0] = abs(RWGNumber_signedTriangles(RWGnumber, 0));
    triangle_numbers[1] = abs(RWGNumber_signedTriangles(RWGnumber, 1));
    triangle_signs[0] = 1;
    triangle_signs[1] = -1;
    // the nodes of the RWG
    const int n0 = RWGNumber_nodes(RWGnumber, 0);
    const int n1 = RWGNumber_nodes(RWGnumber, 1);
    const int n2 = RWGNumber_nodes(RWGnumber, 2);
    const int n3 = RWGNumber_nodes(RWGnumber, 3);
    const double r0[3] = { nodesCoord(n0, 0), nodesCoord(n0, 1), nodesCoord(n0, 2)};
    const double r1[3] = { nodesCoord(n1, 0), nodesCoord(n1, 1), nodesCoord(n1, 2)};
    const double r2[3] = { nodesCoord(n2, 0), nodesCoord(n2, 1), nodesCoord(n2, 2)};
    const double r3[3] = { nodesCoord(n3, 0), nodesCoord(n3, 1), nodesCoord(n3, 2)};
    test_RWGs.push_back(RWG(RWGnumber, triangle_numbers, triangle_signs, r0, r1, r2, r3));
  }
  // triangles
  std::vector< Dictionary<int, int> > srcTriangleToRWG, testTriangleToRWG;
  srcTriangleToRWG.reserve(N_RWG_src*2), testTriangleToRWG.reserve(N_RWG_test*2);
  for (unsigned int i=0 ; i<src_RWGs.size() ; ++i) {
    srcTriangleToRWG.push_back(Dictionary<int, int> (src_RWGs[i].triangleNumbers[0], src_RWGs[i].number));
    srcTriangleToRWG.push_back(Dictionary<int, int> (src_RWGs[i].triangleNumbers[1], src_RWGs[i].number));
  }
  for (unsigned int i=0 ; i<test_RWGs.size() ; ++i) {
    testTriangleToRWG.push_back(Dictionary<int, int> (test_RWGs[i].triangleNumbers[0], test_RWGs[i].number));
    testTriangleToRWG.push_back(Dictionary<int, int> (test_RWGs[i].triangleNumbers[1], test_RWGs[i].number));
  }
  sort(srcTriangleToRWG.begin(), srcTriangleToRWG.end());
  sort(testTriangleToRWG.begin(), testTriangleToRWG.end());
  std::vector<Triangle> triangles_src, triangles_test;
  constructVectorTriangles(triangles_src, src_RWGs, srcTriangleToRWG);
  constructVectorTriangles(triangles_test, test_RWGs, testTriangleToRWG);

  // Z_CFIE computation
  // def of k, mu_i, eps_i
  const bool tds_approx = ((TDS_APPROX != 0) && (Z_s != 0.0));
  const bool tE_tmp = (TENETHNH(0)!=0.0), nE_tmp = ((TENETHNH(1)!=0.0) && (!tds_approx)), tH_tmp = ((TENETHNH(2)!=0.0) && (!tds_approx)), nH_tmp = ((TENETHNH(3)!=0.0) && (!tds_approx));
  int N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R;
  const std::complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu), k_square = k*k;
  const std::complex<double> z_pq_factor = 1.0/(I*w*eps);

  Z_tE_J = 0.0; Z_nE_J = 0.0;
  Z_tH_J = 0.0; Z_nH_J = 0.0;
  for (unsigned int r=0 ; r<triangles_test.size() ; ++r) { // loop on the observation RWGs
    // computation of the triangle-to-triangle terms
    double *n_hat;
    n_hat = triangles_test[r].n_hat;
    double IT_r_square; // n_hat_X_r_p_dot_IT_r;
    double IT_r[3], IT_n_hat_X_r[3];
    IT_fm_fn (IT_r_square, IT_r, triangles_test[r]); // serves for <f_m ; f_n> and <f_m ; n x f_n>
    cross3D(IT_n_hat_X_r, n_hat, IT_r); // serves for <f_m ; n x f_n>
    // the RWGs concerned by the test triangle
    std::vector<int> RWGsIndexes_test(triangles_test[r].RWGIndexes);
    std::vector<int> triangleTest_indexesInRWGs(triangles_test[r].indexesInRWGs);
    std::vector<double> triangleTest_signsInRWGs(triangles_test[r].signInRWG);
    // we now start the loop on the src triangles
    for (unsigned int s=0 ; s<triangles_src.size() ; ++s) { // loop on the source RWGs
      // the RWGs concerned by the source triangle
      std::vector<int> RWGsIndexes_src(triangles_src[s].RWGIndexes);
      std::vector<int> triangleSrc_indexesInRWGs(triangles_src[s].indexesInRWGs);
      std::vector<double> triangleSrc_signsInRWGs(triangles_src[s].signInRWG);

      double r_grav_obs_r_grav_src[3] = {triangles_test[r].r_grav[0] - triangles_src[s].r_grav[0], triangles_test[r].r_grav[1] - triangles_src[s].r_grav[1], triangles_test[r].r_grav[2] - triangles_src[s].r_grav[2]};
      double R_os = sqrt(dot3D(r_grav_obs_r_grav_src, r_grav_obs_r_grav_src));
      const bool IS_TOUCH = ((R_os - triangles_test[r].R_max - triangles_src[s].R_max) <= 0.0);
      const bool IS_NEAR = ((R_os - 1.5 * triangles_test[r].R_max - 1.5 * triangles_src[s].R_max) <= 0.0);
      const bool IS_SAME_TR = (triangles_test[r].number == triangles_src[s].number);

      if (FULL_PRECISION!=0) {
        EXTRACT_1_R = (EXTRACT_R = 0); N_points_o = (N_points_s = 6);
        if (IS_SAME_TR) {
          EXTRACT_1_R = (EXTRACT_R = 1); N_points_o = (N_points_s = 13);
        }
        else if (IS_TOUCH) {
          EXTRACT_1_R = (EXTRACT_R = 1); N_points_o = (N_points_s = 9);
        }
        else if (IS_NEAR) {
          EXTRACT_1_R = (EXTRACT_R = 1); N_points_o = (N_points_s = 6);
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

      // declaration of the scalars and vectors needed in the integrations
      std::complex<double> ITo_ITs_G, ITo_r_dot_ITs_G_rprime, ITo_n_hat_X_r_dot_ITs_G_rprime, IDTo_l_hat_dot_r_ITs_G;
      std::complex<double> ITo_n_hat_X_r_dot_r_X_ITs_grad_G;
      std::complex<double> ITo_r_ITs_G[3], ITo_ITs_G_rprime[3], IDTo_l_hat_ITs_G[3];
      std::complex<double> ITo_ITs_grad_G[3], ITo_r_X_ITs_grad_G[3], ITo_n_hat_X_r_X_ITs_grad_G[3];

      ITo_ITs_free (ITo_ITs_G, ITo_r_ITs_G, ITo_ITs_G_rprime, ITo_r_dot_ITs_G_rprime, ITo_n_hat_X_r_dot_ITs_G_rprime, ITo_ITs_grad_G, ITo_r_X_ITs_grad_G, ITo_n_hat_X_r_dot_r_X_ITs_grad_G, ITo_n_hat_X_r_X_ITs_grad_G, triangles_test[r], triangles_src[s], k, N_points_o, N_points_s, EXTRACT_1_R, EXTRACT_R);

      const std::complex<double> ITo_n_hat_X_r_ITs_G[3] = {n_hat[1]*ITo_r_ITs_G[2] - n_hat[2]*ITo_r_ITs_G[1],
                                                           n_hat[2]*ITo_r_ITs_G[0] - n_hat[0]*ITo_r_ITs_G[2],
                                                           n_hat[0]*ITo_r_ITs_G[1] - n_hat[1]*ITo_r_ITs_G[0]};
      const std::complex<double> n_hat_dot_ITo_r_X_ITs_grad_G(n_hat[0]*ITo_r_X_ITs_grad_G[0] + n_hat[1]*ITo_r_X_ITs_grad_G[1] + n_hat[2]*ITo_r_X_ITs_grad_G[2]);

      if ((IS_TOUCH) && (nE_tmp || nH_tmp)) IDTo_ITs_free(IDTo_l_hat_dot_r_ITs_G, IDTo_l_hat_ITs_G, triangles_test[r], triangles_src[s], k, 3, N_points_s, EXTRACT_1_R, EXTRACT_R);

      for (unsigned int p=0 ; p<RWGsIndexes_test.size() ; ++p) {
        const int index_p = RWGsIndexes_test[p];
        const int local_number_edge_p = test_RWGs[index_p].number;
        const double l_p = test_RWGs[index_p].length;
        const double sign_edge_p = triangleTest_signsInRWGs[p];
        const double C_p = sign_edge_p * l_p * 0.5/triangles_test[r].A;
        double *r_p, n_hat_X_r_p[3];
        if (triangleTest_indexesInRWGs[p]==0) r_p = test_RWGs[index_p].vertexesCoord_0;
        else r_p = test_RWGs[index_p].vertexesCoord_3;
        cross3D(n_hat_X_r_p, n_hat, r_p);
        const int IS_CFIE = testRWGNumber_CFIE_OK(local_number_edge_p);
        const bool tEJ = tE_tmp, nEJ = nE_tmp, tHJ = (tH_tmp * IS_CFIE ), nHJ = (nH_tmp * IS_CFIE);

        // temporary elements for Z_tE_J
        const std::complex<double> p_dot_ITo_ITs_G_rprime(r_p[0]*ITo_ITs_G_rprime[0] + r_p[1]*ITo_ITs_G_rprime[1] + r_p[2]*ITo_ITs_G_rprime[2]);

        // temporary elements for Z_nE_J
        const std::complex<double> n_hat_X_r_p_dot_ITo_ITs_G_rprime(n_hat_X_r_p[0]*ITo_ITs_G_rprime[0] + n_hat_X_r_p[1]*ITo_ITs_G_rprime[1] + n_hat_X_r_p[2]*ITo_ITs_G_rprime[2]);
        const std::complex<double> n_hat_X_r_p_dot_ITo_ITs_grad_G(n_hat_X_r_p[0]*ITo_ITs_grad_G[0] + n_hat_X_r_p[1]*ITo_ITs_grad_G[1] + n_hat_X_r_p[2]*ITo_ITs_grad_G[2]);

        // temporary elements for Z_tH_J
        const std::complex<double> r_p_X_ITo_ITs_grad_G[3] = {r_p[1]*ITo_ITs_grad_G[2] - r_p[2]*ITo_ITs_grad_G[1],
                                                              r_p[2]*ITo_ITs_grad_G[0] - r_p[0]*ITo_ITs_grad_G[2],
                                                              r_p[0]*ITo_ITs_grad_G[1] - r_p[1]*ITo_ITs_grad_G[0]};
        const std::complex<double> ITo_r_rp_X_ITs_grad_G[3] = {ITo_r_X_ITs_grad_G[0]-r_p_X_ITo_ITs_grad_G[0],
                                                               ITo_r_X_ITs_grad_G[1]-r_p_X_ITo_ITs_grad_G[1],
                                                               ITo_r_X_ITs_grad_G[2]-r_p_X_ITo_ITs_grad_G[2]};
        const double n_hat_X_rp_dot_IT_r(n_hat_X_r_p[0]*IT_r[0] + n_hat_X_r_p[1]*IT_r[1] + n_hat_X_r_p[2]*IT_r[2]);

        // temporary elements for Z_nH_J
        const std::complex<double> n_hat_X_r_p_X_ITo_ITs_grad_G[3] = {n_hat_X_r_p[1]*ITo_ITs_grad_G[2] - n_hat_X_r_p[2]*ITo_ITs_grad_G[1],
                                                                      n_hat_X_r_p[2]*ITo_ITs_grad_G[0] - n_hat_X_r_p[0]*ITo_ITs_grad_G[2],
                                                                      n_hat_X_r_p[0]*ITo_ITs_grad_G[1] - n_hat_X_r_p[1]*ITo_ITs_grad_G[0]};
        const std::complex<double> ITo_n_hat_X_r_rp_X_ITs_grad_G[3] = {ITo_n_hat_X_r_X_ITs_grad_G[0] - n_hat_X_r_p_X_ITo_ITs_grad_G[0], 
                                                                       ITo_n_hat_X_r_X_ITs_grad_G[1] - n_hat_X_r_p_X_ITo_ITs_grad_G[1],
                                                                       ITo_n_hat_X_r_X_ITs_grad_G[2] - n_hat_X_r_p_X_ITo_ITs_grad_G[2]};
        const std::complex<double> n_hat_X_r_p_dot_ITo_r_X_ITs_grad_G(n_hat_X_r_p[0]*ITo_r_X_ITs_grad_G[0] + n_hat_X_r_p[1]*ITo_r_X_ITs_grad_G[1] + n_hat_X_r_p[2]*ITo_r_X_ITs_grad_G[2]);

        for (unsigned int q=0 ; q<RWGsIndexes_src.size() ; ++q) {
          const int index_q = RWGsIndexes_src[q];
          const int local_number_edge_q = src_RWGs[index_q].number;
          const double l_q = src_RWGs[index_q].length;
          const double sign_edge_q = triangleSrc_signsInRWGs[q];
          const double C_pq = C_p * sign_edge_q * l_q * 0.5/triangles_src[s].A;
          double *r_q;
          if (triangleSrc_indexesInRWGs[q]==0) r_q = src_RWGs[index_q].vertexesCoord_0;
          else r_q = src_RWGs[index_q].vertexesCoord_3;
          const double rp_dot_rq(r_p[0]*r_q[0] + r_p[1]*r_q[1] + r_p[2]*r_q[2]);
          const double rp_plus_rq_dot_ITr((r_p[0]+r_q[0])*IT_r[0] + (r_p[1]+r_q[1])*IT_r[1] + (r_p[2]+r_q[2])*IT_r[2]);
          std::complex<double> z_pq, D_mn_1, D_mn_2;

          // <f_p ; EFIE> : Z_tE_J computation. Z_tH_M = eps/mu * Z_tE_J
          if (tEJ) {
            D_mn_1 = (-4.0 * C_pq) * ITo_ITs_G;
            D_mn_2 = C_pq * (ITo_r_dot_ITs_G_rprime - (ITo_r_ITs_G[0]*r_q[0] + ITo_r_ITs_G[1]*r_q[1] + ITo_r_ITs_G[2]*r_q[2]) - p_dot_ITo_ITs_G_rprime + rp_dot_rq * ITo_ITs_G);
            z_pq = z_pq_factor * (D_mn_1 + k_square * D_mn_2); // z_pq_factor = 1.0/(I*w*eps)
            if (tds_approx && IS_SAME_TR) z_pq += 0.5 * C_pq * Z_s * ( IT_r_square - rp_plus_rq_dot_ITr + rp_dot_rq * triangles_test[r].A );
            Z_tE_J (local_number_edge_p, local_number_edge_q) += (signSurfObs * signSurfSrc) * z_pq;
          }
          // <n x f_p ; EFIE> : Z_nE_J computation. Z_nH_M = eps/mu * Z_nE_J
          if (nEJ) {
            if (IS_TOUCH) D_mn_1 = 2.0 * C_pq * (-IDTo_l_hat_dot_r_ITs_G + (r_p[0]*IDTo_l_hat_ITs_G[0] + r_p[1]*IDTo_l_hat_ITs_G[1] + r_p[2]*IDTo_l_hat_ITs_G[2]));
            else D_mn_1 = 2.0 * C_pq * ( n_hat_dot_ITo_r_X_ITs_grad_G - n_hat_X_r_p_dot_ITo_ITs_grad_G);
            D_mn_2 = C_pq * (ITo_n_hat_X_r_dot_ITs_G_rprime - (ITo_n_hat_X_r_ITs_G[0]*r_q[0] + ITo_n_hat_X_r_ITs_G[1]*r_q[1] + ITo_n_hat_X_r_ITs_G[2]*r_q[2]) - n_hat_X_r_p_dot_ITo_ITs_G_rprime + (n_hat_X_r_p[0]*r_q[0] + n_hat_X_r_p[1]*r_q[1] + n_hat_X_r_p[2]*r_q[2]) * ITo_ITs_G);
            z_pq = z_pq_factor * ( D_mn_1 + k_square*D_mn_2); // z_pq_factor = 1.0/(I*w*eps)
            Z_nE_J (local_number_edge_p, local_number_edge_q) += (signSurfObs * signSurfSrc) * z_pq;
          }
          // <f_p ; MFIE> : Z_tH_J computation. Z_tE_M = -Z_tH_J
          if (tHJ) {
            if (!IS_SAME_TR) z_pq = (signSurfObs * signSurfSrc * C_pq) * ((r_p[0]-r_q[0])*ITo_r_rp_X_ITs_grad_G[0] + (r_p[1]-r_q[1])*ITo_r_rp_X_ITs_grad_G[1] + (r_p[2]-r_q[2])*ITo_r_rp_X_ITs_grad_G[2]);
            // else we have: z_pq = 0.5 * <f_p ; n x f_q>
            else z_pq = (signSurfObs * 0.5 * C_pq) * ( n_hat_X_rp_dot_IT_r + (r_q[0]*IT_n_hat_X_r[0] + r_q[1]*IT_n_hat_X_r[1] + r_q[2]*IT_n_hat_X_r[2]) - (r_q[0]*n_hat_X_r_p[0] + r_q[1]*n_hat_X_r_p[1] + r_q[2]*n_hat_X_r_p[2]) * triangles_test[r].A );
            Z_tH_J (local_number_edge_p, local_number_edge_q) += z_pq;
          }
          // <n x f_p ; MFIE> : Z_nH_J computation. Z_nE_M = -Z_nH_J
          if (nHJ) {
            if (!IS_SAME_TR) z_pq = (-signSurfObs * signSurfSrc * C_pq) * ( ITo_n_hat_X_r_dot_r_X_ITs_grad_G + (r_q[0]*ITo_n_hat_X_r_rp_X_ITs_grad_G[0] + r_q[1]*ITo_n_hat_X_r_rp_X_ITs_grad_G[1] + r_q[2]*ITo_n_hat_X_r_rp_X_ITs_grad_G[2]) - n_hat_X_r_p_dot_ITo_r_X_ITs_grad_G  );
            // else we have: z_pq = 0.5 * <n x f_p ; n x f_q> = 0.5 * <f_p ; f_q>
            else z_pq = (signSurfObs * 0.5 * C_pq) * ( IT_r_square - rp_plus_rq_dot_ITr + rp_dot_rq * triangles_test[r].A );
            Z_nH_J (local_number_edge_p, local_number_edge_q) += z_pq;
          }

        } // for q
      } // for p
      
    } // for s
  } // for r
}


