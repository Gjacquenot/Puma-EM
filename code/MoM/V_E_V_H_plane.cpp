#include <iostream>
#include <complex>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <vector>
#include <algorithm>

using namespace blitz;

#include "EMConstants.h"
#include "GK_triangle.h"
#include "triangle_int.h"
#include "mesh.h"
#include "V_E_V_H.h"

void V_EJ_HJ_plane (blitz::Array<std::complex<double>, 1> V_tE_J,
                    blitz::Array<std::complex<double>, 1> V_nE_J,
                    blitz::Array<std::complex<double>, 1> V_tH_J,
                    blitz::Array<std::complex<double>, 1> V_nH_J,
                    const blitz::Array<std::complex<double>, 1>& E_0,
                    const blitz::Array<double, 1>& k_hat,
                    const blitz::Array<double, 1>& r_ref,
                    const blitz::Array<int, 1>& numbers_RWG_test,
                    const blitz::Array<int, 1>& RWGNumber_CFIE_OK,
                    const blitz::Array<int, 2>& RWGNumber_signedTriangles,
                    const blitz::Array<double, 2>& RWGNumber_vertexesCoord,
                    const blitz::Array<double, 2>& RWGNumber_oppVertexesCoord,
                    const double w,
                    const complex<double>& eps_r,
                    const complex<double>& mu_r,
                    const int FULL_PRECISION)
{
  // def of k, mu_i, eps_i
  blitz::Range all = blitz::Range::all();
  int N_RWG_test = numbers_RWG_test.size();
  std::complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu);

  // weights and abscissas for triangle integration
  double sum_weigths;
  const double *xi, *eta, *weigths;
  // triangle integration precision. Possible values for N_points: 1, 3, 6, 9, 12, 13
  int N_points;
  if (FULL_PRECISION!=0) N_points = 6; 
  else N_points = 3;
  IT_points (xi, eta, weigths, sum_weigths, N_points);

  blitz::TinyVector<double, 3> kHat, rRef;
  for (int i=0 ; i<3 ; ++i) kHat(i) = k_hat(i);
  for (int i=0 ; i<3 ; ++i) rRef(i) = r_ref(i);

  V_tE_J = 0.0; V_tH_J = 0.0;
  V_nE_J = 0.0; V_nH_J = 0.0;
  std::vector<RWG> test_RWGs;
  test_RWGs.reserve(N_RWG_test);
  for (int i=0 ; i<N_RWG_test ; ++i) {
    const int RWGnumber = numbers_RWG_test(i);
    blitz::TinyVector<int, 2> triangle_numbers, triangle_signs;
    triangle_numbers(0) = abs(RWGNumber_signedTriangles(RWGnumber, 0));
    triangle_numbers(1) = abs(RWGNumber_signedTriangles(RWGnumber, 1));
    triangle_signs(0) = 1;
    triangle_signs(1) = -1;
    // the nodes of the RWG
    blitz::Array<double, 1> r0(3), r1(3), r2(3), r3(3);
    r0 = RWGNumber_oppVertexesCoord(RWGnumber, blitz::Range(0,2));
    r3 = RWGNumber_oppVertexesCoord(RWGnumber, blitz::Range(3,5));
    r1 = RWGNumber_vertexesCoord(RWGnumber, blitz::Range(0,2));
    r2 = RWGNumber_vertexesCoord(RWGnumber, blitz::Range(3,5));
    test_RWGs.push_back(RWG(RWGnumber, triangle_numbers, triangle_signs, r0, r1, r2, r3));
  }
  // triangles
  std::vector< Dictionary<int, int> > testTriangleToRWG;
  testTriangleToRWG.reserve(N_RWG_test*2);
  for (int i=0 ; i<test_RWGs.size() ; ++i) {
    testTriangleToRWG.push_back(Dictionary<int, int>(test_RWGs[i].triangleNumbers(0), test_RWGs[i].number));
    testTriangleToRWG.push_back(Dictionary<int, int>(test_RWGs[i].triangleNumbers(1), test_RWGs[i].number));
  }
  sort(testTriangleToRWG.begin(), testTriangleToRWG.end());
  std::vector<Triangle> triangles_test;
  constructVectorTriangles(triangles_test, test_RWGs, testTriangleToRWG);

  // geometrical entities
  blitz::TinyVector<double, 3> r0, r1, r2, r_obs;
  std::complex<double> ITo_r_dot_H_inc, ITo_r_dot_E_inc, ITo_n_hat_X_r_dot_H_inc, ITo_n_hat_X_r_dot_E_inc;
  blitz::TinyVector<std::complex<double>, 3> ITo_H_inc, ITo_E_inc;
  // computation of H_0
  blitz::TinyVector<std::complex<double>, 3> E0, H0;
  for (int i=0 ; i<3 ; ++i) E0(i) = E_0(i);
  H0 = k_hat(1)*E_0(2)-k_hat(2)*E_0(1),
       k_hat(2)*E_0(0)-k_hat(0)*E_0(2),
       k_hat(0)*E_0(1)-k_hat(1)*E_0(0);
  H0 *= sqrt(eps/mu);

  for (int r=0 ; r<triangles_test.size() ; ++r) { // loop on the observation triangles
      // the RWGs concerned by the test triangle
      std::vector<int> RWGsIndexes_test(triangles_test[r].RWGIndexes);
      std::vector<int> triangleTest_indexesInRWGs(triangles_test[r].indexesInRWGs);
      std::vector<double> triangleTest_signsInRWGs(triangles_test[r].signInRWG);
      ITo_E_inc = 0.0; // TinyVector<complex<double>, 3>
      ITo_H_inc = 0.0; // TinyVector<complex<double>, 3>
      ITo_r_dot_E_inc = 0.0; // complex<double>
      ITo_r_dot_H_inc = 0.0; // complex<double>
      ITo_n_hat_X_r_dot_E_inc = 0.0; // complex<double>
      ITo_n_hat_X_r_dot_H_inc = 0.0; // complex<double>

      r0 = triangles_test[r].r_nodes(0);
      r1 = triangles_test[r].r_nodes(1);
      r2 = triangles_test[r].r_nodes(2);

      // triangle integration
      for (int i=0 ; i<N_points ; ++i) {
        r_obs = r0 * xi[i] + r1 * eta[i] + r2 * (1-xi[i]-eta[i]);
        blitz::TinyVector<double, 3> n_hat_X_r, n_hat = triangles_test[r].n_hat;
        n_hat_X_r = n_hat(1)*r_obs(2)-n_hat(2)*r_obs(1),
                    n_hat(2)*r_obs(0)-n_hat(0)*r_obs(2),
                    n_hat(0)*r_obs(1)-n_hat(1)*r_obs(0);

        blitz::TinyVector<std::complex<double>, 3> H_inc_i, E_inc_i;
        // computation of ITo_E_inc due to a dipole located at r_dip
        E_inc_i = E0 * exp(-I*k * dot(kHat, r_obs - rRef) ) * weigths[i];
        ITo_E_inc += E_inc_i;
        ITo_r_dot_E_inc += dot(r_obs, E_inc_i);
        ITo_n_hat_X_r_dot_E_inc += dot(n_hat_X_r, E_inc_i);

        // computation of ITo_H_inc due to a dipole located at r_dip
        H_inc_i = H0 * exp(-I*k * dot(kHat, r_obs - rRef) ) * weigths[i];
        ITo_H_inc += H_inc_i;
        ITo_r_dot_H_inc += dot(r_obs, H_inc_i);
        ITo_n_hat_X_r_dot_H_inc += dot(n_hat_X_r, H_inc_i);
 
      }
      const double norm_factor = triangles_test[r].A/sum_weigths;

      ITo_E_inc *= norm_factor;
      ITo_H_inc *= norm_factor;
      ITo_r_dot_E_inc *= norm_factor;
      ITo_r_dot_H_inc *= norm_factor;
      ITo_n_hat_X_r_dot_E_inc *= norm_factor;
      ITo_n_hat_X_r_dot_H_inc *= norm_factor;

      for (int p=0 ; p<RWGsIndexes_test.size() ; ++p) {
        const int index_p = RWGsIndexes_test[p];
        const int local_number_edge_p = test_RWGs[index_p].number;
        const double l_p = test_RWGs[index_p].length;
        const double sign_edge_p = triangleTest_signsInRWGs[p];
        const double C_rp = sign_edge_p * l_p * 0.5/triangles_test[r].A;
        blitz::TinyVector<double, 3> r_p;
        if (triangleTest_indexesInRWGs[p]==0) r_p = test_RWGs[index_p].vertexesCoord(0);
        else if (triangleTest_indexesInRWGs[p]==1) r_p = test_RWGs[index_p].vertexesCoord(3);
        const TinyVector<double, 3> n_hat_X_r_p(cross(triangles_test[r].n_hat, r_p));

        V_tE_J (local_number_edge_p) += -C_rp * (ITo_r_dot_E_inc - sum(r_p * ITo_E_inc)); // -<f_m ; E_inc> 
        V_tH_J (local_number_edge_p) += -C_rp * (ITo_r_dot_H_inc - sum(r_p * ITo_H_inc)); // -<f_m ; H_inc>
        V_nE_J (local_number_edge_p) += -C_rp * (ITo_n_hat_X_r_dot_E_inc - sum(n_hat_X_r_p * ITo_E_inc)); // -<n_hat x f_m ; E_inc> 
        V_nH_J (local_number_edge_p) += -C_rp * (ITo_n_hat_X_r_dot_H_inc - sum(n_hat_X_r_p * ITo_H_inc)); // -<n_hat x f_m ; H_inc> 
      }
  }
}


void V_CFIE_plane (blitz::Array<std::complex<float>, 1> V_CFIE,
                   const blitz::Array<std::complex<float>, 1>& CFIE,
                   const blitz::Array<std::complex<double>, 1>& E_0,
                   const blitz::Array<double, 1>& k_hat,
                   const blitz::Array<double, 1>& r_ref,
                   const blitz::Array<int, 1>& numbers_RWG_test,
                   const blitz::Array<int, 1>& RWGNumber_CFIE_OK,
                   const blitz::Array<float, 2>& RWGNumber_trianglesCoord,
                   const double w,
                   const complex<double>& eps_r,
                   const complex<double>& mu_r,
                   const int FULL_PRECISION)
{
  // def of k, mu_i, eps_i
  Range all = Range::all();
  int N_RWG_test = numbers_RWG_test.size();
  complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu);
  const complex<double> tE = CFIE(0), nE = CFIE(1), tH = CFIE(2), nH = CFIE(3);

  double sum_weigths;
  const double *xi, *eta, *weigths;
  int N_points;
  if (FULL_PRECISION!=0) N_points = 6; 
  else N_points = 3;
  IT_points (xi, eta, weigths, sum_weigths, N_points);

  // r_ref is the reference for the phase of the incoming plane wave
  blitz::TinyVector<double, 3> kHat, rRef;
  for (int i=0 ; i<3 ; ++i) kHat(i) = k_hat(i);
  for (int i=0 ; i<3 ; ++i) rRef(i) = r_ref(i);

  complex<double> ITo_r_dot_H_inc, ITo_r_dot_E_inc, ITo_n_hat_X_r_dot_H_inc, ITo_n_hat_X_r_dot_E_inc;

  // computation of H_0
  blitz::TinyVector<std::complex<double>, 3> E0, H0;
  for (int i=0 ; i<3 ; ++i) E0(i) = E_0(i);
  H0 = k_hat(1)*E_0(2)-k_hat(2)*E_0(1),
       k_hat(2)*E_0(0)-k_hat(0)*E_0(2),
       k_hat(0)*E_0(1)-k_hat(1)*E_0(0);
  H0 *= sqrt(eps/mu);

  V_CFIE = 0.0;
  blitz::TinyVector<double, 3> r0, r1, r2, r_obs;
  for (int rwg=0 ; rwg<N_RWG_test ; ++rwg) { // loop on the RWGs
    for (int tr = 0 ; tr<2 ; ++tr) {
      double l_p;
      blitz::Array<double, 1> rt0(3), rt1(3), rt2(3), r_opp(3);
      if (tr==0) {
        for (int i=0; i<3; i++) { 
          rt0(i) = RWGNumber_trianglesCoord(rwg, i);
          rt1(i) = RWGNumber_trianglesCoord(rwg, i+3);
          rt2(i) = RWGNumber_trianglesCoord(rwg, i+6);
        }
        for (int i=0 ; i<3 ; i++) {
          r0(i) = rt0(i);
          r1(i) = rt1(i);
          r2(i) = rt2(i);
        }
        r_opp = rt0;
        l_p = sqrt(sum((rt1-rt2) * (rt1-rt2)));
      }
      else{
        for (int i=0; i<3; i++) { 
          rt0(i) = RWGNumber_trianglesCoord(rwg, i+6);
          rt1(i) = RWGNumber_trianglesCoord(rwg, i+3);
          rt2(i) = RWGNumber_trianglesCoord(rwg, i+9);
        }
        for (int i=0 ; i<3 ; i++) {
          r0(i) = rt0(i);
          r1(i) = rt1(i);
          r2(i) = rt2(i);
        }
        r_opp = rt2;
        l_p = sqrt(sum((rt0-rt1) * (rt0-rt1)));
      }
      Triangle triangle(r0, r1, r2, 0);
      blitz::TinyVector<std::complex<double>, 3> ITo_E_inc, ITo_H_inc;
      ITo_E_inc = 0.0; // Array<complex<double>, 1>
      ITo_H_inc = 0.0; // Array<complex<double>, 1>
      ITo_r_dot_E_inc = 0.0; // complex<double>
      ITo_r_dot_H_inc = 0.0; // complex<double>
      ITo_n_hat_X_r_dot_E_inc = 0.0; // complex<double>
      ITo_n_hat_X_r_dot_H_inc = 0.0; // complex<double>

      r0 = triangle.r_nodes(0);
      r1 = triangle.r_nodes(1);
      r2 = triangle.r_nodes(2);

      for (int i=0 ; i<N_points ; ++i) {
        blitz::TinyVector<double, 3>  r_obs = r0 * xi[i] + r1 * eta[i] + r2 * (1-xi[i]-eta[i]);
        blitz::TinyVector<double, 3> n_hat_X_r;
        n_hat_X_r = triangle.n_hat(1)*r_obs(2)-triangle.n_hat(2)*r_obs(1),
                    triangle.n_hat(2)*r_obs(0)-triangle.n_hat(0)*r_obs(2),
                    triangle.n_hat(0)*r_obs(1)-triangle.n_hat(1)*r_obs(0);

        blitz::TinyVector<std::complex<double>, 3> H_inc_i, E_inc_i;
        // computation of ITo_E_inc due to a dipole located at r_dip
        E_inc_i = E0 * exp(-I*k * dot(kHat, r_obs - rRef) ) * weigths[i];
        ITo_E_inc += E_inc_i;
        ITo_r_dot_E_inc += dot(r_obs, E_inc_i);
        ITo_n_hat_X_r_dot_E_inc += dot(n_hat_X_r, E_inc_i);

        // computation of ITo_H_inc due to a dipole located at r_dip
        H_inc_i = H0 * exp(-I*k * dot(kHat, r_obs - rRef) ) * weigths[i];
        ITo_H_inc += H_inc_i;
        ITo_r_dot_H_inc += dot(r_obs, H_inc_i);
        ITo_n_hat_X_r_dot_H_inc += dot(n_hat_X_r, H_inc_i);
      }
      double norm_factor = triangle.A/sum_weigths;

      ITo_E_inc *= norm_factor;
      ITo_H_inc *= norm_factor;
      ITo_r_dot_E_inc *= norm_factor;
      ITo_r_dot_H_inc *= norm_factor;
      ITo_n_hat_X_r_dot_E_inc *= norm_factor;
      ITo_n_hat_X_r_dot_H_inc *= norm_factor;

      const int local_number_edge_p = numbers_RWG_test(rwg);
      const int sign_edge_p = (tr==0) ? 1 : -1;
      const double C_rp = sign_edge_p * l_p * 0.5/triangle.A;
      blitz::TinyVector<double, 3> r_p;
      for (int i=0; i<3 ; ++i) r_p(i) = r_opp(i);
      const TinyVector<double, 3> n_hat_X_r_p(cross(triangle.n_hat, r_p));

      complex<double> tmpResult(0.0, 0.0);
      tmpResult -= tE * C_rp * (ITo_r_dot_E_inc - sum(r_p * ITo_E_inc)); // -<f_m ; E_inc>
      if (RWGNumber_CFIE_OK(local_number_edge_p) == 1) {
        tmpResult -= nE * C_rp * (ITo_n_hat_X_r_dot_E_inc - sum(n_hat_X_r_p * ITo_E_inc)); // -<n_hat x f_m ; E_inc> 
        tmpResult -= tH * C_rp * (ITo_r_dot_H_inc - sum(r_p * ITo_H_inc)); // -<f_m ; H_inc>
        tmpResult -= nH * C_rp * (ITo_n_hat_X_r_dot_H_inc - sum(n_hat_X_r_p * ITo_H_inc)); // -<n_hat x f_m ; H_inc> 
      }
      V_CFIE(local_number_edge_p) += tmpResult;
    }
  }
}

void local_V_CFIE_plane (blitz::Array<std::complex<float>, 1>& V_CFIE,
                          const blitz::Array<std::complex<double>, 1>& E_0,
                          const blitz::Array<double, 1>& k_hat,
                          const blitz::Array<double, 1>& r_ref,
                          const LocalMesh & local_target_mesh,
                          const double w,
                          const std::complex<double>& eps_r,
                          const std::complex<double>& mu_r,
                          const blitz::Array<std::complex<float>, 1>& CFIE,
                          const int FULL_PRECISION)
{
  // We now compute the excitation vectors
  V_CFIE.resize(local_target_mesh.N_local_RWG);
  V_CFIE_plane (V_CFIE, CFIE, E_0, k_hat, r_ref, local_target_mesh.reallyLocalRWGNumbers, local_target_mesh.localRWGNumber_CFIE_OK, local_target_mesh.localRWGNumber_trianglesCoord, w, eps_r, mu_r, FULL_PRECISION);
}


// void V_E_V_H_horn_BBHA (Array<complex<double>,2>& V_EJ, Array<complex<double>,2>& V_HJ, Array<complex<double>,2>& V_EM, Array<complex<double>,2>& V_HM, const mesh & MESH_TARGET, const mesh & MESH_ANTENNA, const TinyVector<double,3>& r_ant, const TinyVector<double,3>& x_hat_ant, const TinyVector<double,3>& y_hat_ant, const G_EJ_grid & G_EJ_tab_ant_object, const G_HJ_grid & G_HJ_tab_ant_object, const layers_constants & LC, const G_EJ_grid & G_EJ_tab_ant_object_dual, const G_HJ_grid & G_HJ_tab_ant_object_dual, const layers_constants & LC_dual) {

//   int i, j, N_TR_ANT = MESH_ANTENNA.triangles.rows(), E = MESH_TARGET.edges.rows()/2; 
//   double sum_weigths;
//   const double *xi, *eta, *weigths;
//   int N_points = 9;
//   IT_points (xi, eta, weigths, sum_weigths, N_points);

//   double rho_1 = 0.29015095530181, rho_2 = 0.25187444860284, a_1 = 0.238, b_1 = 0.138; // BBHA antenna geometrical parameters
//   double E_0 = sqrt(4.0*sqrt(LC.mu_0/LC.eps_0) / (a_1*b_1)); // so that P_rad = 1 (Balanis eq. 12-51)
//   V_EJ = 0.0;
//   V_HJ = 0.0;
//   V_EM = 0.0;
//   V_HM = 0.0;
//   Array<complex<double>, 2> V_EJ_triangle (E, 3), V_EM_triangle (E, 3), V_HJ_triangle (E, 3), V_HM_triangle (E, 3);
//   Array<complex<double>, 2> V_EJ_j (E, 3), V_EM_j (E, 3), V_HJ_j (E, 3), V_HM_j (E, 3);
//   for (i=0 ; i<N_TR_ANT ; i++) { // loop on antenna triangles
//     cout << "\r" << i*100/N_TR_ANT;
//     flush(cout);
//     V_EJ_triangle = 0.0;
//     V_EM_triangle = 0.0;
//     V_HJ_triangle = 0.0;
//     V_HM_triangle = 0.0;
//     complex<double> term_1;
//     TinyVector<double, 3> r_j, r_J_on_antenna;
//     TinyVector<complex<double>, 3> J_ant, M_ant;
//     triangle T_antenna_i;
//     construct_triangle (T_antenna_i, MESH_ANTENNA, i);
//     for (j=0 ; j<N_points ; j++) { // loop on points of antenna triangle
//       r_J_on_antenna = T_antenna_i.r_nodes (0) * xi[j] + T_antenna_i.r_nodes (1) * eta[j] + T_antenna_i.r_nodes (2) * (1-xi[j]-eta[j]);
//       r_j = r_ant + r_J_on_antenna (0) * x_hat_ant + r_J_on_antenna (1) * y_hat_ant;
//       term_1 = E_0 * cos(M_PI/a_1*r_J_on_antenna (0)) * exp(-I*LC.k_0*(pow2(r_J_on_antenna (0))/rho_2 + pow2(r_J_on_antenna (1))/rho_1)/2.0);
//       J_ant = -1.0/sqrt(LC.mu_0/LC.eps_0) * term_1 * y_hat_ant;
//       V_dipole_body (V_EJ_j, V_HJ_j, MESH_TARGET, G_EJ_tab_ant_object, G_HJ_tab_ant_object, J_ant, r_j, LC);
//       V_EJ_triangle += V_EJ_j * weigths[j];
//       V_HJ_triangle += V_HJ_j * weigths[j];

//       M_ant = term_1 * x_hat_ant;
//       V_dipole_body (V_HM_j, V_EM_j, MESH_TARGET, G_EJ_tab_ant_object_dual, G_HJ_tab_ant_object_dual, M_ant, r_j, LC_dual);
//       V_EM_triangle -= V_EM_j * weigths[j];
//       V_HM_triangle += V_HM_j * weigths[j];
//     }
//     double norm_factor = T_antenna_i.A/sum_weigths;
//     V_EJ += V_EJ_triangle * norm_factor;
//     V_HJ += V_HJ_triangle * norm_factor;
//     V_EM += V_EM_triangle * norm_factor;
//     V_HM += V_HM_triangle * norm_factor;
//   }

// }
