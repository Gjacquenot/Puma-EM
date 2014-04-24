#include <iostream>
#include <complex>
#include <blitz/array.h>
#include <vector>
#include <algorithm>

using namespace std;

#include "EMConstants.h"
#include "GK_triangle.h"
#include "triangle_int.h"
#include "dictionary.h"
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
                    const std::complex<double>& eps_r,
                    const std::complex<double>& mu_r,
                    const int FULL_PRECISION)
{
  // def of k, mu_i, eps_i
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

  double kHat[3], rRef[3];
  for (int i=0 ; i<3 ; ++i) kHat[i] = k_hat(i);
  for (int i=0 ; i<3 ; ++i) rRef[i] = r_ref(i);

  V_tE_J = 0.0; V_tH_J = 0.0;
  V_nE_J = 0.0; V_nH_J = 0.0;
  std::vector<RWG> test_RWGs;
  test_RWGs.reserve(N_RWG_test);
  for (int i=0 ; i<N_RWG_test ; ++i) {
    const int RWGnumber = numbers_RWG_test(i);
    int triangle_numbers[2], triangle_signs[2];
    triangle_numbers[0] = abs(RWGNumber_signedTriangles(RWGnumber, 0));
    triangle_numbers[1] = abs(RWGNumber_signedTriangles(RWGnumber, 1));
    triangle_signs[0] = 1;
    triangle_signs[1] = -1;
    // the nodes of the RWG
    const double r0[3] = { RWGNumber_oppVertexesCoord(RWGnumber, 0), RWGNumber_oppVertexesCoord(RWGnumber, 1), RWGNumber_oppVertexesCoord(RWGnumber, 2)};
    const double r3[3] = { RWGNumber_oppVertexesCoord(RWGnumber, 3), RWGNumber_oppVertexesCoord(RWGnumber, 4), RWGNumber_oppVertexesCoord(RWGnumber, 5)};
    const double r1[3] = { RWGNumber_vertexesCoord(RWGnumber, 0), RWGNumber_vertexesCoord(RWGnumber, 1), RWGNumber_vertexesCoord(RWGnumber, 2)};
    const double r2[3] = { RWGNumber_vertexesCoord(RWGnumber, 3), RWGNumber_vertexesCoord(RWGnumber, 4), RWGNumber_vertexesCoord(RWGnumber, 5)};
    test_RWGs.push_back(RWG(RWGnumber, triangle_numbers, triangle_signs, r0, r1, r2, r3));
  }
  // triangles
  std::vector< Dictionary<int, int> > testTriangleToRWG;
  testTriangleToRWG.reserve(N_RWG_test*2);
  for (unsigned int i=0 ; i<test_RWGs.size() ; ++i) {
    testTriangleToRWG.push_back(Dictionary<int, int>(test_RWGs[i].triangleNumbers[0], test_RWGs[i].number));
    testTriangleToRWG.push_back(Dictionary<int, int>(test_RWGs[i].triangleNumbers[1], test_RWGs[i].number));
  }
  sort(testTriangleToRWG.begin(), testTriangleToRWG.end());
  std::vector<Triangle> triangles_test;
  constructVectorTriangles(triangles_test, test_RWGs, testTriangleToRWG);

  // geometrical entities
  const double *r0, *r1, *r2;
  std::complex<double> ITo_r_dot_H_inc, ITo_r_dot_E_inc, ITo_n_hat_X_r_dot_H_inc, ITo_n_hat_X_r_dot_E_inc;
  std::complex<double> ITo_H_inc[3], ITo_E_inc[3], E_inc_i[3], H_inc_i[3];
  // computation of H_0
  std::complex<double> E0[3], H0[3];
  for (int i=0 ; i<3 ; ++i) E0[i] = E_0(i);
  H0[0] = (k_hat(1)*E_0(2)-k_hat(2)*E_0(1)) * sqrt(eps/mu);
  H0[1] = (k_hat(2)*E_0(0)-k_hat(0)*E_0(2)) * sqrt(eps/mu);
  H0[2] = (k_hat(0)*E_0(1)-k_hat(1)*E_0(0)) * sqrt(eps/mu);

  for (unsigned int r=0 ; r<triangles_test.size() ; ++r) { // loop on the observation triangles
      // the RWGs concerned by the test triangle
      std::vector<int> RWGsIndexes_test(triangles_test[r].RWGIndexes);
      std::vector<int> triangleTest_indexesInRWGs(triangles_test[r].indexesInRWGs);
      std::vector<double> triangleTest_signsInRWGs(triangles_test[r].signInRWG);
      for (int i=0 ; i<3 ; ++i) ITo_E_inc[i] = 0.0; // Vector<complex<double>, 3>
      for (int i=0 ; i<3 ; ++i) ITo_H_inc[i] = 0.0; // Vector<complex<double>, 3>
      ITo_r_dot_E_inc = 0.0; // complex<double>
      ITo_r_dot_H_inc = 0.0; // complex<double>
      ITo_n_hat_X_r_dot_E_inc = 0.0; // complex<double>
      ITo_n_hat_X_r_dot_H_inc = 0.0; // complex<double>

      r0 = triangles_test[r].r_nodes_0;
      r1 = triangles_test[r].r_nodes_1;
      r2 = triangles_test[r].r_nodes_2;

      // triangle integration
      for (int j=0 ; j<N_points ; ++j) {
        double r_obs[3];
        r_obs[0] = r0[0] * xi[j] + r1[0] * eta[j] + r2[0] * (1-xi[j]-eta[j]);
        r_obs[1] = r0[1] * xi[j] + r1[1] * eta[j] + r2[1] * (1-xi[j]-eta[j]);
        r_obs[2] = r0[2] * xi[j] + r1[2] * eta[j] + r2[2] * (1-xi[j]-eta[j]);
        double n_hat_X_r[3];
        cross3D(n_hat_X_r, triangles_test[r].n_hat, r_obs);

        // computation of ITo_E_inc due to a dipole located at r_dip
        const double r_obs_rRef[3] = {r_obs[0] - rRef[0], r_obs[1] - rRef[1], r_obs[2] - rRef[2]};
        const std::complex<double> temp(exp(-I*k * (kHat[0] * r_obs_rRef[0] + kHat[1] * r_obs_rRef[1] + kHat[2] * r_obs_rRef[2]) ));
        for (int m=0 ; m<3 ; m++) E_inc_i[m] = E0[m] * temp * weigths[j];
        ITo_E_inc[0] += E_inc_i[0];
        ITo_E_inc[1] += E_inc_i[1];
        ITo_E_inc[2] += E_inc_i[2];
        ITo_r_dot_E_inc += r_obs[0] * E_inc_i[0] + r_obs[1] * E_inc_i[1] + r_obs[2] * E_inc_i[2];
        ITo_n_hat_X_r_dot_E_inc += n_hat_X_r[0] * E_inc_i[0] + n_hat_X_r[1] * E_inc_i[1] + n_hat_X_r[2] * E_inc_i[2];

        // computation of ITo_H_inc due to a dipole located at r_dip
        for (int m=0 ; m<3 ; m++) H_inc_i[m] = H0[m] * temp * weigths[j];
        ITo_H_inc[0] += H_inc_i[0];
        ITo_H_inc[1] += H_inc_i[1];
        ITo_H_inc[2] += H_inc_i[2];
        ITo_r_dot_H_inc += r_obs[0] * H_inc_i[0] + r_obs[1] * H_inc_i[1] + r_obs[2] * H_inc_i[2];
        ITo_n_hat_X_r_dot_H_inc += n_hat_X_r[0] * H_inc_i[0] + n_hat_X_r[1] * H_inc_i[1] + n_hat_X_r[2] * H_inc_i[2];
      }
      const double norm_factor = triangles_test[r].A/sum_weigths;

      for (int i=0 ; i<3 ; ++i) ITo_E_inc[i] *= norm_factor;
      for (int i=0 ; i<3 ; ++i) ITo_H_inc[i] *= norm_factor;
      ITo_r_dot_E_inc *= norm_factor;
      ITo_r_dot_H_inc *= norm_factor;
      ITo_n_hat_X_r_dot_E_inc *= norm_factor;
      ITo_n_hat_X_r_dot_H_inc *= norm_factor;

      for (unsigned int p=0 ; p<RWGsIndexes_test.size() ; ++p) {
        const int index_p = RWGsIndexes_test[p];
        const int local_number_edge_p = test_RWGs[index_p].number;
        const double l_p = test_RWGs[index_p].length;
        const double sign_edge_p = triangleTest_signsInRWGs[p];
        const double C_rp = sign_edge_p * l_p * 0.5/triangles_test[r].A;
        double r_p[3], n_hat_X_r_p[3];
        if (triangleTest_indexesInRWGs[p]==0) {
          for (int i=0 ; i<3 ; ++i) r_p[i] = test_RWGs[index_p].vertexesCoord_0[i];
        }
        else if (triangleTest_indexesInRWGs[p]==1) {
          for (int i=0 ; i<3 ; ++i) r_p[i] = test_RWGs[index_p].vertexesCoord_3[i];
        }
        cross3D(n_hat_X_r_p, triangles_test[r].n_hat, r_p);

        V_tE_J (local_number_edge_p) += -C_rp * (ITo_r_dot_E_inc - (r_p[0]*ITo_E_inc[0] + r_p[1]*ITo_E_inc[1] + r_p[2]*ITo_E_inc[2])); // -<f_m ; E_inc> 
        V_tH_J (local_number_edge_p) += -C_rp * (ITo_r_dot_H_inc - (r_p[0]*ITo_H_inc[0] + r_p[1]*ITo_H_inc[1] + r_p[2]*ITo_H_inc[2])); // -<f_m ; H_inc>
        V_nE_J (local_number_edge_p) += -C_rp * (ITo_n_hat_X_r_dot_E_inc - (n_hat_X_r_p[0]*ITo_E_inc[0] + n_hat_X_r_p[1]*ITo_E_inc[1] + n_hat_X_r_p[2]*ITo_E_inc[2])); // -<n_hat x f_m ; E_inc> 
        V_nH_J (local_number_edge_p) += -C_rp * (ITo_n_hat_X_r_dot_H_inc - (n_hat_X_r_p[0]*ITo_H_inc[0] + n_hat_X_r_p[1]*ITo_H_inc[1] + n_hat_X_r_p[2]*ITo_H_inc[2])); // -<n_hat x f_m ; H_inc> 
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
                   const std::complex<double>& eps_r,
                   const std::complex<double>& mu_r,
                   const int FULL_PRECISION)
{
  // def of k, mu_i, eps_i
  int N_RWG_test = numbers_RWG_test.size();
  std::complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu);
  const std::complex<double> tE = CFIE(0), nE = CFIE(1), tH = CFIE(2), nH = CFIE(3);

  double sum_weigths;
  const double *xi, *eta, *weigths;
  int N_points;
  if (FULL_PRECISION!=0) N_points = 6; 
  else N_points = 3;
  IT_points (xi, eta, weigths, sum_weigths, N_points);

  // r_ref is the reference for the phase of the incoming plane wave
  double kHat[3], rRef[3];
  for (int i=0 ; i<3 ; ++i) kHat[i] = k_hat(i);
  for (int i=0 ; i<3 ; ++i) rRef[i] = r_ref(i);

  // computation of H_0
  double r0[3], r1[3], r2[3], *r_opp;
  std::complex<double> ITo_r_dot_H_inc, ITo_r_dot_E_inc, ITo_n_hat_X_r_dot_H_inc, ITo_n_hat_X_r_dot_E_inc;
  std::complex<double> H_inc_i[3], ITo_H_inc[3], E_inc_i[3], ITo_E_inc[3];
  // computation of H_0
  std::complex<double> E0[3], H0[3];
  for (int i=0 ; i<3 ; ++i) E0[i] = E_0(i);
  H0[0] = (k_hat(1)*E_0(2)-k_hat(2)*E_0(1)) * sqrt(eps/mu);
  H0[1] = (k_hat(2)*E_0(0)-k_hat(0)*E_0(2)) * sqrt(eps/mu);
  H0[2] = (k_hat(0)*E_0(1)-k_hat(1)*E_0(0)) * sqrt(eps/mu);

  V_CFIE = 0.0;
  for (int rwg=0 ; rwg<N_RWG_test ; ++rwg) { // loop on the RWGs
    for (int tr = 0 ; tr<2 ; ++tr) {
      double l_p;
      if (tr==0) {
        for (int i=0; i<3; i++) { 
          r0[i] = RWGNumber_trianglesCoord(rwg, i);
          r1[i] = RWGNumber_trianglesCoord(rwg, i+3);
          r2[i] = RWGNumber_trianglesCoord(rwg, i+6);
        }
        r_opp = r0;
        double r1_r2[3] = {r1[0]-r2[0], r1[1]-r2[1], r1[2]-r2[2]};
        l_p = sqrt(dot3D(r1_r2, r1_r2));
      }
      else{
        for (int i=0 ; i<3 ; i++) {
          r0[i] = RWGNumber_trianglesCoord(rwg, i+6);
          r1[i] = RWGNumber_trianglesCoord(rwg, i+3);
          r2[i] = RWGNumber_trianglesCoord(rwg, i+9);
        }
        r_opp = r2;
        double r0_r1[3] = {r0[0]-r1[0], r0[1]-r1[1], r0[2]-r1[2]};
        l_p = sqrt(dot3D(r0_r1, r0_r1));
      }
      Triangle triangle(r0, r1, r2, 0);
      for (int i=0 ; i<3 ; ++i) ITo_E_inc[i] = 0.0; // Vector<complex<double>, 3>
      for (int i=0 ; i<3 ; ++i) ITo_H_inc[i] = 0.0; // Vector<complex<double>, 3>
      ITo_r_dot_E_inc = 0.0; // complex<double>
      ITo_r_dot_H_inc = 0.0; // complex<double>
      ITo_n_hat_X_r_dot_E_inc = 0.0; // complex<double>
      ITo_n_hat_X_r_dot_H_inc = 0.0; // complex<double>

      for (int j=0 ; j<N_points ; ++j) {
        double r_obs[3];
        r_obs[0] = r0[0] * xi[j] + r1[0] * eta[j] + r2[0] * (1-xi[j]-eta[j]);
        r_obs[1] = r0[1] * xi[j] + r1[1] * eta[j] + r2[1] * (1-xi[j]-eta[j]);
        r_obs[2] = r0[2] * xi[j] + r1[2] * eta[j] + r2[2] * (1-xi[j]-eta[j]);
        double n_hat_X_r[3];
        cross3D(n_hat_X_r, triangle.n_hat, r_obs);

        // computation of ITo_E_inc due to a dipole located at r_dip
        const double r_obs_rRef[3] = {r_obs[0] - rRef[0], r_obs[1] - rRef[1], r_obs[2] - rRef[2]};
        const std::complex<double> temp(exp(-I*k * (kHat[0] * r_obs_rRef[0] + kHat[1] * r_obs_rRef[1] + kHat[2] * r_obs_rRef[2]) ));
        for (int m=0 ; m<3 ; m++) E_inc_i[m] = E0[m] * temp * weigths[j];
        ITo_E_inc[0] += E_inc_i[0];
        ITo_E_inc[1] += E_inc_i[1];
        ITo_E_inc[2] += E_inc_i[2];
        ITo_r_dot_E_inc += r_obs[0] * E_inc_i[0] + r_obs[1] * E_inc_i[1] + r_obs[2] * E_inc_i[2];
        ITo_n_hat_X_r_dot_E_inc += n_hat_X_r[0] * E_inc_i[0] + n_hat_X_r[1] * E_inc_i[1] + n_hat_X_r[2] * E_inc_i[2];

        // computation of ITo_H_inc due to a dipole located at r_dip
        for (int m=0 ; m<3 ; m++) H_inc_i[m] = H0[m] * temp * weigths[j];
        ITo_H_inc[0] += H_inc_i[0];
        ITo_H_inc[1] += H_inc_i[1];
        ITo_H_inc[2] += H_inc_i[2];
        ITo_r_dot_H_inc += r_obs[0] * H_inc_i[0] + r_obs[1] * H_inc_i[1] + r_obs[2] * H_inc_i[2];
        ITo_n_hat_X_r_dot_H_inc += n_hat_X_r[0] * H_inc_i[0] + n_hat_X_r[1] * H_inc_i[1] + n_hat_X_r[2] * H_inc_i[2];
      }
      double norm_factor = triangle.A/sum_weigths;

      for (int i=0 ; i<3 ; ++i) ITo_E_inc[i] *= norm_factor;
      for (int i=0 ; i<3 ; ++i) ITo_H_inc[i] *= norm_factor;
      ITo_r_dot_E_inc *= norm_factor;
      ITo_r_dot_H_inc *= norm_factor;
      ITo_n_hat_X_r_dot_E_inc *= norm_factor;
      ITo_n_hat_X_r_dot_H_inc *= norm_factor;

      const int local_number_edge_p = numbers_RWG_test(rwg);
      const int sign_edge_p = (tr==0) ? 1 : -1;
      const double C_rp = sign_edge_p * l_p * 0.5/triangle.A;
      double *r_p;
      r_p = r_opp;
      double n_hat_X_r_p[3];
      cross3D(n_hat_X_r_p, triangle.n_hat, r_p);

      complex<double> tmpResult(0.0, 0.0);
      tmpResult -= tE * C_rp * (ITo_r_dot_E_inc - (r_p[0]*ITo_E_inc[0] + r_p[1]*ITo_E_inc[1] + r_p[2]*ITo_E_inc[2])); // -<f_m ; E_inc>
      if (RWGNumber_CFIE_OK(local_number_edge_p) == 1) {
        tmpResult -= nE * C_rp * (ITo_n_hat_X_r_dot_E_inc - (n_hat_X_r_p[0]*ITo_E_inc[0] + n_hat_X_r_p[1]*ITo_E_inc[1] + n_hat_X_r_p[2]*ITo_E_inc[2])); // -<n_hat x f_m ; E_inc> 
        tmpResult -= tH * C_rp * (ITo_r_dot_H_inc - (r_p[0]*ITo_H_inc[0] + r_p[1]*ITo_H_inc[1] + r_p[2]*ITo_H_inc[2])); // -<f_m ; H_inc>
        tmpResult -= nH * C_rp * (ITo_n_hat_X_r_dot_H_inc - (n_hat_X_r_p[0]*ITo_H_inc[0] + n_hat_X_r_p[1]*ITo_H_inc[1] + n_hat_X_r_p[2]*ITo_H_inc[2])); // -<n_hat x f_m ; H_inc> 
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


// void V_E_V_H_horn_BBHA (Array<complex<double>,2>& V_EJ, Array<complex<double>,2>& V_HJ, Array<complex<double>,2>& V_EM, Array<complex<double>,2>& V_HM, const mesh & MESH_TARGET, const mesh & MESH_ANTENNA, const Vector<double,3>& r_ant, const Vector<double,3>& x_hat_ant, const Vector<double,3>& y_hat_ant, const G_EJ_grid & G_EJ_tab_ant_object, const G_HJ_grid & G_HJ_tab_ant_object, const layers_constants & LC, const G_EJ_grid & G_EJ_tab_ant_object_dual, const G_HJ_grid & G_HJ_tab_ant_object_dual, const layers_constants & LC_dual) {

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
//     Vector<double, 3> r_j, r_J_on_antenna;
//     Vector<complex<double>, 3> J_ant, M_ant;
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
