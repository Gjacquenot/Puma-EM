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

void G_EJ_G_HJ (std::vector<std::vector< std::complex<double> > >& G_EJ,
                std::vector<std::vector< std::complex<double> > >& G_HJ,
                const double r_dip[],
                const double r_obs[],
                const std::complex<double>& eps,
                const std::complex<double>& mu,
                const std::complex<double>& k)
/**
 * This function computes the homogeneous space Green's functions
 * due to an elementary electric current element.
 *
 * By reciprocity, we have that G_EM = -G_HJ and G_HM = eps/mu * G_EJ.
 */
{
  const double r_obs_r_dip[3] = {r_obs[0]-r_dip[0], r_obs[1]-r_dip[1], r_obs[2]-r_dip[2]};
  const double R = sqrt(dot3D(r_obs_r_dip, r_obs_r_dip));
  const std::complex<double> kRsquare = k*k*R*R;
  const std::complex<double> term_1 = 1.0 + 1.0/(I*k*R);
  const std::complex<double> term_2 = 1.0/R * term_1 + I*k/2.0 * (term_1 - 1.0/kRsquare);
  const std::complex<double> exp_ikR = exp(-I*k*R), exp_ikR_R = sqrt(mu/eps)/(2.0*M_PI) * exp_ikR/R;
  const double x_xp = r_obs_r_dip[0], y_yp = r_obs_r_dip[1], z_zp = r_obs_r_dip[2];
  const double ONE_R_R = 1.0/(R*R);
  const double x_xp_R_square = (x_xp*x_xp) * ONE_R_R;
  const double y_yp_R_square = (y_yp*y_yp) * ONE_R_R;
  const double z_zp_R_square = (z_zp*z_zp) * ONE_R_R;
  G_EJ [0][1] = term_2*y_yp/R*x_xp/R * exp_ikR_R; 
  G_EJ [0][2] = term_2*z_zp/R*x_xp/R * exp_ikR_R; 
  G_EJ [1][2] = term_2*z_zp/R*y_yp/R * exp_ikR_R;
  G_EJ [1][0] = G_EJ [0][1];
  G_EJ [2][0] = G_EJ [0][2];
  G_EJ [2][1] = G_EJ [1][2];
  G_EJ [0][0] = (x_xp_R_square*term_1/R-(1.0-x_xp_R_square)*I*k/2.0*(term_1 - 1.0/kRsquare)) * exp_ikR_R;
  G_EJ [1][1] = (y_yp_R_square*term_1/R-(1.0-y_yp_R_square)*I*k/2.0*(term_1 - 1.0/kRsquare)) * exp_ikR_R;
  G_EJ [2][2] = (z_zp_R_square*term_1/R-(1.0-z_zp_R_square)*I*k/2.0*(term_1 - 1.0/kRsquare)) * exp_ikR_R;

  std::complex<double> G_i = exp_ikR/(4.0*M_PI) * (1.0+I*k*R)/(R*R*R);
  G_HJ [0][1] = (z_zp) * G_i;
  G_HJ [1][0] = -G_HJ [0][1];
  G_HJ [2][0] = (y_yp) * G_i;
  G_HJ [2][1] = -1.0*(x_xp) * G_i;
  G_HJ [0][2] = -G_HJ [2][0];
  G_HJ [1][2] = -G_HJ [2][1];
  for (int i=0 ; i<3 ; i++) G_HJ[i][i] = 0.0;
}

void V_EJ_HJ_dipole (std::vector<std::complex<double> >& V_tE_J,
                     std::vector<std::complex<double> >& V_nE_J,
                     std::vector<std::complex<double> >& V_tH_J,
                     std::vector<std::complex<double> >& V_nH_J,
                     const std::complex<double> JDip[],
                     const double rDip[],
                     const blitz::Array<int, 1>& numbers_RWG_test,
                     const blitz::Array<int, 1>& RWGNumber_CFIE_OK,
                     const blitz::Array<int, 2>& RWGNumber_signedTriangles,
                     const blitz::Array<double, 2>& RWGNumber_vertexesCoord,
                     const blitz::Array<double, 2>& RWGNumber_oppVertexesCoord,
                     const double w,
                     const std::complex<double>& eps_r,
                     const std::complex<double>& mu_r,
                     const int FULL_PRECISION)
/**
 * This function computes the 4 excitation vectors of the MoM due to an
 * elementary electric current element JDip located at rDip.
 *
 * By reciprocity, we have that V_EM = -V_HJ and V_HM = eps/mu * V_EJ.
 */
{
  // def of k, mu_i, eps_i
  int N_RWG_test = numbers_RWG_test.size();
  std::complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu);
  // triangle integration precision. Possible values for N_points: 1, 3, 6, 9, 12, 13
  int N_points_far, N_points_near;
  if (FULL_PRECISION!=0) {
    N_points_far = 6; 
    N_points_near = 9;
  }
  else {
    N_points_far = 1; 
    N_points_near = 3;
  }

  for (unsigned int i=0; i<V_tE_J.size(); i++) {
    V_tE_J[i] = 0.0; V_tH_J[i] = 0.0;
    V_nE_J[i] = 0.0; V_nH_J[i] = 0.0;
  }
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
  const double *r0, *r1, *r2, *rGrav;
  std::complex<double> ITo_r_dot_H_inc, ITo_r_dot_E_inc, ITo_n_hat_X_r_dot_H_inc, ITo_n_hat_X_r_dot_E_inc;
  std::complex<double> H_inc_i[3], ITo_H_inc[3], E_inc_i[3], ITo_E_inc[3];
  std::vector< std::vector < std::complex<double> > > G_EJ, G_HJ;
  G_EJ.resize(3);
  G_HJ.resize(3);
  for (int i=0; i<3; i++) {
    G_EJ[i].resize(3);
    G_HJ[i].resize(3);
  }


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

      r0 = triangles_test[r].r_nodes[0];
      r1 = triangles_test[r].r_nodes[1];
      r2 = triangles_test[r].r_nodes[2];
      rGrav = triangles_test[r].r_grav;

      // weights and abscissas for triangle integration
      const double rGrav_rDip[3] = {rGrav[0] - rDip[0], rGrav[1] - rDip[1], rGrav[2] - rDip[2]};
      const double R_os = sqrt(dot3D(rGrav_rDip, rGrav_rDip));
      const bool IS_NEAR = (R_os - 1.5 * triangles_test[r].R_max <= 0.0);
      int N_points = N_points_far;
      if (IS_NEAR) N_points = N_points_near;
      double sum_weigths;
      const double *xi, *eta, *weigths;
      IT_points (xi, eta, weigths, sum_weigths, N_points);
      // triangle integration
      for (int j=0 ; j<N_points ; ++j) {
        double r_obs[3];
        r_obs[0] = r0[0] * xi[j] + r1[0] * eta[j] + r2[0] * (1-xi[j]-eta[j]);
        r_obs[1] = r0[1] * xi[j] + r1[1] * eta[j] + r2[1] * (1-xi[j]-eta[j]);
        r_obs[2] = r0[2] * xi[j] + r1[2] * eta[j] + r2[2] * (1-xi[j]-eta[j]);
        double n_hat_X_r[3];
        cross3D(n_hat_X_r, triangles_test[r].n_hat, r_obs);
        G_EJ_G_HJ (G_EJ, G_HJ, rDip, r_obs, eps, mu, k);

        // computation of ITo_E_inc due to a dipole located at r_dip
        for (int m=0 ; m<3 ; m++) E_inc_i[m] = (G_EJ [m][0] * JDip[0] + G_EJ [m][1] * JDip[1] + G_EJ [m][2] * JDip[2]) * weigths[j];
        ITo_E_inc[0] += E_inc_i[0];
        ITo_E_inc[1] += E_inc_i[1];
        ITo_E_inc[2] += E_inc_i[2];
        ITo_r_dot_E_inc += r_obs[0] * E_inc_i[0] + r_obs[1] * E_inc_i[1] + r_obs[2] * E_inc_i[2];
        ITo_n_hat_X_r_dot_E_inc += n_hat_X_r[0] * E_inc_i[0] + n_hat_X_r[1] * E_inc_i[1] + n_hat_X_r[2] * E_inc_i[2];

        // computation of ITo_H_inc due to a dipole located at r_dip
        for (int m=0 ; m<3 ; m++) H_inc_i[m] = (G_HJ [m][0] * JDip[0] + G_HJ [m][1] * JDip[1] + G_HJ [m][2] * JDip[2]) * weigths[j];
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

        V_tE_J [local_number_edge_p] += -C_rp * (ITo_r_dot_E_inc - (r_p[0]*ITo_E_inc[0] + r_p[1]*ITo_E_inc[1] + r_p[2]*ITo_E_inc[2])); // -<f_m ; E_inc> 
        V_tH_J [local_number_edge_p] += -C_rp * (ITo_r_dot_H_inc - (r_p[0]*ITo_H_inc[0] + r_p[1]*ITo_H_inc[1] + r_p[2]*ITo_H_inc[2])); // -<f_m ; H_inc>
        V_nE_J [local_number_edge_p] += -C_rp * (ITo_n_hat_X_r_dot_E_inc - (n_hat_X_r_p[0]*ITo_E_inc[0] + n_hat_X_r_p[1]*ITo_E_inc[1] + n_hat_X_r_p[2]*ITo_E_inc[2])); // -<n_hat x f_m ; E_inc> 
        V_nH_J [local_number_edge_p] += -C_rp * (ITo_n_hat_X_r_dot_H_inc - (n_hat_X_r_p[0]*ITo_H_inc[0] + n_hat_X_r_p[1]*ITo_H_inc[1] + n_hat_X_r_p[2]*ITo_H_inc[2])); // -<n_hat x f_m ; H_inc> 
      }
  }
}

void V_CFIE_dipole (blitz::Array<std::complex<float>, 1> V_CFIE,
                    const blitz::Array<std::complex<float>, 1>& CFIE,
                    const blitz::Array<std::complex<double>, 1>& J_dip,
                    const blitz::Array<double, 1>& r_dip,
                    const blitz::Array<int, 1>& numbers_RWG_test,
                    const blitz::Array<int, 1>& RWGNumber_CFIE_OK,
                    const blitz::Array<float, 2>& RWGNumber_trianglesCoord,
                    const double w,
                    const std::complex<double>& eps_r,
                    const std::complex<double>& mu_r,
                    const int FULL_PRECISION)
/**
 * This function computes the CFIE excitation vectors of the MoM due to an
 * elementary electric current element J_dip located at r_dip.
 *
 */
{
  // def of k, mu_i, eps_i
  int N_RWG_test = numbers_RWG_test.size();
  std::complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu);
  const complex<double> tE = CFIE(0), nE = CFIE(1), tH = CFIE(2), nH = CFIE(3);
  // triangle integration precision. Possible values for N_points: 1, 3, 6, 9, 12, 13
  int N_points_far, N_points_near;
  if (FULL_PRECISION!=0) {
    N_points_far = 6; 
    N_points_near = 9;
  }
  else {
    N_points_far = 1; 
    N_points_near = 3;
  }

  // geometrical entities
  double r0[3], r1[3], r2[3], *r_opp;
  std::complex<double> ITo_r_dot_H_inc, ITo_r_dot_E_inc, ITo_n_hat_X_r_dot_H_inc, ITo_n_hat_X_r_dot_E_inc;
  std::complex<double> H_inc_i[3], ITo_H_inc[3], E_inc_i[3], ITo_E_inc[3];
  std::vector< std::vector < std::complex<double> > > G_EJ, G_HJ;
  G_EJ.resize(3);
  G_HJ.resize(3);
  for (int i=0; i<3; i++) {
    G_EJ[i].resize(3);
    G_HJ[i].resize(3);
  }
  double rDip[3];
  std::complex<double> JDip[3];
  for (int i=0 ; i<3 ; ++i) rDip[i] = r_dip(i);
  for (int i=0 ; i<3 ; ++i) JDip[i] = J_dip(i);

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
        for (int i=0; i<3; i++) { 
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

      // weights and abscissas for triangle integration
      const double rGrav_rDip[3] = {triangle.r_grav[0] - rDip[0], triangle.r_grav[1] - rDip[1], triangle.r_grav[2] - rDip[2]};
      double R_os = sqrt (dot3D(rGrav_rDip, rGrav_rDip));
      bool IS_NEAR = (R_os - 1.5 * triangle.R_max <= 0.0);
      int N_points = N_points_far;
      if (IS_NEAR) N_points = N_points_near;
      double sum_weigths;
      const double *xi, *eta, *weigths;
      IT_points (xi, eta, weigths, sum_weigths, N_points);
      // triangle integration
      for (int j=0 ; j<N_points ; ++j) {
        double r_obs[3];
        r_obs[0] = r0[0] * xi[j] + r1[0] * eta[j] + r2[0] * (1-xi[j]-eta[j]);
        r_obs[1] = r0[1] * xi[j] + r1[1] * eta[j] + r2[1] * (1-xi[j]-eta[j]);
        r_obs[2] = r0[2] * xi[j] + r1[2] * eta[j] + r2[2] * (1-xi[j]-eta[j]);
        double n_hat_X_r[3];
        cross3D(n_hat_X_r, triangle.n_hat, r_obs);
        G_EJ_G_HJ (G_EJ, G_HJ, rDip, r_obs, eps, mu, k);

        // computation of ITo_E_inc due to a dipole located at r_dip
        for (int m=0 ; m<3 ; m++) E_inc_i[m] = (G_EJ [m][0] * JDip[0] + G_EJ [m][1] * JDip[1] + G_EJ [m][2] * JDip[2]) * weigths[j];
        ITo_E_inc[0] += E_inc_i[0];
        ITo_E_inc[1] += E_inc_i[1];
        ITo_E_inc[2] += E_inc_i[2];
        ITo_r_dot_E_inc += (r_obs[0] * E_inc_i[0] + r_obs[1] * E_inc_i[1] + r_obs[2] * E_inc_i[2]);
        ITo_n_hat_X_r_dot_E_inc += (n_hat_X_r[0] * E_inc_i[0] + n_hat_X_r[1] * E_inc_i[1] + n_hat_X_r[2] * E_inc_i[2]);

        // computation of ITo_H_inc due to a dipole located at r_dip
        for (int m=0 ; m<3 ; m++) H_inc_i[m] = (G_HJ [m][0] * JDip[0] + G_HJ [m][1] * JDip[1] + G_HJ [m][2] * JDip[2]) * weigths[j];
        ITo_H_inc[0] += H_inc_i[0];
        ITo_H_inc[1] += H_inc_i[1];
        ITo_H_inc[2] += H_inc_i[2];
        ITo_r_dot_H_inc += (r_obs[0] * H_inc_i[0] + r_obs[1] * H_inc_i[1] + r_obs[2] * H_inc_i[2]);
        ITo_n_hat_X_r_dot_H_inc += (n_hat_X_r[0] * H_inc_i[0] + n_hat_X_r[1] * H_inc_i[1] + n_hat_X_r[2] * H_inc_i[2]);
      }
      const double norm_factor = triangle.A/sum_weigths;

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
      tmpResult -= nE * C_rp * (ITo_n_hat_X_r_dot_E_inc - (n_hat_X_r_p[0]*ITo_E_inc[0] + n_hat_X_r_p[1]*ITo_E_inc[1] + n_hat_X_r_p[2]*ITo_E_inc[2])); // -<n_hat x f_m ; E_inc> 
      if (RWGNumber_CFIE_OK(local_number_edge_p) == 1) {
        tmpResult -= tH * C_rp * (ITo_r_dot_H_inc - (r_p[0]*ITo_H_inc[0] + r_p[1]*ITo_H_inc[1] + r_p[2]*ITo_H_inc[2])); // -<f_m ; H_inc>
        tmpResult -= nH * C_rp * (ITo_n_hat_X_r_dot_H_inc - (n_hat_X_r_p[0]*ITo_H_inc[0] + n_hat_X_r_p[1]*ITo_H_inc[1] + n_hat_X_r_p[2]*ITo_H_inc[2])); // -<n_hat x f_m ; H_inc> 
      }
      V_CFIE(local_number_edge_p) += tmpResult;
    }
  }
}

void local_V_CFIE_dipole (blitz::Array<std::complex<float>, 1>& V_CFIE,
                          const blitz::Array<std::complex<double>, 1>& J_dip,
                          const blitz::Array<double, 1>& r_dip,
                          const LocalMesh & local_target_mesh,
                          const double w,
                          const std::complex<double>& eps_r,
                          const std::complex<double>& mu_r,
                          const blitz::Array<std::complex<float>, 1>& CFIE,
                          const int FULL_PRECISION)
{
  // We now compute the excitation vectors
  V_CFIE.resize(local_target_mesh.N_local_RWG);
  V_CFIE_dipole (V_CFIE, CFIE, J_dip, r_dip, local_target_mesh.reallyLocalRWGNumbers, local_target_mesh.localRWGNumber_CFIE_OK, local_target_mesh.localRWGNumber_trianglesCoord, w, eps_r, mu_r, FULL_PRECISION);
}


void V_CFIE_dipole_array (blitz::Array<std::complex<float>, 1> V_CFIE,
                          const blitz::Array<std::complex<float>, 1>& CFIE,
                          const blitz::Array<std::complex<double>, 2>& J_dip,
                          const blitz::Array<double, 2>& r_dip,
                          const blitz::Array<int, 1>& numbers_RWG_test,
                          const blitz::Array<int, 1>& RWGNumber_CFIE_OK,
                          const blitz::Array<float, 2>& RWGNumber_trianglesCoord,
                          const double w,
                          const std::complex<double>& eps_r,
                          const std::complex<double>& mu_r,
                          const char CURRENT_TYPE,
                          const int FULL_PRECISION)
/**
 * This function computes the CFIE excitation vectors of the MoM due to
 * multiple elementary electric current elements J_dip located at r_dip.
 * This allows for more complex sources to be modeled.
 *
 */
{
  // def of k, mu_i, eps_i
  int N_RWG_test = numbers_RWG_test.size(), N_dipoles = J_dip.rows();
  if (r_dip.size() != J_dip.size()) {
    cout << "Error in V_CFIE_dipole_array: J_dip and R_dip don't have the same size. Aborting." << endl;
    exit(1);
  }
  std::complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu);
  const complex<double> tE = CFIE(0), nE = CFIE(1), tH = CFIE(2), nH = CFIE(3);
  // triangle integration precision. Possible values for N_points: 1, 3, 6, 9, 12, 13
  int N_points_far, N_points_near;
  if (FULL_PRECISION!=0) {
    N_points_far = 6; 
    N_points_near = 9;
  }
  else {
    N_points_far = 1; 
    N_points_near = 3;
  }

  // geometrical entities
  double r0[3], r1[3], r2[3], *r_opp;
  std::complex<double> ITo_r_dot_H_inc, ITo_r_dot_E_inc, ITo_n_hat_X_r_dot_H_inc, ITo_n_hat_X_r_dot_E_inc;
  std::complex<double> H_inc_i[3], ITo_H_inc[3], E_inc_i[3], ITo_E_inc[3];
  std::vector< std::vector < std::complex<double> > > G_EJ, G_HJ;
  G_EJ.resize(3);
  G_HJ.resize(3);
  for (int i=0; i<3; i++) {
    G_EJ[i].resize(3);
    G_HJ[i].resize(3);
  }
  double rDip[3];
  std::complex<double> JDip[3];
  for (int i=0 ; i<3 ; ++i) rDip[i] = r_dip(i);
  for (int i=0 ; i<3 ; ++i) JDip[i] = J_dip(i);

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
        for (int i=0; i<3; i++) { 
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

      // we first select the weights and abscissas for triangle integration
      bool IS_NEAR;
      int N_points = N_points_far;
      for (int dipNumber = 0 ; dipNumber < N_dipoles ; dipNumber++) {
        for (int i=0 ; i<3 ; ++i) rDip[i] = r_dip(dipNumber, i);
        const double rGrav_rDip[3] = {triangle.r_grav[0] - rDip[0], triangle.r_grav[1] - rDip[1], triangle.r_grav[2] - rDip[2]};
        double R_os = sqrt (dot3D(rGrav_rDip, rGrav_rDip));
        IS_NEAR = (R_os - 1.5 * triangle.R_max <= 0.0);
        if (IS_NEAR) N_points = N_points_near;
        if (IS_NEAR) break;
      }
      double sum_weigths;
      const double *xi, *eta, *weigths;
      IT_points (xi, eta, weigths, sum_weigths, N_points);
      // now the loop on the sources
      for (int dipNumber = 0 ; dipNumber < N_dipoles ; dipNumber++) {
        for (int i=0 ; i<3 ; ++i) rDip[i] = r_dip(dipNumber, i);
        for (int i=0 ; i<3 ; ++i) JDip[i] = J_dip(dipNumber, i);
        // triangle integration
        for (int j=0 ; j<N_points ; ++j) {
          double r_obs[3];
          r_obs[0] = r0[0] * xi[j] + r1[0] * eta[j] + r2[0] * (1-xi[j]-eta[j]);
          r_obs[1] = r0[1] * xi[j] + r1[1] * eta[j] + r2[1] * (1-xi[j]-eta[j]);
          r_obs[2] = r0[2] * xi[j] + r1[2] * eta[j] + r2[2] * (1-xi[j]-eta[j]);
          double n_hat_X_r[3];
          cross3D(n_hat_X_r, triangle.n_hat, r_obs);
          G_EJ_G_HJ (G_EJ, G_HJ, rDip, r_obs, eps, mu, k);
          // we check if we have electric or magnetic dipole 
          if (CURRENT_TYPE == 'M') {
            // then G_EJ should be filled with G_EM values
            // and G_HJ with G_HM values. Since G_EM = -G_HJ, and G_HM = eps/mu * G_EJ
            // it is a simple interchange and scaling of G_EJ and G_HJ values.
            std::vector<std::vector<std::complex<double> > > G_tmp;
            G_tmp.resize(3);
            for (int i=0; i<3; i++) G_tmp[i].resize(3);
            for (int i=0; i<3; i++) {
              for (int ii=0; ii<3; ii++) G_tmp[i][ii] = G_EJ[i][ii]; // we "save" G_EJ values;
            }
            for (int i=0; i<3; i++) {
              for (int ii=0; ii<3; ii++) G_EJ[i][ii] = -G_HJ[i][ii]; // G_EJ becomes G_EM
            }
            for (int i=0; i<3; i++) {
              for (int ii=0; ii<3; ii++) G_HJ[i][ii] = eps/mu * G_tmp[i][ii]; // G_HJ becomes G_HM
            }
          }

          // computation of ITo_E_inc due to a dipole located at r_dip
          for (int m=0 ; m<3 ; m++) E_inc_i[m] = (G_EJ [m][0] * JDip[0] + G_EJ [m][1] * JDip[1] + G_EJ [m][2] * JDip[2]) * weigths[j];
          ITo_E_inc[0] += E_inc_i[0];
          ITo_E_inc[1] += E_inc_i[1];
          ITo_E_inc[2] += E_inc_i[2];
          ITo_r_dot_E_inc += (r_obs[0] * E_inc_i[0] + r_obs[1] * E_inc_i[1] + r_obs[2] * E_inc_i[2]);
          ITo_n_hat_X_r_dot_E_inc += (n_hat_X_r[0] * E_inc_i[0] + n_hat_X_r[1] * E_inc_i[1] + n_hat_X_r[2] * E_inc_i[2]);

          // computation of ITo_H_inc due to a dipole located at r_dip
          for (int m=0 ; m<3 ; m++) H_inc_i[m] = (G_HJ [m][0] * JDip[0] + G_HJ [m][1] * JDip[1] + G_HJ [m][2] * JDip[2]) * weigths[j];
          ITo_H_inc[0] += H_inc_i[0];
          ITo_H_inc[1] += H_inc_i[1];
          ITo_H_inc[2] += H_inc_i[2];
          ITo_r_dot_H_inc += (r_obs[0] * H_inc_i[0] + r_obs[1] * H_inc_i[1] + r_obs[2] * H_inc_i[2]);
          ITo_n_hat_X_r_dot_H_inc += (n_hat_X_r[0] * H_inc_i[0] + n_hat_X_r[1] * H_inc_i[1] + n_hat_X_r[2] * H_inc_i[2]);

        } // loop on the triangle integration points
      } // loop on the dipoles
 
      const double norm_factor = triangle.A/sum_weigths;

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

      std::complex<double> tmpResult(0.0, 0.0);
      tmpResult -= tE * C_rp * (ITo_r_dot_E_inc - (r_p[0]*ITo_E_inc[0] + r_p[1]*ITo_E_inc[1] + r_p[2]*ITo_E_inc[2])); // -<f_m ; E_inc>
      tmpResult -= nE * C_rp * (ITo_n_hat_X_r_dot_E_inc - (n_hat_X_r_p[0]*ITo_E_inc[0] + n_hat_X_r_p[1]*ITo_E_inc[1] + n_hat_X_r_p[2]*ITo_E_inc[2])); // -<n_hat x f_m ; E_inc> 
      if (RWGNumber_CFIE_OK(local_number_edge_p) == 1) {
        tmpResult -= tH * C_rp * (ITo_r_dot_H_inc - (r_p[0]*ITo_H_inc[0] + r_p[1]*ITo_H_inc[1] + r_p[2]*ITo_H_inc[2])); // -<f_m ; H_inc>
        tmpResult -= nH * C_rp * (ITo_n_hat_X_r_dot_H_inc - (n_hat_X_r_p[0]*ITo_H_inc[0] + n_hat_X_r_p[1]*ITo_H_inc[1] + n_hat_X_r_p[2]*ITo_H_inc[2])); // -<n_hat x f_m ; H_inc> 
      }
      V_CFIE(local_number_edge_p) += tmpResult;
    }
  }
}

void local_V_CFIE_dipole_array (blitz::Array<std::complex<float>, 1>& V_CFIE,
                                const blitz::Array<std::complex<double>, 2>& J_dip,
                                const blitz::Array<double, 2>& r_dip,
                                const LocalMesh & local_target_mesh,
                                const double w,
                                const std::complex<double>& eps_r,
                                const std::complex<double>& mu_r,
                                const blitz::Array<std::complex<float>, 1>& CFIE,
                                const char CURRENT_TYPE,
                                const int FULL_PRECISION)
{
  // We now compute the excitation vectors
  V_CFIE.resize(local_target_mesh.N_local_RWG);
  V_CFIE_dipole_array (V_CFIE, CFIE, J_dip, r_dip, local_target_mesh.reallyLocalRWGNumbers, local_target_mesh.localRWGNumber_CFIE_OK, local_target_mesh.localRWGNumber_trianglesCoord, w, eps_r, mu_r, CURRENT_TYPE, FULL_PRECISION);
}

/*****************************************
 * computation of the observation fields *
 *****************************************/
void compute_E_obs(blitz::Array<std::complex<float>, 1>& E_obs,
                   blitz::Array<std::complex<float>, 1>& H_obs,
                   const blitz::Array<double, 1>& r_obs,
                   const blitz::Array<std::complex<float>, 1>& ZI,
                   const blitz::Array<int, 1>& numbers_RWG_test,
                   const blitz::Array<float, 2>& RWGNumber_trianglesCoord,
                   const double w,
                   const std::complex<double>& eps_r,
                   const std::complex<double>& mu_r,
                   const int FULL_PRECISION)
{
  // def of k, mu_i, eps_i
  int N_RWG_test = numbers_RWG_test.size();
  std::complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu);
  // triangle integration precision. Possible values for N_points: 1, 3, 6, 9, 12, 13
  int N_points_far, N_points_near;
  if (FULL_PRECISION!=0) {
    N_points_far = 6; 
    N_points_near = 9;
  }
  else {
    N_points_far = 1; 
    N_points_near = 3;
  }

  // geometrical entities
  double r0[3], r1[3], r2[3], *r_opp;
  std::vector< std::vector < std::complex<double> > > G_EJ, G_HJ;
  G_EJ.resize(3);
  G_HJ.resize(3);
  for (int i=0; i<3; i++) {
    G_EJ[i].resize(3);
    G_HJ[i].resize(3);
  }
  double rObs[3];
  for (int i=0 ; i<3 ; ++i) rObs[i] = r_obs(i);

  E_obs = 0.0;   // E_obs = int_S G_EJ dot J + int_S G_EM dot M
  H_obs = 0.0;   // H_obs = int_S G_HJ dot J + int_S G_HM dot M
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
        for (int i=0; i<3; i++) { 
          r0[i] = RWGNumber_trianglesCoord(rwg, i+6);
          r1[i] = RWGNumber_trianglesCoord(rwg, i+3);
          r2[i] = RWGNumber_trianglesCoord(rwg, i+9);
        }
        r_opp = r2;
        double r0_r1[3] = {r0[0]-r1[0], r0[1]-r1[1], r0[2]-r1[2]};
        l_p = sqrt(dot3D(r0_r1, r0_r1));
      }
      Triangle triangle(r0, r1, r2, 0);

      // weights and abscissas for triangle integration
      const double rGrav_rObs[3] = {triangle.r_grav[0] - rObs[0], triangle.r_grav[1] - rObs[1], triangle.r_grav[2] - rObs[2]};
      double R_os = sqrt (dot3D(rGrav_rObs, rGrav_rObs));
      bool IS_NEAR = (R_os - 1.5 * triangle.R_max <= 0.0);
      int N_points = N_points_far;
      if (IS_NEAR) N_points = N_points_near;
      double sum_weigths;
      const double *xi, *eta, *weigths;
      IT_points (xi, eta, weigths, sum_weigths, N_points);
      // triangle integration
      std::complex<double> ITo_E_obs[3] = {0.0, 0.0, 0.0};
      std::complex<double> ITo_H_obs[3] = {0.0, 0.0, 0.0};
      for (int j=0 ; j<N_points ; ++j) {
        double r_src[3];
        r_src[0] = r0[0] * xi[j] + r1[0] * eta[j] + r2[0] * (1-xi[j]-eta[j]);
        r_src[1] = r0[1] * xi[j] + r1[1] * eta[j] + r2[1] * (1-xi[j]-eta[j]);
        r_src[2] = r0[2] * xi[j] + r1[2] * eta[j] + r2[2] * (1-xi[j]-eta[j]);
        G_EJ_G_HJ (G_EJ, G_HJ, rObs, r_src, eps, mu, k);

        // computation of E_obs
        const double fn[3] = {r_src[0]-r_opp[0], r_src[1]-r_opp[1], r_src[2]-r_opp[2]};
        for (int m=0 ; m<3 ; m++) ITo_E_obs[m] += (G_EJ [m][0] * fn[0] + G_EJ [m][1] * fn[1] + G_EJ [m][2] * fn[2]) * weigths[j];

        // computation of H_obs
        for (int m=0 ; m<3 ; m++) ITo_H_obs[m] += (G_HJ [m][0] * fn[0] + G_HJ [m][1] * fn[1] + G_HJ [m][2] * fn[2]) * weigths[j];
      }
      const double norm_factor = triangle.A/sum_weigths;
      for (int i=0 ; i<3 ; ++i) ITo_E_obs[i] *= norm_factor;
      for (int i=0 ; i<3 ; ++i) ITo_H_obs[i] *= norm_factor;

      const double sign_edge_p = (tr==0) ? 1.0 : -1.0;
      const double C_rp = sign_edge_p * l_p * 0.5/triangle.A;
      const int local_number_edge_p = numbers_RWG_test(rwg);
      const std::complex<double> i_rwg = ZI(local_number_edge_p);

      for (int i=0 ; i<3 ; ++i) E_obs(i) += ITo_E_obs[i] * (C_rp*i_rwg);
      for (int i=0 ; i<3 ; ++i) H_obs(i) += ITo_H_obs[i] * (C_rp*i_rwg);
    }
  }
}


void local_compute_E_obs (blitz::Array<std::complex<float>, 1>& E_obs,
                          blitz::Array<std::complex<float>, 1>& H_obs,
                          const blitz::Array<double, 1>& r_obs,
                          const blitz::Array<std::complex<float>, 1>& ZI,
                          const LocalMesh & local_target_mesh,
                          const double w,
                          const std::complex<double>& eps_r,
                          const std::complex<double>& mu_r,
                          const int FULL_PRECISION)
{
  // We now compute the excitation vectors
  E_obs.resize(3);
  H_obs.resize(3);
  compute_E_obs(E_obs, H_obs, r_obs, ZI, local_target_mesh.reallyLocalRWGNumbers, local_target_mesh.localRWGNumber_trianglesCoord, w, eps_r, mu_r, FULL_PRECISION);
}

/*void V_EJ_HJ_dipole_alternative (blitz::Array<std::complex<double>, 1> V_tE_J,
                                 blitz::Array<std::complex<double>, 1> V_nE_J,
                                 blitz::Array<std::complex<double>, 1> V_tH_J,
                                 blitz::Array<std::complex<double>, 1> V_nH_J,
                                 const blitz::Array<std::complex<double>, 1>& J_dip,
                                 const blitz::Array<double, 1>& r_dip,
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
  // This alternative tries to get rid of the singularities associated with a dipole source
  // located near the surface of interest. This function can be used for computing the fields
  // near surface currents distributions by reciprocity.

  // def of k, mu_i, eps_i
  Range all = Range::all();
  int N_RWG_test = numbers_RWG_test.size(), EXTRACT_R, EXTRACT_1_R, N_points;
  complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu);
  // triangle integration precision. Possible values for N_points: 1, 3, 6, 9, 12, 13
  int N_points_far, N_points_near;
  if (FULL_PRECISION!=0) {
    N_points_far = 6; 
    N_points_near = 9;
  }
  else {
    N_points_far = 1; 
    N_points_near = 3;
  }

  // transformation of some input Arrays into Vectors...
  double rDip[3];
  std::complex<double> JDip[3];
  for (int i=0; i<3; i++) rDip[i] = r_dip(i);
  for (int i=0; i<3; i++) JDip[i] = J_dip(i);

  complex<double> ITo_G;
  std::complex<double> ITo_G_r[3], n_hat_X_ITo_G_r[3], ITo_G_r_rDip[3], ITo_grad_G[3], rDip_X_ITo_grad_G[3], r_p_X_ITo_grad_G[3], ITo_n_hat_X_r_X_grad_G[3], n_hat_X_r_p_X_ITo_grad_G[3];

  V_tE_J = 0.0; V_tH_J = 0.0;
  V_nE_J = 0.0; V_nH_J = 0.0;
  vector<HalfRWG> test_halfRWGs;
  for (int i=0 ; i<N_RWG_test ; ++i) {
    for (int j=0 ; j<2 ; ++j) {
      const int RWGnumber = numbers_RWG_test(i);
      const int triangle_number = abs(RWGNumber_signedTriangles(RWGnumber, j)), triangle_sign = (j==0) ? 1 : -1;
      blitz::Array<double, 1> r_opp(3);
      r_opp = (j==0) ? RWGNumber_oppVertexesCoord(RWGnumber, blitz::Range(0,2)) : RWGNumber_oppVertexesCoord(RWGnumber, blitz::Range(3,5));
      const double edge_length = RWGNumber_edgeLength(RWGnumber);
      test_halfRWGs.push_back(HalfRWG(RWGnumber, triangle_number, triangle_sign, edge_length, r_opp));
    }
  }
  sort(test_halfRWGs.begin(), test_halfRWGs.end());
  // vector of triangles construction
  vector<Triangle> triangles_test;
  constructVectorTriangles(triangles_test, test_halfRWGs, testTriangle_vertexesCoord);

  for (int r=0 ; r<triangles_test.size() ; ++r) { // loop on the observation triangles

    const Triangle T(triangles_test[r]);
    const double rGrav_rDip[3] = {triangle.r_grav[0] - rDip[0], triangle.r_grav[1] - rDip[1], triangle.r_grav[2] - rDip[2]};
    double R_os = sqrt (dot3D(rGrav_rDip, rGrav_rDip));
    bool IS_NEAR = (R_os - 1.5 * triangles_test[r].R_max <= 0.0);
    EXTRACT_1_R = (EXTRACT_R = 0); N_points = N_points_far;
    if (IS_NEAR) {
      EXTRACT_1_R = (EXTRACT_R = 1); N_points = N_points_near;
    }
    V_EH_ITo_free (ITo_G, ITo_G_r_rDip, ITo_grad_G, ITo_n_hat_X_r_X_grad_G, rDip, T, k, N_points, EXTRACT_1_R, EXTRACT_R);
    //ITs_free (ITs_G, ITs_G_r_rDip, ITs_grad_G, rDip, T, k, N_points, EXTRACT_1_R, EXTRACT_R);
    ITo_G_r = ITo_G_r_rDip + rDip * ITo_G;
    n_hat_X_ITo_G_r = T.n_hat(1)*ITo_G_r(2)-T.n_hat(2)*ITo_G_r(1),
                      T.n_hat(2)*ITo_G_r(0)-T.n_hat(0)*ITo_G_r(2),
                      T.n_hat(0)*ITo_G_r(1)-T.n_hat(1)*ITo_G_r(0);
    rDip_X_ITo_grad_G = rDip(1)*ITo_grad_G(2)-rDip(2)*ITo_grad_G(1),
                        rDip(2)*ITo_grad_G(0)-rDip(0)*ITo_grad_G(2),
                        rDip(0)*ITo_grad_G(1)-rDip(1)*ITo_grad_G(0);
    
    for (int p=0 ; p<triangles_test[r].halfRWGIndexes.size() ; ++p) {
      const int index_p = triangles_test[r].halfRWGIndexes[p];
      const int number_edge_p = test_halfRWGs[index_p].number;
      const int local_number_edge_p = test_halfRWGs[index_p].number;
      const double l_p = test_halfRWGs[index_p].length;
      const int sign_edge_p = test_halfRWGs[index_p].triangleSign;
      const double C_rp = sign_edge_p * l_p * 0.5/triangles_test[r].A;
      const double r_p[3](test_halfRWGs[index_p].r_opp_vertexesCoord);
      const double n_hat_X_r_p[3](cross(triangles_test[r].n_hat, r_p));
      r_p_X_ITo_grad_G = r_p(1)*ITo_grad_G(2)-r_p(2)*ITo_grad_G(1),
                         r_p(2)*ITo_grad_G(0)-r_p(0)*ITo_grad_G(2),
                         r_p(0)*ITo_grad_G(1)-r_p(1)*ITo_grad_G(0);

      n_hat_X_r_p_X_ITo_grad_G = n_hat_X_r_p(1)*ITo_grad_G(2)-n_hat_X_r_p(2)*ITo_grad_G(1),
                                 n_hat_X_r_p(2)*ITo_grad_G(0)-n_hat_X_r_p(0)*ITo_grad_G(2),
                                 n_hat_X_r_p(0)*ITo_grad_G(1)-n_hat_X_r_p(1)*ITo_grad_G(0);

      V_tE_J (local_number_edge_p) += -1.0/(I*w*eps*4.0*M_PI) * C_rp * dot(JDip, k*k*(ITo_G_r - r_p*ITo_G) + 2.0 * ITo_grad_G); // -<f_m ; E_inc> 
      V_nE_J (local_number_edge_p) += -1.0/(I*w*eps*4.0*M_PI) * C_rp * dot(JDip, k*k*(n_hat_X_ITo_G_r - n_hat_X_r_p*ITo_G) ); // -<n x f_m ; E_inc> 
      V_tH_J (local_number_edge_p) += C_rp/(4.0*M_PI) * dot(JDip, rDip_X_ITo_grad_G - r_p_X_ITo_grad_G); // -<f_m ; H_inc> 
      V_nH_J (local_number_edge_p) += C_rp/(4.0*M_PI) * dot(JDip, ITo_n_hat_X_r_X_grad_G - n_hat_X_r_p_X_ITo_grad_G); // -<f_m ; H_inc> 
    } // end for (p=0 ; p<3 ; p++)
  } // end for (m=0 ; m<N_triangles ; m++)
}
*/

