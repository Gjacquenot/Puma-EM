#include <iostream>
#include <complex>
#include <blitz/array.h>
#include <vector>
#include <algorithm>

using namespace std;

#include "EMConstants.h"
#include "GK_triangle.h"
#include "triangle_int.h"
#include "mesh.h"
#include "V_E_V_H.h"

void V_CFIE_slot (blitz::Array<std::complex<float>, 1> V_CFIE,
                  const blitz::Array<std::complex<float>, 1>& CFIE,
                  const std::complex<double>& E_0,
                  const blitz::Array<double, 1>& l_hat,
                  const blitz::Array<double, 1>& r_ref,
                  const double slot_length,
                  const blitz::Array<int, 1>& numbers_RWG_test,
                  const blitz::Array<int, 1>& RWGNumber_CFIE_OK,
                  const blitz::Array<float, 2>& RWGNumber_trianglesCoord,
                  const double w,
                  const std::complex<double>& eps_r,
                  const std::complex<double>& mu_r,
                  const int FULL_PRECISION)
/**
 * This function computes the CFIE excitation vectors of the MoM due to a
 * slot which center is located at r_ref, having length slot_length, 
 * directed following l_hat, whose electric field has amplitude E_0.
 *
 */
{
  // def of k, mu_i, eps_i
  blitz::Range all = blitz::Range::all();
  int N_RWG_test = numbers_RWG_test.size();
  std::complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu);
  const complex<double> tE = CFIE(0), nE = CFIE(1), tH = CFIE(2), nH = CFIE(3);

  // geometrical entities
  double r0[3], r1[3], r2[3], *r_opp;
  std::complex<double> ITo_r_dot_H_inc, ITo_r_dot_E_inc, ITo_n_hat_X_r_dot_H_inc, ITo_n_hat_X_r_dot_E_inc;
  std::complex<double> H_inc_i[3], ITo_H_inc[3], E_inc_i[3], ITo_E_inc[3];
  blitz::Array<std::complex<double>, 2> G_EJ (3, 3), G_HJ (3, 3);
  blitz::Array<std::complex<double>, 2> G_EM (3, 3), G_HM (3, 3);
  double rRef[3];
  std::complex<double> E0[3];
  for (int i=0 ; i<3 ; ++i) rRef[i] = r_ref(i);
  for (int i=0 ; i<3 ; ++i) E0[i] = E_0(i);

  V_CFIE = 0.0;

  for (int rwg=0 ; rwg<N_RWG_test ; ++rwg) { // loop on the RWGs
    for (int tr = 0 ; tr<2 ; ++tr) {
      double l_p;
      blitz::Array<double, 1> rt0(3), rt1(3), rt2(3), r_opp(3);
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
      int N_points = 6;
      if (IS_NEAR) N_points = 9;
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
        for (int m=0 ; m<3 ; m++) E_inc_i[m] = (G_EJ (m, 0) * JDip[0] + G_EJ (m, 1) * JDip[1] + G_EJ (m, 2) * JDip[2]) * weigths[j];
        ITo_E_inc[0] += E_inc_i[0];
        ITo_E_inc[1] += E_inc_i[1];
        ITo_E_inc[2] += E_inc_i[2];
        ITo_r_dot_E_inc += (r_obs[0] * E_inc_i[0] + r_obs[1] * E_inc_i[1] + r_obs[2] * E_inc_i[2]);
        ITo_n_hat_X_r_dot_E_inc += (n_hat_X_r[0] * E_inc_i[0] + n_hat_X_r[1] * E_inc_i[1] + n_hat_X_r[2] * E_inc_i[2]);

        // computation of ITo_H_inc due to a dipole located at r_dip
        for (int m=0 ; m<3 ; m++) H_inc_i[m] = (G_HJ (m, 0) * JDip[0] + G_HJ (m, 1) * JDip[1] + G_HJ (m, 2) * JDip[2]) * weigths[j];
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

      std::complex<double> tmpResult(0.0, 0.0);
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

void local_V_CFIE_slot (blitz::Array<std::complex<float>, 1>& V_CFIE,
                        const std::complex<double> E_0,
                        const blitz::Array<double, 1>& l_hat,
                        const blitz::Array<double, 1>& r_ref,
                        const double slot_length,
                        const LocalMesh & local_target_mesh,
                        const double w,
                        const std::complex<double>& eps_r,
                        const std::complex<double>& mu_r,
                        const blitz::Array<std::complex<float>, 1>& CFIE,
                        const int FULL_PRECISION)
{
  // slot excitation vector
  V_CFIE.resize(local_target_mesh.N_local_RWG);
  V_CFIE_slot (V_CFIE, CFIE, E_0, l_hat, r_ref, slot_length, local_target_mesh.reallyLocalRWGNumbers, local_target_mesh.localRWGNumber_CFIE_OK, local_target_mesh.localRWGNumber_trianglesCoord, w, eps_r, mu_r, FULL_PRECISION);
}

