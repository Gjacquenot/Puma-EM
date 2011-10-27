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
  blitz::TinyVector<double, 3> r0, r1, r2, r_obs;
  std::complex<double> ITo_r_dot_H_inc, ITo_r_dot_E_inc, ITo_n_hat_X_r_dot_H_inc, ITo_n_hat_X_r_dot_E_inc;
  blitz::TinyVector<std::complex<double>, 3> H_inc_i, ITo_H_inc, E_inc_i, ITo_E_inc;
  blitz::Array<std::complex<double>, 2> G_EJ (3, 3), G_HJ (3, 3);
  blitz::Array<std::complex<double>, 2> G_EM (3, 3), G_HM (3, 3);
  blitz::TinyVector<double, 3> rRef;
  blitz::TinyVector<std::complex<double>, 3> E0;
  for (int i=0 ; i<3 ; ++i) rRef(i) = r_ref(i);
  for (int i=0 ; i<3 ; ++i) E0(i) = E_0(i);

  V_CFIE = 0.0;

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
      ITo_E_inc = 0.0; // TinyVector<complex<double>, 3>
      ITo_H_inc = 0.0; // TinyVector<complex<double>, 3>
      ITo_r_dot_E_inc = 0.0; // complex<double>
      ITo_r_dot_H_inc = 0.0; // complex<double>
      ITo_n_hat_X_r_dot_E_inc = 0.0; // complex<double>
      ITo_n_hat_X_r_dot_H_inc = 0.0; // complex<double>

      r0 = triangle.r_nodes(0);
      r1 = triangle.r_nodes(1);
      r2 = triangle.r_nodes(2);
      // weights and abscissas for triangle integration
      double R_os = sqrt (dot (triangle.r_grav - rDip, triangle.r_grav - rDip));
      bool IS_NEAR = (R_os - 1.5 * triangle.R_max <= 0.0);
      int N_points = 6;
      if (IS_NEAR) N_points = 9;
      double sum_weigths;
      const double *xi, *eta, *weigths;
      IT_points (xi, eta, weigths, sum_weigths, N_points);
      // triangle integration
      for (int j=0 ; j<N_points ; ++j) {
        r_obs = r0 * xi[j] + r1 * eta[j] + r2 * (1-xi[j]-eta[j]);
        blitz::TinyVector<double, 3> n_hat_X_r(cross(triangle.n_hat, r_obs));
        G_EJ_G_HJ (G_EJ, G_HJ, rDip, r_obs, eps, mu, k);

        // computation of ITo_E_inc due to a dipole located at r_dip
        for (int m=0 ; m<3 ; m++) E_inc_i (m) = (G_EJ (m, 0) * JDip (0) + G_EJ (m, 1) * JDip (1) + G_EJ (m, 2) * JDip (2)) * weigths[j];
        ITo_E_inc += E_inc_i;
        ITo_r_dot_E_inc += sum(r_obs * E_inc_i);
        ITo_n_hat_X_r_dot_E_inc += sum(n_hat_X_r * E_inc_i);

        // computation of ITo_H_inc due to a dipole located at r_dip
        for (int m=0 ; m<3 ; m++) H_inc_i (m) = (G_HJ (m, 0) * JDip (0) + G_HJ (m, 1) * JDip (1) + G_HJ (m, 2) * JDip (2)) * weigths[j];
        ITo_H_inc += H_inc_i;
        ITo_r_dot_H_inc += dot(r_obs, H_inc_i);
        ITo_n_hat_X_r_dot_H_inc += dot(n_hat_X_r, H_inc_i);
      }
      const double norm_factor = triangle.A/sum_weigths;

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

