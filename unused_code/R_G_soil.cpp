#include <iostream>
#include <complex>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

using namespace blitz;

const complex<double> I (0.0, 1.0);

#include "layers_constants.h"
#include "K_G_grid.h"
#include "mesh.h"
#include "triangle_struct.h"
#include "interp_K_G_grid.h"
#include "vector_functions.h"
#include "GK_triangle.h"

void E_H_J_dip (TinyVector<complex<double>, 3>& E_J, TinyVector<complex<double>, 3>& H_J, const G_EJ_grid & G_EJ_tab, const G_HJ_grid & G_HJ_tab, const TinyVector<complex<double>,3>& J_src, const TinyVector<double,3>& r_src, const TinyVector<double,3>& r_obs, const layers_constants & LC) {

  int j;
  Array<complex<double>, 2> G_EJ (3, 3), G_HJ (3, 3);
  G_EJ = 0.0;
  G_HJ = 0.0;
  // computation of E_J due to a dipole located at r_src
  interp_G_EJ_mm (G_EJ, G_EJ_tab, r_obs, r_src);
  for (j=0 ; j<3 ; j++) E_J (j) = (G_EJ (j, 0) * J_src (0) + G_EJ (j, 1) * J_src (1) + G_EJ (j, 2) * J_src (2));

  // computation of ITo_H_inc due to a dipole located at r_ant
  interp_G_HJ_mm (G_HJ, G_HJ_tab, r_obs, r_src);
  for (j=0 ; j<3 ; j++) H_J (j) = (G_HJ (j, 0) * J_src (0) + G_HJ (j, 1) * J_src (1) + G_HJ (j, 2) * J_src (2));

}

void E_H_J_M_cosine (TinyVector<complex<double>, 3>& E_J_M, TinyVector<complex<double>, 3>& H_J_M, const mesh & MESH_ANTENNA, const TinyVector<double,3>& r_ant, const TinyVector<double,3>& x_hat_ant, const TinyVector<double,3>& y_hat_ant, const TinyVector<double,3>& r_obs, const G_EJ_grid & G_EJ_tab_ant_ant, const G_HJ_grid & G_HJ_tab_ant_ant, const layers_constants & LC, const G_EJ_grid & G_EJ_tab_ant_ant_dual, const G_HJ_grid & G_HJ_tab_ant_ant_dual, const layers_constants & LC_dual) {

  int i, j, N_TR_ANT = MESH_ANTENNA.triangles.rows();
  double sum_weigths;
  const double *xi, *eta, *weigths;
  int N_points = 9;
  IT_points (xi, eta, weigths, sum_weigths, N_points);

  double rho_1 = 0.29015095530181, rho_2 = 0.25187444860284, a_1 = 0.238, b_1 = 0.138; // BBHA antenna geometrical parameters
  double E_0 = sqrt(4.0*sqrt(LC.mu_0/LC.eps_0) / (a_1*b_1)); // so that P_rad = 1 (Balanis eq. 12-51)

  E_J_M = 0.0;
  H_J_M = 0.0;
  TinyVector<complex<double>, 3> E_J_M_triangle, H_J_M_triangle, E_J_j, H_J_j, E_M_j, H_M_j;
  for (i=0 ; i<N_TR_ANT ; i++) { // loop on antenna triangles
    complex<double> term_1;
    TinyVector<double, 3> r_j, r_J_on_antenna;
    TinyVector<complex<double>, 3> J_ant, M_ant;
    E_J_M_triangle = 0.0;
    H_J_M_triangle = 0.0;
    triangle T_antenna_i;
    construct_triangle (T_antenna_i, MESH_ANTENNA, i);
    for (j=0 ; j<N_points ; j++) { // loop on points of antenna triangle
      r_J_on_antenna = T_antenna_i.r_nodes (0) * xi[j] + T_antenna_i.r_nodes (1) * eta[j] + T_antenna_i.r_nodes (2) * (1-xi[j]-eta[j]);
      r_j = r_ant + r_J_on_antenna (0) * x_hat_ant + r_J_on_antenna (1) * y_hat_ant;
      term_1 = E_0 * cos(M_PI/a_1*r_J_on_antenna (0)) * exp(-I*LC.k_0*(pow2(r_J_on_antenna (0))/rho_2 + pow2(r_J_on_antenna (1))/rho_1)/2.0);
      J_ant = -1.0/sqrt(LC.mu_0/LC.eps_0) * term_1 * y_hat_ant;
      M_ant = term_1 * x_hat_ant;
      E_H_J_dip (E_J_j, H_J_j, G_EJ_tab_ant_ant, G_HJ_tab_ant_ant, J_ant, r_j, r_obs, LC);
      E_H_J_dip (H_M_j, E_M_j, G_EJ_tab_ant_ant_dual, G_HJ_tab_ant_ant_dual, M_ant, r_j, r_obs, LC_dual);
     
      E_J_M_triangle += (E_J_j - E_M_j) * weigths[j];
      H_J_M_triangle += (H_J_j + H_M_j) * weigths[j];
      
    }
    double norm_factor = T_antenna_i.A/sum_weigths;
    E_J_M += E_J_M_triangle * norm_factor;
    H_J_M += H_J_M_triangle * norm_factor;
  }

}

complex<double> R_horn_soil (const mesh & MESH_ANTENNA, const TinyVector<double,3>& r_ant, const TinyVector<double,3>& x_hat_ant, const TinyVector<double,3>& y_hat_ant, const G_EJ_grid & G_EJ_tab_ant_ant, const G_HJ_grid & G_HJ_tab_ant_ant, const layers_constants & LC, const G_EJ_grid & G_EJ_tab_ant_ant_dual, const G_HJ_grid & G_HJ_tab_ant_ant_dual, const layers_constants & LC_dual) {

  int i, j, N_TR_ANT = MESH_ANTENNA.triangles.rows();
  double sum_weigths;
  const double *xi, *eta, *weigths;
  int N_points = 9;
  IT_points (xi, eta, weigths, sum_weigths, N_points);

  double rho_1 = 0.29015095530181, rho_2 = 0.25187444860284, a_1 = 0.238, b_1 = 0.138; // BBHA antenna geometrical parameters
  double E_0 = sqrt(4.0*sqrt(LC.mu_0/LC.eps_0) / (a_1*b_1)); // so that P_rad = 1 (Balanis eq. 12-51)

  TinyVector<complex<double>, 3> E_2_j, H_2_j;
  complex<double> E_2_J_ant_triangle, H_2_M_ant_triangle;
  complex<double> result = 0.0;
  for (i=0 ; i<N_TR_ANT ; i++) { // loop on antenna triangles
    cout << "\r" << i*100/N_TR_ANT;
    flush(cout);
    complex<double> term_1;
    TinyVector<double, 3> r_j, r_J_on_antenna;
    TinyVector<complex<double>, 3> J_ant, M_ant;
    E_2_J_ant_triangle = 0.0;
    H_2_M_ant_triangle = 0.0;
    triangle T_antenna_i;
    construct_triangle (T_antenna_i, MESH_ANTENNA, i);
    for (j=0 ; j<N_points ; j++) { // loop on points of antenna triangle
      r_J_on_antenna = T_antenna_i.r_nodes (0) * xi[j] + T_antenna_i.r_nodes (1) * eta[j] + T_antenna_i.r_nodes (2) * (1-xi[j]-eta[j]);
      r_j = r_ant + r_J_on_antenna (0) * x_hat_ant + r_J_on_antenna (1) * y_hat_ant;
      term_1 = E_0 * cos(M_PI/a_1*r_J_on_antenna (0)) * exp(-I*LC.k_0*(pow2(r_J_on_antenna (0))/rho_2 + pow2(r_J_on_antenna (1))/rho_1)/2.0);
      J_ant = -1.0/sqrt(LC.mu_0/LC.eps_0) * term_1 * y_hat_ant;
      M_ant = term_1 * x_hat_ant;
      E_H_J_M_cosine (E_2_j, H_2_j, MESH_ANTENNA, r_ant, x_hat_ant, y_hat_ant, r_j, G_EJ_tab_ant_ant, G_HJ_tab_ant_ant, LC, G_EJ_tab_ant_ant_dual, G_HJ_tab_ant_ant_dual, LC_dual);
      E_2_J_ant_triangle += dot (E_2_j, J_ant) * weigths[j];
      H_2_M_ant_triangle += dot (H_2_j, M_ant) * weigths[j];
    }
    double norm_factor = T_antenna_i.A/sum_weigths;
    result += -0.5 * (E_2_J_ant_triangle - H_2_M_ant_triangle) * norm_factor; // see AP_1
  }
  return result;
}
