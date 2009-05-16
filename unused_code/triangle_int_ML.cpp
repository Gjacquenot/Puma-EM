#include <iostream>
#include <complex>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

using namespace blitz;

const complex<double> I (0.0, 1.0);

#include "K_G_grid.h"
#include "triangle_struct.h"
#include "interp_K_G_grid.h"
#include "GK_triangle.h"

void ITs_KA_ML (Array<complex<double>, 2>& ITs_KA, TinyVector<complex<double>,3>& ITs_KA_dot_rprime, const KA_grid & KA_tab, const TinyVector<double,3>& r, const triangle & Ts, const int N_points) {

  int j, line;
  double sum_weigths, norm_factor;
  TinyVector<double, 3> rprime;
  Array<complex<double>, 2> KA (3, 3);

  const double *xi, *eta, *weigths;
  IT_points (xi, eta, weigths, sum_weigths, N_points);

  ITs_KA_dot_rprime = 0.0; // TinyVector<complex<double>,3>
  ITs_KA = 0.0; // Array<complex<double>, 2>
  for (j=0 ; j<N_points ; j++) {
    rprime = Ts.r_nodes (0)*xi[j] + Ts.r_nodes (1)*eta[j] + Ts.r_nodes (2)*(1-xi[j]-eta[j]);
    interp_KA_mm (KA, KA_tab, r, rprime);
    KA *= weigths[j];
    ITs_KA += KA;
    for (line=0 ; line<3 ; line++) ITs_KA_dot_rprime (line) +=  KA (line,0) * rprime (0) + KA (line,1) * rprime (1) + KA (line,2) * rprime (2);
  }
  norm_factor = Ts.A/sum_weigths;
  ITs_KA *= norm_factor;
  ITs_KA_dot_rprime *= norm_factor;
}

void ITo_ITs_KA_ML (Array<complex<double>, 2>& ITo_ITs_KA, TinyVector<complex<double>, 3>& ITo_r_dot_ITs_KA, TinyVector<complex<double>, 3>& ITo_ITs_KA_dot_rprime, complex<double>& ITo_r_dot_ITs_KA_dot_rprime, TinyVector<complex<double>, 3>& ITo_n_hat_X_r_dot_ITs_KA, complex<double>& ITo_n_hat_X_r_dot_ITs_KA_dot_rprime, const KA_grid & KA_tab, const triangle & To, const triangle & Ts, const complex<double> k, const int N_points_o, const int N_points_s) {

  int j, ind;
  double sum_weigths, norm_factor;
  TinyVector<double, 3> r, n_hat_X_r;
  TinyVector<complex<double>, 3> ITs_KA_dot_rprime_j, r_dot_ITs_KA_j, n_hat_X_r_dot_ITs_KA_j;
  Array<complex<double>, 2> ITs_KA_j (3, 3);

  const double *xi, *eta, *weigths;
  IT_points (xi, eta, weigths, sum_weigths, N_points_o);

  ITo_ITs_KA = 0.0; // Array<complex<double>, 2>
  ITo_r_dot_ITs_KA = 0.0; // TinyVector<complex<double>, 3>
  ITo_ITs_KA_dot_rprime = 0.0; // TinyVector<complex<double>, 3>
  ITo_r_dot_ITs_KA_dot_rprime = 0.0; // complex<double>
  ITo_n_hat_X_r_dot_ITs_KA = 0.0; // TinyVector<complex<double>, 3>
  ITo_n_hat_X_r_dot_ITs_KA_dot_rprime = 0.0; // complex<double>
  for (j=0 ; j<N_points_o ; j++) {
    r = To.r_nodes (0)*xi[j] + To.r_nodes (1)*eta[j] + To.r_nodes (2)*(1-xi[j]-eta[j]);
    n_hat_X_r = cross(To.n_hat, r);
    ITs_KA_ML (ITs_KA_j, ITs_KA_dot_rprime_j, KA_tab, r, Ts, N_points_s);

    ITs_KA_j *= weigths[j];
    ITo_ITs_KA += ITs_KA_j;

    for (ind=0 ; ind<3 ; ind++) r_dot_ITs_KA_j (ind) = r(0) * ITs_KA_j(0,ind) + r(1) * ITs_KA_j(1,ind) + r(2) * ITs_KA_j(2,ind);
    ITo_r_dot_ITs_KA += r_dot_ITs_KA_j;

    ITs_KA_dot_rprime_j *= weigths[j];
    ITo_ITs_KA_dot_rprime += ITs_KA_dot_rprime_j;

    ITo_r_dot_ITs_KA_dot_rprime += dot (r, ITs_KA_dot_rprime_j);

    for (ind=0 ; ind<3 ; ind++) n_hat_X_r_dot_ITs_KA_j (ind) = n_hat_X_r(0) * ITs_KA_j(0,ind) + n_hat_X_r(1) * ITs_KA_j(1,ind) + n_hat_X_r(2) * ITs_KA_j(2,ind);
    ITo_n_hat_X_r_dot_ITs_KA += n_hat_X_r_dot_ITs_KA_j;

    ITo_n_hat_X_r_dot_ITs_KA_dot_rprime += dot (n_hat_X_r, ITs_KA_dot_rprime_j);
  }
  norm_factor = To.A/sum_weigths;
  ITo_ITs_KA *= norm_factor;
  ITo_r_dot_ITs_KA *= norm_factor;
  ITo_ITs_KA_dot_rprime *= norm_factor;
  ITo_r_dot_ITs_KA_dot_rprime *= norm_factor;
  ITo_n_hat_X_r_dot_ITs_KA *= norm_factor;
  ITo_n_hat_X_r_dot_ITs_KA_dot_rprime *= norm_factor;
}

void ITs_Kphi_ML (complex<double>& ITs_K_phi, TinyVector<complex<double>, 3>& ITs_grad_K_phi, const Kphi_grid & Kphi_tab, const TinyVector<double,3>& r, const triangle & Ts, const int N_points_s) {

  int j, line;
  double sum_weigths, norm_factor;
  TinyVector<double, 3> rprime;
  complex<double> K_phi;
  TinyVector<complex<double>, 3> grad_K_phi;

  const double *xi, *eta, *weigths;
  IT_points (xi, eta, weigths, sum_weigths, N_points_s);

  ITs_K_phi = 0.0; // complex<double>
  ITs_grad_K_phi = 0.0; // TinyVector<complex<double>, 3>
  for (j=0 ; j<N_points_s ; j++) {
    rprime = Ts.r_nodes (0)*xi[j] + Ts.r_nodes (1)*eta[j] + Ts.r_nodes (2)*(1-xi[j]-eta[j]);
    interp_Kphi_mm (K_phi, Kphi_tab, r, rprime);
    interp_grad_K_phi (grad_K_phi, Kphi_tab, r, rprime);
    ITs_K_phi += K_phi * weigths[j];
    ITs_grad_K_phi += grad_K_phi * weigths[j];
  }
  norm_factor = Ts.A/sum_weigths;
  ITs_K_phi *= norm_factor;
  ITs_grad_K_phi *= norm_factor;
}

void ITo_ITs_Kphi_ML (complex<double>& ITo_ITs_K_phi, TinyVector<complex<double>, 3>& ITo_ITs_grad_K_phi, complex<double>& ITo_n_hat_X_r_dot_ITs_grad_K_phi, const Kphi_grid & Kphi_tab, const triangle & To, const triangle & Ts, const complex<double> k, const int N_points_o, const int N_points_s) {

  int j, ind;
  double sum_weigths, norm_factor;
  TinyVector<double, 3> r, n_hat_X_r;
  complex<double> ITs_K_phi_j;
  TinyVector<complex<double>, 3> ITs_grad_K_phi_j;

  const double *xi, *eta, *weigths;
  IT_points (xi, eta, weigths, sum_weigths, N_points_o);

  ITo_ITs_K_phi = 0.0; // complex<double>
  ITo_ITs_grad_K_phi = 0.0; // TinyVector<complex<double>, 3>
  ITo_n_hat_X_r_dot_ITs_grad_K_phi = 0.0; // complex<double>
  for (j=0 ; j<N_points_o ; j++) {
    r = To.r_nodes (0)*xi[j] + To.r_nodes (1)*eta[j] + To.r_nodes (2)*(1-xi[j]-eta[j]);
    n_hat_X_r = cross(To.n_hat, r);
    ITs_Kphi_ML (ITs_K_phi_j, ITs_grad_K_phi_j, Kphi_tab, r, Ts, N_points_s);
    ITo_ITs_K_phi += ITs_K_phi_j * weigths[j];

    ITs_grad_K_phi_j *= weigths[j];
    ITo_ITs_grad_K_phi += ITs_grad_K_phi_j;
    ITo_n_hat_X_r_dot_ITs_grad_K_phi += dot (n_hat_X_r, ITs_grad_K_phi_j);
  }
  norm_factor = To.A/sum_weigths;
  ITo_ITs_K_phi *= norm_factor;
  ITo_ITs_grad_K_phi *= norm_factor;
  ITo_n_hat_X_r_dot_ITs_grad_K_phi *= norm_factor;
}

void ITs_G_HJ_ML (Array<complex<double>, 2>& ITs_G_HJ, TinyVector<complex<double>,3>& ITs_G_HJ_dot_rprime, const G_HJ_grid & G_HJ_tab, const TinyVector<double,3>& r, const triangle & Ts, const int N_points) {

  int j, line;
  double sum_weigths, norm_factor;
  TinyVector<double, 3> rprime;
  Array<complex<double>, 2> G_HJ (3, 3);

  const double *xi, *eta, *weigths;
  IT_points (xi, eta, weigths, sum_weigths, N_points);

  ITs_G_HJ_dot_rprime = 0.0; // TinyVector<complex<double>,3>
  ITs_G_HJ = 0.0; // Array<complex<double>, 2>
  for (j=0 ; j<N_points ; j++) {
    rprime = Ts.r_nodes (0)*xi[j] + Ts.r_nodes (1)*eta[j] + Ts.r_nodes (2)*(1-xi[j]-eta[j]);
    interp_G_HJ_mm (G_HJ, G_HJ_tab, r, rprime);
    G_HJ *= weigths[j];
    ITs_G_HJ += G_HJ;
    for (line=0 ; line<3 ; line++) ITs_G_HJ_dot_rprime (line) +=  G_HJ (line,0) * rprime (0) + G_HJ (line,1) * rprime (1) + G_HJ (line,2) * rprime (2);
  }
  norm_factor = Ts.A/sum_weigths;
  ITs_G_HJ *= norm_factor;
  ITs_G_HJ_dot_rprime *= norm_factor;
}

void ITo_ITs_G_HJ_ML (Array<complex<double>, 2>& ITo_ITs_G_HJ, TinyVector<complex<double>, 3>& ITo_r_dot_ITs_G_HJ, TinyVector<complex<double>, 3>& ITo_ITs_G_HJ_dot_rprime, complex<double>& ITo_r_dot_ITs_G_HJ_dot_rprime, TinyVector<complex<double>, 3>& ITo_n_hat_X_r_dot_ITs_G_HJ, complex<double>& ITo_n_hat_X_r_dot_ITs_G_HJ_dot_rprime, const G_HJ_grid & G_HJ_tab, const triangle & To, const triangle & Ts, const complex<double> k, const int N_points_o, const int N_points_s) {

  int j, ind;
  double sum_weigths, norm_factor;
  TinyVector<double, 3> r, n_hat_X_r;
  complex<double> r_dot_ITs_G_HJ_dot_rprime_j, n_hat_X_r_dot_ITs_G_HJ_dot_rprime_j;
  TinyVector<complex<double>, 3> ITs_G_HJ_dot_rprime_j, r_dot_ITs_G_HJ_j, n_hat_X_r_dot_ITs_G_HJ_j;
  Array<complex<double>, 2> ITs_G_HJ_j (3, 3);

  const double *xi, *eta, *weigths;
  IT_points (xi, eta, weigths, sum_weigths, N_points_o);

  ITo_ITs_G_HJ = 0.0; // Array<complex<double>, 2>
  ITo_r_dot_ITs_G_HJ = 0.0; // TinyVector<complex<double>, 3>
  ITo_ITs_G_HJ_dot_rprime = 0.0; // TinyVector<complex<double>, 3>
  ITo_r_dot_ITs_G_HJ_dot_rprime = 0.0; // complex<double>
  ITo_n_hat_X_r_dot_ITs_G_HJ = 0.0; // TinyVector<complex<double>, 3>
  ITo_n_hat_X_r_dot_ITs_G_HJ_dot_rprime = 0.0; // complex<double>
    for (j=0 ; j<N_points_o ; j++) {
      r = To.r_nodes (0)*xi[j] + To.r_nodes (1)*eta[j] + To.r_nodes (2)*(1-xi[j]-eta[j]);
      n_hat_X_r = cross (To.n_hat, r);
      ITs_G_HJ_ML (ITs_G_HJ_j, ITs_G_HJ_dot_rprime_j, G_HJ_tab, r, Ts, N_points_s);

      ITs_G_HJ_j *= weigths[j];
      ITo_ITs_G_HJ += ITs_G_HJ_j;

      for (ind=0 ; ind<3 ; ind++) r_dot_ITs_G_HJ_j (ind) = r(0) * ITs_G_HJ_j(0,ind) + r(1) * ITs_G_HJ_j(1,ind) + r(2) * ITs_G_HJ_j(2,ind);
      ITo_r_dot_ITs_G_HJ += r_dot_ITs_G_HJ_j;

      ITs_G_HJ_dot_rprime_j *= weigths[j];
      ITo_ITs_G_HJ_dot_rprime += ITs_G_HJ_dot_rprime_j;
     
      r_dot_ITs_G_HJ_dot_rprime_j = dot (r, ITs_G_HJ_dot_rprime_j);
      ITo_r_dot_ITs_G_HJ_dot_rprime += r_dot_ITs_G_HJ_dot_rprime_j;

      for (ind=0 ; ind<3 ; ind++) n_hat_X_r_dot_ITs_G_HJ_j (ind) = n_hat_X_r(0) * ITs_G_HJ_j(0,ind) + n_hat_X_r(1) * ITs_G_HJ_j(1,ind) + n_hat_X_r(2) * ITs_G_HJ_j(2,ind);
      ITo_n_hat_X_r_dot_ITs_G_HJ += n_hat_X_r_dot_ITs_G_HJ_j;

      n_hat_X_r_dot_ITs_G_HJ_dot_rprime_j = dot (n_hat_X_r, ITs_G_HJ_dot_rprime_j);
      ITo_n_hat_X_r_dot_ITs_G_HJ_dot_rprime += n_hat_X_r_dot_ITs_G_HJ_dot_rprime_j;
    }
  norm_factor = To.A/sum_weigths;
  ITo_ITs_G_HJ *= norm_factor;
  ITo_r_dot_ITs_G_HJ *= norm_factor;
  ITo_ITs_G_HJ_dot_rprime *= norm_factor;
  ITo_r_dot_ITs_G_HJ_dot_rprime *= norm_factor;
  ITo_n_hat_X_r_dot_ITs_G_HJ *= norm_factor;
  ITo_n_hat_X_r_dot_ITs_G_HJ_dot_rprime *= norm_factor;
}
