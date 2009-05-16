#ifndef ALPHA_COMPUTATION_H
#define ALPHA_COMPUTATION_H

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

std::complex<double> alpha_computation (const double & theta, /**< INPUT: angle \f$ \theta \f$ */
                                        const double & phi, /**< INPUT: angle \f$ \phi \f$ */
                                        const blitz::Array<std::complex<double>, 1>& h2_sph, /**< INPUT: spherical Hankel function */
                                        const blitz::TinyVector<double, 3>& r_mn, /**< INPUT: \f$ \mathbf{r}_{mn} = \mathbf{r}_m - \mathbf{r}_n \f$ */
                                        const int L, /**< the expansion number */
                                        const int L_prime, /**< the second expansion number */
                                        const std::complex<double>& k) /**< INPUT: the wavenumber */
/**
 * The \f$ \alpha \f$ translation coefficients are computed as follows:
 * 
 * \f[ 
 * \alpha^{L, L'}_{\mathbf{r}_{mn}} \left( \mathbf{\hat{s}} \right) = \sum \limits_{0 \leq l \leq L} \left(-j\right)^l \left( 2l + 1 \right) h_l^{\left(2\right)} \! \left( k r_{mn} \right) P_l\! \left( \mathbf{\hat{s}} \cdot \mathbf{\hat{r}}_{mn} \right) +  \sum \limits_{L+1 \leq l \leq L'} C_l^{L, L'} \left(-j\right)^L \left( 2L + 1 \right) h_L^{\left(2\right)} \! \left( k r_{mn} \right) P_l\! \left( \mathbf{\hat{s}} \cdot \mathbf{\hat{r}}_{mn} \right)
 * \f]
 *
 * where \f$ h_l^{\left(2\right)} \f$ is the spherical Hankel function of the second kind, \f$ P_l \f$ the Legendre polynomial of order \f$ l \f$, \f$ \mathbf{r}_{mn}  = \mathbf{r}_m - \mathbf{r}_n \f$ is the vector joining center of cube \f$ n \f$ to the center of cube \f$ m \f$, \f$ k \f$ is the wavenumber, \f$ \mathbf{\hat{s}} = \left[ \sin \theta \cos \phi, \sin \theta \sin \phi, \cos \theta \right]\f$  and \f[ C_l^{L, L'} = \cos \left( \frac{l - L}{L'-L} \, \frac{\pi}{2}   \right) \f] is a smoothing coefficient.
 *
 *
 */
{
  const double norm_r_mn = sqrt(dot(r_mn, r_mn)), sin_theta = sin(theta), cos_theta = cos(theta);
  blitz::TinyVector<double, 3> r_mn_hat, s_hat; 
  blitz::Array<double, 1> P_Leg(L_prime+1);
  blitz::Array<std::complex<double>, 1> coeff(L_prime+1);
  r_mn_hat = r_mn/norm_r_mn;
  s_hat = sin_theta*cos(phi), sin_theta*sin(phi), cos_theta;
  P_Legendre (P_Leg, dot(s_hat, r_mn_hat));
  for (int j=0 ; j<L+1 ; ++j) coeff(j) = pow(-I, j) * (2*j + 1.0) * h2_sph(j);
  const std::complex<double> const_coeff_lissage = pow(-I, L) * (2*L + 1.0) * h2_sph(L);
  for (int j=L+1 ; j<L_prime+1 ; ++j) coeff(j) = const_coeff_lissage * pow2(cos((j-L) * M_PI/2.0 / (L_prime-L)));
  return (-I * k / (16.0*M_PI*M_PI) * sum(coeff * P_Leg));
}

template <typename T>
void IT_theta_IT_phi_alpha_C (blitz::Array<std::complex<T>, 2>& alpha, /**< OUTPUT: 2D array \f$ \alpha \f$ */
                              const blitz::TinyVector<double, 3>& r_mn, /**< INPUT: \f$ r_{mn} = r_m - r_n \f$ */
                              const std::complex<double>& k, /**< INPUT: the wavenumber */
                              const int L, /**< the expansion number */
                              const int L_prime, /**< the second expansion number */
                              const blitz::Array<T, 1>& Xtheta, /**< 1D array of \f$ \theta \f$ angles */
                              const blitz::Array<T, 1>& Xphi) /**< 1D array of \f$ \phi \f$ angles */
{
  int kode = 1, M = 2, N = 1, nz, ierr;
  const double norm_r_mn = sqrt(dot(r_mn, r_mn));
  blitz::Array<std::complex<double>, 1> h2_sph(L+1);
  std::complex<double> z(k * norm_r_mn);
  for (int i=0 ; i<L+1 ; i++) {
    double fnu = i + 0.5;
    zbesh(z, fnu, kode, M, N, h2_sph(i), nz, ierr);
  }
  h2_sph *= sqrt(M_PI/(2.0 * z));
  for (int index_theta=0 ; index_theta<Xtheta.size() ; index_theta++) {
    for (int index_phi=0 ; index_phi<Xphi.size() ; index_phi++) {
      double theta = static_cast<double>(Xtheta(index_theta)), phi = static_cast<double>(Xphi(index_phi));
      alpha(index_theta, index_phi) = static_cast< complex<T> > ( alpha_computation (theta, phi, h2_sph, r_mn, L, L_prime, k) );
    }
  }
}

template <typename T>
void IT_theta_IT_phi_alpha_C2 (blitz::Array<std::complex<T>, 1>& alpha, /**< OUTPUT: 2D array \f$ \alpha \f$ */
                              const blitz::TinyVector<double, 3>& r_mn, /**< INPUT: \f$ r_{mn} = r_m - r_n \f$ */
                              const std::complex<double>& k, /**< INPUT: the wavenumber */
                              const int L, /**< the expansion number */
                              const int L_prime, /**< the second expansion number */
                              const blitz::Array<T, 2>& thetasPhis)
{
  int kode = 1, M = 2, N = 1, nz, ierr;
  const double norm_r_mn = sqrt(dot(r_mn, r_mn));
  blitz::Array<std::complex<double>, 1> h2_sph(L+1);
  std::complex<double> z(k * norm_r_mn);
  for (int i=0 ; i<L+1 ; i++) {
    double fnu = i + 0.5;
    zbesh(z, fnu, kode, M, N, h2_sph(i), nz, ierr);
  }
  h2_sph *= sqrt(M_PI/(2.0 * z));
  const int N_directions = alpha.size();
  for (int i=0 ; i<N_directions ; ++i) {
    double theta = static_cast<double>(thetasPhis(i, 0)), phi = static_cast<double>(thetasPhis(i, 1));
    alpha(i) = static_cast< std::complex<T> > ( alpha_computation (theta, phi, h2_sph, r_mn, L, L_prime, k) );
  }
}

template <typename T>
void IT_theta_IT_phi_alpha (blitz::Array<std::complex<T>, 2>& alpha, 
                            const blitz::Array<std::complex<double>, 1>& h2_sph, 
                            const blitz::TinyVector<double, 3>& r_mn, 
                            const std::complex<double>& k, 
                            const blitz::Array<double, 1>& Xtheta,
                            const blitz::Array<double, 1>& Xphi)
{
  const double norm_r_mn = sqrt(dot(r_mn, r_mn));
  const int L = h2_sph.size()-1, L_prime = h2_sph.size()-1;
  for (int index_theta=0 ; index_theta<Xtheta.size() ; index_theta++) {
    for (int index_phi=0 ; index_phi<Xphi.size() ; index_phi++) {
      alpha(index_theta, index_phi) = static_cast< std::complex<T> > ( alpha_computation (Xtheta(index_theta), Xphi(index_phi), h2_sph, r_mn, L, L_prime, k) );
    }
  }
}

#endif
