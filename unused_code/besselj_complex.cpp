#include <iostream>
#include <complex>
#include <blitz/array.h>

using namespace blitz;

const complex<double> I (0.0, 1.0);

/* this method comes from the paper of : Cornelis F. du Toit
 * "The numerical computations of Bessel functions of the
 * first and second kind for integer orders and complex arguments"
 * IEEE Transactions on Antennas and Propagation, Vol. 38, No. 9,
 * september 1990, pp. 1341-1349
 */
 
/***************************/
/* transitional function f */
/***************************/
double transition_function (complex<double> z, double n) {
  double f;
  double abs_z = abs(z);
  double cond = n-abs_z;
  if (cond <= -5) f = 1;
  else if (cond >= 5) f = 0;
  else f = pow2 (cos((n-abs_z+5)*M_PI/20));
  return f;
 }

/*********************************************/
/* Envelope of the asymptotic form of J_n(z) */
/*********************************************/
double envelope_bessel (complex<double> z, double n) {
  complex<double> beta = log((n/z) + I*sqrt(1.0-pow2 (n/z)))/I;
  complex<double> theta = z * sin(beta);
  double sigma = imag(theta)>=0?1:-1; 
  double epsilon = exp(sigma * imag(theta - beta * n));
  double f = (epsilon + transition_function(z,n)/epsilon)/(sqrt(2*M_PI*abs(z)) * pow (pow8 (abs(sin(beta))) + 1e-8, 1.0/16.0));
  return f;
}

int q_calculation (complex<double> z, double n, float err) {
  double a = ceil(1.1*abs(z)), q_low, q_high, c;
  double err_q = err *  envelope_bessel(z, 0)/2;
  double test_1 = envelope_bessel(z, a) - err_q;
  if (test_1 > 0) {
    q_high = ceil(a*1.414);
    while ((envelope_bessel(z, q_high) - err_q) >= 0) q_high = ceil(q_high*1.414);
    q_low = a;
  }
  else {
    q_low = floor(a/1.414);
    while ((envelope_bessel(z, q_low) - err_q) <= 0) q_low = floor(q_low/1.414);
    q_high = a;
  }
  c = ceil((q_low+q_high)/2);
  while ((q_high - q_low) > 3) {
    if ((envelope_bessel(z, c) - err_q) < 0) q_high = c;
    else q_low = c;
  c = ceil((q_low+q_high)/2);
  }
  int q = (int) c;
  return q;
}


/************************************************/
/* P_n(z) used in the asymptotic form of J_n(z) */
/************************************************/
complex<double> P_n (complex<double> z, double n, float err) {
  complex<double> eight_z_inv = 1.0/(z*8.0);
  double four_n_square = 4*n*n, PROD = 1.0, FACT = 1.0;
  complex<double> SUM_TERM, SUM (1.0, 0.0);
  double k = 1;
  do {
    PROD = PROD * (four_n_square - pow2 (2*(2*k-1) - 1)) * (four_n_square - pow2 (4*k - 1));
    FACT = FACT * 2*k * (2*k-1);
    SUM_TERM = pow(eight_z_inv, 2*k) * pow(-1, k)*PROD/FACT;
    SUM = SUM + SUM_TERM;
    k++;
  } while (abs(SUM_TERM) > err*abs(SUM));
  complex<double> P_n = SUM;
  return P_n;
  }


/************************************************/
/* Q_n(z) used in the asymptotic form of J_n(z) */
/************************************************/
complex<double> Q_n (complex<double> z, double n, float err) {
  complex<double> eight_z_inv = 1.0/(z*8.0);
  double four_n_square = 4*n*n, PROD = (four_n_square - 1), FACT = 1.0;
  complex<double> SUM_TERM, SUM = eight_z_inv * PROD;
  double k = 1;
  do {
    PROD = PROD * (four_n_square - pow2 (4*k - 1)) * (four_n_square - pow2 (2*(2*k+1) - 1));
    FACT = FACT * 2*k * (2*k+1);
    SUM_TERM = pow(eight_z_inv, 2*k+1) *  pow(-1, k) * PROD/FACT;
    SUM = SUM + SUM_TERM;
    k++;
  } while (abs(SUM_TERM) > err*abs(SUM));
  complex<double> Q_n = SUM;
  return Q_n;
  }

complex<double> besselj_complex (complex<double> z, double n, float err) {
  int i, k;
  complex<double> J_n;

  if (abs(z) > 0) {
    if (abs(z)<40) { // we use the recurrence relation to compute J(z)
      int q = q_calculation (z, n, err);
      Array<complex<double>, 1> B(q+2);
      B(q+1) = 0.0;
      B(q) = 1.0;

      for (i=q-1 ; i>-1 ; i=i-1) B(i) = B(i+1)*(2.0*(i+1))/z - B(i+2);

      /***********************************************/
      /* calculation of the normalization constant S */
      /***********************************************/
      complex<double> S, SUM (0.0, 0.0);
      if (imag(z)<=1) {
        for (k=1 ; k<=(q/2) ; k++) SUM = SUM + B(2*k);
        S = B(0) + SUM * 2.0;
      }
      else {
        for (k=1 ; k<=(q/2) ; k++) SUM = SUM + B(2*k) * pow(-1.0, k);
        S = (B(0) + SUM * 2.0)/cos(z);
      }
      J_n = B((int) n)/S;
    }
    else { /* if abs(z) >= 50 */
      J_n = sqrt(M_2_PI/z) * (P_n(z, n, err)*cos(z - n*M_PI_2 - M_PI_4) - Q_n(z, n, err) * sin(z - n*M_PI_2 - M_PI_4));
    }
  }
  else { /* if abs(z) = 0 */
    if (n == 0) J_n = 1.0;
    else J_n = 0.0;
  }
  return J_n;
}
