#include <complex>
#include <blitz/array.h>

using namespace std;

#include "interpolation.h"

template <typename T>
int find_index(const T x, const blitz::Array<T, 1>& x_i) {
  // a log(n) routine for finding the index of a x in an array x_i of size n
  // the array is supposed to be monotonically increasing
  int ind_inf, ind_sup, ind_mid;
  const int n = x_i.size();
  if ( (x < x_i(0)) || (x > x_i(n-1)) ) {
    cout << "find_index(): x outside range. x_i = " << x_i << ", x = " << x << endl;
    exit(1);
  }
  else {
    ind_inf = 0;
    ind_sup = n-1;
    while(ind_sup-ind_inf > 1) {
      ind_mid = (ind_sup+ind_inf)/2;
      if (x <= x_i(ind_mid)) ind_sup = ind_mid;
      else ind_inf = ind_mid;
    }
    return ind_inf;
  }
}


/**************************************************************/
/*********** regular step Lagrange interpolation **************/
/**************************************************************/

void Lagrangian_regular_matrix(const int n, 
                               blitz::Array<double, 2>& M) 
{
  switch (n) {
    case 1: M.resize(n+1, n+1);
            M = 1.0, -1.0,
                0.0,  1.0;
             break;
    case 2: M.resize(n+1, n+1);
            M = 1.0, -1.5, 0.5,
                0.0,  2.0, -1.0,
                0.0, -0.5, 0.5;
            break;
    case 3: M.resize(n+1, n+1);
            M = 1.0, -11.0/6.0, 1.0, -1.0/6.0,
                0.0, 3.0, -2.5, 0.5,
                0.0, -1.5, 2.0, -0.5,
                0.0, 1.0/3.0, -0.5, 1.0/6.0;
            break;
    case 4: M.resize(n+1, n+1);
            M = 1.0, -25.0/12.0, 35.0/24.0, -5.0/12.0, 1.0/24.0,
                0.0, 4.0, -13.0/3.0, 1.5, -1.0/6.0,
                0.0, -3.0, 4.75, -2.0, 0.25,
                0.0, 4.0/3.0, -7.0/3.0, 7.0/6.0, -1.0/6.0,
                0.0, -0.25, 11.0/24.0, -0.25, 1.0/24.0;
            break;
    default: cout << "Lagrangian_regular_matrix(): order too high for Lagrange interpolation" << endl; exit(1);
  }
}

void Lagrange_vector_fixedstep_interpolation(blitz::Array<std::complex<double>, 1>& y, 
                                             const blitz::Array<double, 1>& x, 
                                             const blitz::Array<double, 1>& xi, 
                                             const blitz::Array<std::complex<double>, 1>& yi, 
                                             const blitz::Array<double, 2>& M) 
{
  const int Ni = xi.size(), N = x.size(), n = M.rows()-1; // n is the number of DeltaXi
  int startInd, i, j, k;
  double h = xi(1) - xi(0), s;
  if ((n > Ni-1) || (n < 1)) {
    cout << "Lagrange_vector_fixedstep_interpolation(): order n of interpolation too high or too low given the number of abscissas. n = " << n << ". Exiting..." << endl;
    exit(1);
  }
  else {
    blitz::Array<double, 1> S(n+1), M_prod_S(n+1);
    for (k=0 ; k<N ; k++) {
      i = (int) floor((x(k)-xi(0))/h);
      startInd = i-n/2;
      if (startInd<0) startInd = 0;
      else if (startInd+n>Ni-1) startInd = (Ni-1) - n;
      s = ( x(k)-xi(startInd) ) / h;
      for (j=0 ; j<n+1 ; j++) S(j) = pow(s, j);
      for (j=0 ; j<n+1 ; j++) M_prod_S(j) = sum(M(j, blitz::Range::all()) * S);
      y(k) = blitz::sum(yi(blitz::Range(startInd, startInd+n)) * M_prod_S);
    }
  }
}


void decimateAbscissa (blitz::Array<double, 1>& x, 
                       const blitz::Array<double, 1>& xi, 
                       const int decimFact)
{
  const int Nxi = xi.size(), Nx = decimFact*(Nxi-1) + 1;
  x.resize(Nx);
  const double Dx = (xi(Nxi-1)-xi(0))/(Nx-1);
  for (int j=0 ; j<Nx ; j++) x(j) = xi(0) + j*Dx;
  x(blitz::Range(0, blitz::toEnd, decimFact)) = xi; // we get rid off the roundoff errors at periodic points
}

void decimate_2D (blitz::Array<std::complex<double>, 2> Y, 
                  const blitz::Array<std::complex<double>, 2>& Y_i, 
                  const blitz::Array<double, 1>& X1_i, // the abscissas following 1st dimension
                  const blitz::Array<double, 1>& X2_i, // the abscissas following 2nd dimension
                  const int n) 
{
  unsigned int N_X1_i = Y_i.extent(0), N_X2_i = Y_i.extent(1);
  if ( (N_X1_i != X1_i.size()) || (N_X2_i != X2_i.size()) ) {
    cout << "decimate_2D() : (N_X1_i != X1_i.size()) || (N_X2_i != X2_i.size())" << endl;
    exit(1);
  }
  Y = 0.0;
  Y(blitz::Range(0, blitz::toEnd, 2), blitz::Range(0, blitz::toEnd, 2)) = Y_i;
  blitz::Array<complex<double>, 1> y_tmp; // a temporary array

  // we now define the theta and phi arrays (interpolation abscissas)
  blitz::Array<double, 1> X1, X2; 
  decimateAbscissa (X1, X1_i, 2);
  decimateAbscissa (X2, X2_i, 2);

  blitz::Array<double, 2> M;
  Lagrangian_regular_matrix(n, M);
  // we first interpolate following phi
  y_tmp.resize(N_X2_i-1);
  for (unsigned int j=0 ; j<2*N_X1_i - 1 ; j=j+2) {
    Lagrange_vector_fixedstep_interpolation(y_tmp, X2(blitz::Range(1, blitz::toEnd, 2)), X2(blitz::Range(0, blitz::toEnd, 2)), Y(j, blitz::Range(0, blitz::toEnd, 2)), M);
    Y(j, blitz::Range(1, blitz::toEnd, 2)) = y_tmp;
  }
  // we then interpolate following X1. Pay attention to the fact
  // that we now have (2*N_X2_i - 1) values following X1!!
  // Therefore we do not use Y_i for interpolation anymore.
  y_tmp.resize(N_X1_i-1);
  for (unsigned int j=0 ; j<2*N_X2_i-1 ; j++) {
    Lagrange_vector_fixedstep_interpolation(y_tmp, X1(blitz::Range(1, blitz::toEnd, 2)), X1(blitz::Range(0, blitz::toEnd, 2)), Y(blitz::Range(0, blitz::toEnd, 2), j), M);
    Y(blitz::Range(1, blitz::toEnd, 2), j) = y_tmp;
  }
}

/***************************************************************************/
/****************** Lagrange Interpolation Matrices ************************/
/***************************************************************************/

int findStartInd(const float x, const blitz::Array<float, 1>& xi, const int NOrder, const int CYCLIC)
{
  int Nxi = xi.size(), startInd, ind;
  if (CYCLIC<=0) {
    if (x<xi(0)) ind = 0;
    else if (x>xi(Nxi-1)) ind = Nxi-1;
    else ind = find_index(x, xi);
    startInd = ind - NOrder/2;
    if (startInd<0) startInd = 0;
    else if (startInd+NOrder>Nxi-1) startInd = (Nxi-1) - NOrder;
  }
  else {
    if (x<xi(0)) ind = -1;
    else if (x>xi(Nxi-1)) ind = Nxi-1;
    else ind = find_index(x, xi);
    startInd = ind - NOrder/2;
  }
  return startInd;
}

void xiTmpConstruction(blitz::Array<float, 1>& xiTmp,
                       const blitz::Array<float, 1>& xi,
                       const float a, // lim inf of the interval
                       const float b, // lim sup of the interval
                       const int startInd,
                       const int NOrder,
                       const int CYCLIC,
                       const int INCLUDED_BOUNDARIES)
{
  int i, ind, Nxi = xi.size();
  if ( (CYCLIC<=0) || ((startInd>-1) && (startInd+NOrder<Nxi)) ) {
    for (i=0 ; i<NOrder+1 ; i++) xiTmp(i) = xi(startInd + i);
  }
  else {
    if (INCLUDED_BOUNDARIES==0) { // (xi(0)!=a) || (xi(Nxi-1)!=b)
      for (i=0 ; i<NOrder+1 ; i++) {
        ind = startInd + i;
        if (ind<0) xiTmp(i) = xi(Nxi + ind) + a - b;
        else if (ind>Nxi-1) xiTmp(i) = xi(ind-Nxi) + b - a;
        else xiTmp(i) = xi(ind);
      }
    }
    else { // (xi(0)==a) && (xi(Nxi-1)==b)
      for (i=0 ; i<NOrder+1 ; i++) {
        ind = startInd + i;
        if (ind<0) xiTmp(i) = xi(Nxi + ind - 1) + a - b; // we "jump" over xi[Nxi-1] (because xi[Nxi-1]==b)
        else if (ind>Nxi-1) xiTmp(i) = xi(ind-Nxi+1) + b - a; // we "jump" over xi[0] (because xi[0]==a)
        else xiTmp(i) = xi(ind);
      }
    }
  }
}

void indexesConstruction(blitz::Array<int, 1> indexesTmp,
                         const int startInd,
                         const int NOrder)
// pretty simple right now but could become more complicated
{
  for (int i=0 ; i<NOrder+1 ; i++) indexesTmp(i) = startInd + i;
}

void index2DtoIndex1D(int & index1D,
                      float & sign,
                      const int lIndex,
                      const int cIndex,
                      const int Nl,
                      const int Nc,
                      const int BOUNDARIES_THETA)
/**
 * This function is very delicate, as it will perform the index transformation from 2-D 
 * into an index suitable for a 2-D flattened into 1-D array. Nl and Nc are the dimensions 
 * of the 2-D array. Nc MUST be an even number.
 *
 * 1. Interpolation following theta
 * ================================
 *
 * 1.1 interpolation using theta < 0
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Since we perform interpolation over a sphere, we must pay attention
 * to the fact that, for 0 <= theta <= pi:
 *
 *     F(-theta, phi) = -F(theta, phi+pi)      if phi <= pi
 *                    = -F(theta, phi-pi)      if phi > pi
 *
 * from (equation (31) from "Optimal Interpolation of Radiated Fields over a Sphere", IEEE AP november 1991)
 *
 * a requirement that, in terms of indexes, translates into:
 *
 * 1) BOUNDARIES_THETA == 0, i>=1
 *
 *    F(-i, j) = -F(i-1, j+Nc/2)                if j <= Nc/2
 *             = -F(i-1, j-Nc/2)                if j >  Nc/2
 *
 * 2) BOUNDARIES_THETA == 1, i>=1
 *
 *    F(-i, j) = -F(i, j+Nc/2)                  if j <= Nc/2
 *             = -F(i, j-Nc/2)                  if j >  Nc/2
 *
 * One will note the apparition of the minus (-) sign that occurs when theta (or index i) < 0.
 * That minus sign is supposed to multiply the function sample. Therefore this sign should be
 * included into the "coefficients" matrix.
 *
 * 1.2 interpolation using theta > pi
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * If theta > pi, we have that:
 *
 *    F(theta, phi) = F(theta-2*pi, phi) = -F(2*pi-theta, phi + pi)
 *
 * which, in terms of indexes, is translated into:
 *
 *    i > Nl-1: F(i, j) = F(i - 2*Nl + 2, j)      if BOUNDARIES_THETA == 1
 *                      = F(i - 2*Nl, j)          if BOUNDARIES_THETA == 0
 *
 * The resulting row indexes are negative and we therefore must go back to the case of negative 
 * row indexes explained in 1.1.
 *
 * 2. Interpolation following phi
 * ==============================
 * For interpolation following phi, the rules are somewhat similar, except that in this case
 * values of the function over the full interval [0..2*pi] are given. Note that we do not consider
 * the case BOUNDARIES_Phi == 1, since "TRAP" and "PONCELET" methods do not include them...
 *
 *    j < 0   : F(i, j) = F(i, j + Nc)
 *    j > Nc-1: F(i, j) = F(i, j - Nc)
 */
{ // coefficients for all the elements of vector x
  int lIndexTmp = lIndex, cIndexTmp = cIndex;
  sign = 1.0;
  if (lIndexTmp > Nl - 1) lIndexTmp = (BOUNDARIES_THETA == 1) ? lIndex - 2*Nl + 2 : lIndex - 2*Nl;
  if (lIndexTmp < 0) {
    lIndexTmp = (BOUNDARIES_THETA == 1) ? abs(lIndexTmp) : abs(lIndexTmp) - 1;
    sign = -1.0;
    cIndexTmp = (cIndex>Nc/2) ? cIndex-Nc/2 : cIndex+Nc/2;
  }
  if (cIndexTmp < 0) cIndexTmp += Nc;
  else if (cIndexTmp > Nc-1) cIndexTmp -= Nc;
  index1D = lIndexTmp + cIndexTmp * Nl;
}


void LagrangeInterpolationCoeffs(blitz::Array<float, 1> coeffs,
                                 const float x,
                                 const blitz::Array<float, 1>& x_i,
                                 const int PERIODIC)
{
  int i, j, N_points = x_i.size();
  //coeffs.resize(N_points);
  float B, C;
  for (i=0 ; i<N_points ; i++) {
    B = C = 1.0;
    for (j=0 ; j<N_points ; j++) {
      if (j!=i) {
        if (PERIODIC==1) {
          B *= sin( 0.5 * (x - x_i(j)) );
          C *= sin( 0.5 * (x_i(i) - x_i(j)) );
        }
        else {
          B *= (x - x_i(j));
          C *= (x_i(i) - x_i(j));
        }
      }
    }
    coeffs(i) = B/C;
  }
}

LagrangeFastInterpolator2D::LagrangeFastInterpolator2D (void) {}

LagrangeFastInterpolator2D::LagrangeFastInterpolator2D(const blitz::Array<float, 1>& x,
                                                       const blitz::Array<float, 1>& xi,
                                                       const float axi, // lim inf of xi interval
                                                       const float bxi, // lim sup of xi interval
                                                       const int INCLUDED_BOUNDARIES_xi,
                                                       const int NOrderXi,
                                                       const int PERIODIC_xi,
                                                       const int CYCLIC_xi,
                                                       const blitz::Array<float, 1>& y,
                                                       const blitz::Array<float, 1>& yi,
                                                       const float ayi, // lim inf of yi interval
                                                       const float byi, // lim sup of yi interval
                                                       const int INCLUDED_BOUNDARIES_yi,
                                                       const int NOrderYi,
                                                       const int PERIODIC_yi,
                                                       const int CYCLIC_yi)
/**
 * This interpolator works on 2-D Arrays which are reshaped as 1-D column vector.
 * The reason for this is that we will need the interpolator for anterpolation,
 * and this operation is much easier to perform on a 1-D column array.
 *
 * Basically, this interpolator is a 2-D array which will be matrixmultiplied by a 
 * 1-D column vector to produce an interpolated 1-D column vector.
 *
 * The 1-D column vector is constructed from the 2-D Array by piling its columns 
 * one under the other, i.e.:
 *
 *      V_2D(i,j) = V_1D(i + j*Nxi)
 */
{
  blitz::Range all = blitz::Range::all();
  const int Nx = x.size(), Nxi = xi.size(), Ny = y.size(), Nyi = yi.size();
  if ((NOrderXi > Nxi-1) || (NOrderXi < 1)) {
    cout << "LagrangeFastInterpolator2D::LagrangeFastInterpolator2D: NOrderXi of interpolation too high or too low given the number of abscissas xi. NOrderXi = " << NOrderXi << ". Exiting..." << endl;
    exit(1);
  }
  if ((NOrderYi > Nyi-1) || (NOrderYi < 1)) {
    cout << "LagrangeFastInterpolator2D::LagrangeFastInterpolator2D: NOrderYi of interpolation too high or too low given the number of abscissas yi. NOrderYi = " << NOrderYi << ". Exiting..." << endl;
    exit(1);
  }
  coefficientsForLinesInterp.resize(Nx * Ny, (NOrderYi+1));
  coefficientsForColumnsInterp.resize(Nx * Nyi, (NOrderXi+1));
  indexesForLinesInterp.resize(Nx * Ny, (NOrderYi+1));
  indexesForColumnsInterp.resize(Nx * Nyi, (NOrderXi+1));
  coefficientsForLinesInterp = 0.0;
  coefficientsForColumnsInterp = 0.0;
  indexesForLinesInterp = -1;
  indexesForColumnsInterp = -1;

  int startIndXi, startIndYi, index1D, i, j, p, q;
  blitz::Array<float, 1> xiTmp(NOrderXi+1), yiTmp(NOrderYi+1);
  blitz::Array<float, 2> coeffsXi(Nx, NOrderXi+1), coeffsYi(Ny, NOrderYi+1);
  blitz::Array<int, 2> indexesXi(Nx, NOrderXi+1), indexesYi(Ny, NOrderYi+1);
  // coefficients for all the elements of vector x
  for (i=0 ; i<Nx; i++) {
    startIndXi = findStartInd(x(i), xi, NOrderXi, CYCLIC_xi);
    xiTmpConstruction(xiTmp, xi, axi, bxi, startIndXi, NOrderXi, CYCLIC_xi, INCLUDED_BOUNDARIES_xi);
    LagrangeInterpolationCoeffs(coeffsXi(i, all), x(i), xiTmp, PERIODIC_xi);
    indexesConstruction(indexesXi(i, all), startIndXi, NOrderXi);
  }
  // coefficients for all the elements of vector y
  for (i=0 ; i<Ny; i++) {
    startIndYi = findStartInd(y(i), yi, NOrderYi, CYCLIC_yi);
    xiTmpConstruction(yiTmp, yi, ayi, byi, startIndYi, NOrderYi, CYCLIC_yi, INCLUDED_BOUNDARIES_yi);
    LagrangeInterpolationCoeffs(coeffsYi(i, all), y(i), yiTmp, PERIODIC_yi);
    indexesConstruction(indexesYi(i, all), startIndYi, NOrderYi);
  }
  // Lines and Columns interpolators construction
  float sign;
  // construction of Lines interpolator
  for (i=0 ; i<Nx; i++) {
    for (j=0 ; j<Ny; j++) {
      for (q=0 ; q<(NOrderYi+1) ; q++) {
        index2DtoIndex1D(index1D, sign, i, indexesYi(j, q), Nx, Nyi, INCLUDED_BOUNDARIES_yi);
        coefficientsForLinesInterp(i + j*Nx, q) = coeffsYi(j, q) * sign;
        indexesForLinesInterp(i + j*Nx, q) = index1D;
      }
    }
  }
  // construction of Columns interpolator
  for (i=0 ; i<Nx; i++) {
    for (j=0 ; j<Nyi; j++) {
      for (p=0 ; p<(NOrderXi+1) ; p++) {
        index2DtoIndex1D(index1D, sign, indexesXi(i, p), j, Nxi, Nyi, INCLUDED_BOUNDARIES_xi);
        coefficientsForColumnsInterp(i + j*Nx, p) = coeffsXi(i, p) * sign;
        indexesForColumnsInterp(i + j*Nx, p) = index1D;
      }
    }
  }
}

LagrangeFastInterpolator2D::LagrangeFastInterpolator2D (const LagrangeFastInterpolator2D & lfi) // copy
{
  coefficientsForLinesInterp.resize(lfi.getNCoefficientsForLinesInterp(), lfi.getNOrderCoefficientsForLinesInterp());
  coefficientsForColumnsInterp.resize(lfi.getNCoefficientsForColumnsInterp(), lfi.getNOrderCoefficientsForColumnsInterp());
  indexesForLinesInterp.resize(lfi.getNCoefficientsForLinesInterp(), lfi.getNOrderCoefficientsForLinesInterp());
  indexesForColumnsInterp.resize(lfi.getNCoefficientsForColumnsInterp(), lfi.getNOrderCoefficientsForColumnsInterp());
  coefficientsForLinesInterp = lfi.getCoefficientsForLinesInterp();
  coefficientsForColumnsInterp = lfi.getCoefficientsForColumnsInterp();
  indexesForLinesInterp = lfi.getIndexesForLinesInterp();
  indexesForColumnsInterp = lfi.getIndexesForColumnsInterp();
}

LagrangeFastInterpolator2D::~LagrangeFastInterpolator2D()
{
  coefficientsForLinesInterp.free();
  coefficientsForColumnsInterp.free();
  indexesForLinesInterp.free();
  indexesForColumnsInterp.free();
}

void LagrangeFastInterpolator2D::setLfi2D(const LagrangeFastInterpolator2D & lfi)
{
  coefficientsForLinesInterp.resize(lfi.getNCoefficientsForLinesInterp(), lfi.getNOrderCoefficientsForLinesInterp());
  coefficientsForColumnsInterp.resize(lfi.getNCoefficientsForColumnsInterp(), lfi.getNOrderCoefficientsForColumnsInterp());
  indexesForLinesInterp.resize(lfi.getNCoefficientsForLinesInterp(), lfi.getNOrderCoefficientsForLinesInterp());
  indexesForColumnsInterp.resize(lfi.getNCoefficientsForColumnsInterp(), lfi.getNOrderCoefficientsForColumnsInterp());
  coefficientsForLinesInterp = lfi.getCoefficientsForLinesInterp();
  coefficientsForColumnsInterp = lfi.getCoefficientsForColumnsInterp();
  indexesForLinesInterp = lfi.getIndexesForLinesInterp();
  indexesForColumnsInterp = lfi.getIndexesForColumnsInterp();
}

void interpolate2Dlfi(blitz::Array<std::complex<float>, 1> YInterp,
                      const blitz::Array<std::complex<float>, 1>& Y,
                      const LagrangeFastInterpolator2D& lfi2D)
{
  const blitz::Array<float, 2> coefficientsForLinesInterp(lfi2D.getCoefficientsForLinesInterp());
  const blitz::Array<float, 2> coefficientsForColumnsInterp(lfi2D.getCoefficientsForColumnsInterp());
  const blitz::Array<int, 2> indexesForLinesInterp(lfi2D.getIndexesForLinesInterp());
  const blitz::Array<int, 2> indexesForColumnsInterp(lfi2D.getIndexesForColumnsInterp());
  blitz::Array<std::complex<float>, 1> YInterpTmp(coefficientsForColumnsInterp.extent(0));
  YInterp = 0.0;
  YInterpTmp = 0.0;
  // interpolation following dimension 0
  for (int i=0 ; i<coefficientsForColumnsInterp.extent(0) ; i++) {
    for (int j=0 ; j<coefficientsForColumnsInterp.extent(1) ; j++) {
      YInterpTmp(i) += coefficientsForColumnsInterp(i, j) * Y(indexesForColumnsInterp(i, j));
    }
  }
  // interpolation following dimension 1
  for (int i=0 ; i<coefficientsForLinesInterp.extent(0) ; i++) {
    for (int j=0 ; j<coefficientsForLinesInterp.extent(1) ; j++) {
      YInterp(i) += coefficientsForLinesInterp(i, j) * YInterpTmp(indexesForLinesInterp(i, j));
    }
  }
}

void anterpolate2Dlfi(blitz::Array<std::complex<float>, 1> YAnterp,
                      const blitz::Array<std::complex<float>, 1>& Y,
                      const LagrangeFastInterpolator2D& lfi2D)
{
  const blitz::Array<float, 2> coefficientsForLinesInterp(lfi2D.getCoefficientsForLinesInterp());
  const blitz::Array<float, 2> coefficientsForColumnsInterp(lfi2D.getCoefficientsForColumnsInterp());
  const blitz::Array<int, 2> indexesForLinesInterp(lfi2D.getIndexesForLinesInterp());
  const blitz::Array<int, 2> indexesForColumnsInterp(lfi2D.getIndexesForColumnsInterp());
  blitz::Array<std::complex<float>, 1> YAnterpTmp(coefficientsForColumnsInterp.extent(0));
  YAnterp = 0.0;
  YAnterpTmp = 0.0;
  for (int i=0 ; i<coefficientsForLinesInterp.extent(0) ; i++) {
    for (int j=0 ; j<coefficientsForLinesInterp.extent(1) ; j++) {
      YAnterpTmp(indexesForLinesInterp(i, j)) += coefficientsForLinesInterp(i, j) * Y(i);
    }
  }
 for (int i=0 ; i<coefficientsForColumnsInterp.extent(0) ; i++) {
    for (int j=0 ; j<coefficientsForColumnsInterp.extent(1) ; j++) {
      YAnterp(indexesForColumnsInterp(i, j)) += coefficientsForColumnsInterp(i, j) * YAnterpTmp(i);
    }
  }
}

