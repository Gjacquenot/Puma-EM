#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <complex>
#include <blitz/array.h>

using namespace std;

/*void Lagrangian_regular_matrix(const int n,
                               blitz::Array<double, 2>& M);

void Lagrange_vector_fixedstep_interpolation(blitz::Array<complex<double>, 1>& y, 
                                             const blitz::Array<double, 1>& x, 
                                             const blitz::Array<double, 1>& xi, 
                                             const blitz::Array<complex<double>, 1>& yi, 
                                             const blitz::Array<double, 2>& M);

void decimate_2D (blitz::Array<complex<double>, 2> Y, 
                  const blitz::Array<complex<double>, 2>& Y_i, 
                  const blitz::Array<double, 1>& X1_i, // the abscissas following 1st dimension
                  const blitz::Array<double, 1>& X2_i, // the abscissas following 2nd dimension
                  const int n);
*/
class LagrangeFastInterpolator2D {
    blitz::Array<float, 2> coefficientsForLinesInterp;
    blitz::Array<float, 2> coefficientsForColumnsInterp;
    blitz::Array<int, 2> indexesForLinesInterp;
    blitz::Array<int, 2> indexesForColumnsInterp;

  public:
    // constructors and destructor
    LagrangeFastInterpolator2D(void);
    LagrangeFastInterpolator2D(const blitz::Array<float, 1>& /*x*/,
                               const blitz::Array<float, 1>& /*xi*/,
                               const float /*axi*/, // lim inf of xi interval
                               const float /*bxi*/, // lim sup of xi interval
                               const int /*INCLUDED_BOUNDARIES_xi*/,
                               const int /*NOrderLines*/,
                               const int /*PERIODIC_xi*/,
                               const int /*CYCLIC_xi*/,
                               const blitz::Array<float, 1>& /*y*/,
                               const blitz::Array<float, 1>& /*yi*/,
                               const float /*ayi*/, // lim inf of yi interval
                               const float /*byi*/, // lim sup of yi interval
                               const int /*INCLUDED_BOUNDARIES_yi*/,
                               const int /*NOrderColumns*/,
                               const int /*PERIODIC_yi*/,
                               const int /*CYCLIC_yi*/);
    LagrangeFastInterpolator2D(const LagrangeFastInterpolator2D &); // copy constructor
    ~LagrangeFastInterpolator2D();

    // functions
    int getNCoefficientsForLinesInterp() const {return coefficientsForLinesInterp.extent(0);}
    int getNOrderCoefficientsForLinesInterp() const {return coefficientsForLinesInterp.extent(1);}
    int getNCoefficientsForColumnsInterp() const {return coefficientsForColumnsInterp.extent(0);}
    int getNOrderCoefficientsForColumnsInterp() const {return coefficientsForColumnsInterp.extent(1);}
    blitz::Array<float, 2> getCoefficientsForLinesInterp() const {return coefficientsForLinesInterp;}
    blitz::Array<float, 2> getCoefficientsForColumnsInterp() const {return coefficientsForColumnsInterp;}
    blitz::Array<int, 2> getIndexesForLinesInterp() const {return indexesForLinesInterp;}
    blitz::Array<int, 2> getIndexesForColumnsInterp() const {return indexesForColumnsInterp;}
    void setLfi2D(const LagrangeFastInterpolator2D &);
};

void interpolate2Dlfi(blitz::Array<complex<float>, 1> Y_interp,
                      const blitz::Array<complex<float>, 1>& Y,
                      const LagrangeFastInterpolator2D& lfi2D);

void anterpolate2Dlfi(blitz::Array<complex<float>, 1> Y_anterp,
                      const blitz::Array<complex<float>, 1>& Y,
                      const LagrangeFastInterpolator2D& lfi2D);

#endif
