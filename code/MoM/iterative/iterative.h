#ifndef ITERATIVE_H
#define ITERATIVE_H

#include <iostream>
#include <fstream>
#include <string>
#include <blitz/array.h>
#include <complex>
#include <mpi.h>

using namespace std;

#include "../lapack/ztrsm_interface.h"

const complex<float> ZZERO(0.0, 0.0);

/** \file iterative.h
 *  The \c gmres routine uses the "functor" technique, which will allow 
 *  a user to define his own matrix-vector product to be used within \c gmres.
 *  Please note that all the functions are defined in \c gmres.h
 *
 *  There is also a \c classicalMatrixMatvec routine which will be the default matrix-vector
 *  multiplication routine used in \c gmres, so that it does not need to be written by the user.
 *
 *  The default way to use \c gmres is as follows:
 *
 *  \code
 *  {
 *    ....
 *    blitz::Array<T, 2> A(N,N); // T is the type of the blitz::Array
 *    ....
 *    A = some operations
 *    ....
 *    classicalMatrixMatvec<T> Aptr(A); // this class consists in (1) a pointer to A and (2) a matvec
 *                                      // routine which is simply "matrixMultiply" defined in gmres.h
 *    // we now create the functor "matvec" that encapsulates Aptr
 *    TSpecificFunctor<classicalMatrixMatvec<T>, T> matvec(&Aptr, &classicalMatrixMatvec<T>::matvec);
 *    // gmres call
 *    gmres(x, error, iter, flag, matvec, b, tol, restrt, max_it);
 *    ....
 *  }
 *  \endcode
 *
 *  If one wants to use its own matrix-vector multiplication routine
 *  (for example a custom sparse-matrix multiplication routine), he needs to do:
 *
 *  \code
 *  {
 *    .....
 *    sparseMatrix<T> Asparse(N,N); // this user-created object has its own matvec routine
 *    .....
 *    Asparse = some operations
 *    ....
 *    // we now create the functor "matvec" that encapsulates Asparse
 *    TSpecificFunctor<classicalMatrixMatvec<T>, T> matvec(&Asparse, &sparseMatrix<T>::matvec);
 *    // gmres call
 *    gmres(x, error, iter, flag, matvec, b, tol, restrt, max_it);
 *    ....
 *  }
 *  \endcode
 */

// redefinition of conjugate functions. This is necessary, since the "conj"
// built-in function only takes complex<T> and exits an error if its argument
// is real.
inline float conjScalar(const float & x)
{
  return x;
}

inline double conjScalar(const double & x)
{
  return x;
}

template <class T>
inline std::complex<T> conjScalar(const std::complex<T> & x)
{
  return conj(x);
}

inline blitz::Array<float, 1>
conjArray(const blitz::Array<float, 1>& x)
{
  return x;
}

inline blitz::Array<double, 1>
conjArray(const blitz::Array<double, 1>& x)
{
  return x;
}

template <class T>
inline blitz::Array<std::complex<T>, 1>
conjArray(const blitz::Array<std::complex<T>, 1>& x)
{
  blitz::Array<std::complex<T>, 1> y(x.size());
  y = conj(x);
  return y;
}

template <class T>
double norm2 (const blitz::Array<T, 1>& x)
{
  return sqrt(abs(sum(x * conjArray(x))));
}

template <class T>
double squareNorm2 (const blitz::Array<T, 1>& x)
{
  return real(sum(x * conjArray(x)));
}

template <class T>
void checkMatrixVectorDimensions(const int M, /**< required number of A lines */
                                 const int N, /**< required number of A columns */
                                 const blitz::Array<T, 2>& A, /**< a 2D array */
                                 const blitz::Array<T, 1>& x) /**< the vector */
/** check dimensions in view of performing \f$ A x \f$ or solving \f$ A x = y \f$ */
{
  if ( A.extent(0) != M ) {
    std::cout << "checkMatrixVectorDimensions():" << std::endl;
    std::cout << "problem with A dimensions: A.extent(0) = " << A.extent(0) 
              << std::endl;
    std::cout << "M = " << M << std::endl;
    exit(1);
  }
  if ( A.extent(1) != N ) {
    std::cout << "checkMatrixVectorDimensions():" << std::endl;
    std::cout << "problem with A dimensions: A.extent(1) = " << A.extent(1)
              << std::endl;
    std::cout << "N = " << N << std::endl;
    exit(1);
  }
  if ( x.extent(0) != N ) {
    std::cout << "checkMatrixVectorDimensions():" << std::endl;
    std::cout << "problem with x dimensions: x.extent(0) = " << x.extent(0)
              << std::endl;
    std::cout << "N = " << N << std::endl;
    exit(1);
  }
}

template <class T>
blitz::Array<T, 1> matrixMultiply (const blitz::Array<T, 2>& A, /**< a 2D array */
                                   const blitz::Array<T, 1>& x) /**< the vector */
/** performs the matrix-vector multiplication \f$ A x \f$ */
{
  blitz::Range all = blitz::Range::all();
  int i;
  const int M = A.extent(0), N = x.extent(0);
  // check dimensions
  checkMatrixVectorDimensions(M, N, A, x);
  blitz::Array<T, 1> result(M);
  result = 0.0;
  for (i=0 ; i<M ; i++) result(i) += sum(A(i, all) * x(all));
  return result;
}

template <class T>
void triangleUpSolve(blitz::Array<T, 1>& x, /**< output 1D array */
                     const blitz::Array<T, 2>& U, /**< upper triangular matrix */
                     const blitz::Array<T, 1>& y) /**< input 1D array */
/**
 *  This function computes the solution of a linear system
 *
 *     \f[ U x = y \f]
 *
 *  where \f$ U \f$ is a upper triangular matrix, that is, \f$ U(i>j, j) = 0 \f$.
 */
{
  int N = y.extent(0);
  // check dimensions: U must be square
  checkMatrixVectorDimensions(N, N, U, y);
  for (int i=N-1 ; i>-1 ; i--) {
    // x(i) = 1/U(i, i) * ( y(i) - sum(U(i, i+1:N-1) * x(i+1:N-1)) )
    x(i) = y(i);
    for (int j=i+1 ; j<N ; j++) {
      x(i) -= U(i, j) * x(j);
    }
    x(i) /= U(i, i);
  }
  // use of blas level 3 ztrsm routine...NOT WORKING YET
/*  blitz::Array<std::complex<double>, 2> A(N, N, fortranArray);
  for (int i=0 ; i<N ; ++i) {
    for (int j=0 ; j<N ; ++j) A(i, j) = U(i, j);
  }
  char side = 'L', uplo = 'U', transa = 'N', diag = 'N';
  std::complex<double> alpha(1.0, 0.0);
  blitz::Array<std::complex<double>, 2> B(shape(N,N), fortranArray);
  B = 0.0;
  for (int i=0 ; i<N ; ++i) B(i, i) = 1.;
  // we use the ztrsm blas level 3 routine
  //ztrsm(side, uplo, transa, diag, N, N, alpha, A, N, B, N);
  blitz::Array<std::complex<double>, 2> E(N, N, fortranArray);
  blitz::Range all = blitz::Range::all();
  for (int i=0 ; i<N ; ++i) {
    for (int j=0 ; j<N ; ++j) E(i, j) = sum(A(i, all) * B(all, j));
  }
    int ierror, num_procs = MPI::COMM_WORLD.Get_size(), my_id = MPI::COMM_WORLD.Get_rank();
  if (my_id==0) {
    cout << "OK HERE " << endl;
    //cout << "E = " << E << endl;
    cout << "A = " << A << endl;
    //cout << "B = " << B << endl;
  }
  flush(cout);
    ierror = MPI_Barrier(MPI::COMM_WORLD);

/*  x = 0.0;
  for (int i=0 ; i<N ; i++) {
    for (int j=0 ; j<N ; ++j) x(i) += static_cast< std::complex<float> >(B(i, j)) * y(j);
  }*/
}

template <typename T>
void rotmat(T & cs,
            T & sn,
            T & a,
            T & b)
/** compute the Givens rotation matrix parameters for a and b */
{
  const T ONE = static_cast<T>(1.0);
  T temp;
  if (abs(a)==0.0) {
    cs = 0.0;
    sn = 1.0;
    a = b;
  }
  T SCALE = abs(a) + abs(b);
  T NORM = SCALE*sqrt(pow2((abs(a/SCALE)))+pow2((abs(b/SCALE))));
  T ALPHA = a/abs(a);
  cs = abs(a)/NORM;
  sn = ALPHA*conj(b)/NORM;
  a = ALPHA*NORM;

/*  else if ( abs(b) > abs(a) ) {
    temp = -a/b;
    sn = ONE / sqrt(ONE + temp*temp);
    cs = temp * sn;
  }
  else {
    temp = -b/a;
    cs = ONE / sqrt(ONE + temp*temp);
    sn = temp * cs;
  }*/
}

// We make use of functor mechanism for letting the user define its own
// matrix-vector multiplication method. See http://www.newty.de/fpt/functor.html

/** we give a default class that is a classical matrix-vector multiplication */
template <typename T>
class classicalMatrixMatvec
{
    const blitz::Array<T, 2> * matPtr; /**< a pointer towards a 2D \c blitz Array  */

  public:
    /// constructor
    classicalMatrixMatvec(const blitz::Array<T, 2>& A) {matPtr = &A;};
    ~classicalMatrixMatvec(void){};
    /// matvec function
    blitz::Array<T, 1> classicalMatvec(const blitz::Array<T, 1>& x) {return matrixMultiply(*matPtr, x);};
};

// matvec functor
template <class T, class TClassA = classicalMatrixMatvec<T> >
class MatvecFunctor
{
  private:
    blitz::Array<T, 1> (TClassA::*funcpt)(const blitz::Array<T, 1>&); /**< pointer to member function */
    TClassA* pt2Object; /**< pointer to object */

  public:
    /** constructor - takes pointer to an object and pointer to a member
      * and stores them in the two private variables
      */
    MatvecFunctor(TClassA* _pt2Object, blitz::Array<T, 1>(TClassA::*_funcpt)(const blitz::Array<T, 1>&)) { pt2Object = _pt2Object; funcpt = _funcpt; };
    ~MatvecFunctor(void){};

    // operator ()
    blitz::Array<T, 1> operator()(const blitz::Array<T, 1>& x)
    { return (*pt2Object.*funcpt)(x); }; // execute member function
};

/** we give a default class that corresponds to point-diagonal preconditioning */
template <typename T>
class DiagPrecond
{
    const blitz::Array<T, 1> * precondPtr; /**< a pointer towards a 1D \c blitz Array  */

  public:
    /// constructor
    DiagPrecond(const blitz::Array<T, 1>& M) {precondPtr = &M;};
    ~DiagPrecond(void){};
    /// noPrecond solution
    blitz::Array<T, 1> precondSolve(const blitz::Array<T, 1>& x) {return *precondPtr * x;};
};

// precond functor
template <class T, class TClassB = DiagPrecond<T> >
class PrecondFunctor
{
  private:
    blitz::Array<T, 1> (TClassB::*funcpt)(const blitz::Array<T, 1>&); /**< pointer to member function */
    TClassB* pt2Object; /**< pointer to object */

  public:
    /** constructor - takes pointer to an object and pointer to a member
      * and stores them in the two private variables
      */
    PrecondFunctor(TClassB* _pt2Object, blitz::Array<T, 1>(TClassB::*_funcpt)(const blitz::Array<T, 1>&)) { pt2Object = _pt2Object; funcpt = _funcpt; };
    ~PrecondFunctor(void){};

    // operator ()
    blitz::Array<T, 1> operator()(const blitz::Array<T, 1>& x)
    { return (*pt2Object.*funcpt)(x); }; // execute member function
};

// left preconditioned GMRES
template <typename T, typename TClassA, typename TClassB>
void gmres(blitz::Array<T, 1>& x, /**< OUTPUT: converged solution */
           double & error, /**< OUTPUT: the error */
           int & iter, /**< OUTPUT: number of iterations needed */
           int & flag, /**< OUTPUT: success flag: 0 if OK */
           MatvecFunctor<T, TClassA> matvec, /**< INPUT: matvec functor */
           PrecondFunctor<T, TClassB> psolve, /**< INPUT: precond functor */
           const blitz::Array<T, 1>& b, /**< INPUT: right-hand side */
           const double tol, /**< INPUT: tolerance on solution */
           const int RESTRT, /**< INPUT: restart number */
           const int MAXITER, /**< INPUT: max number of iterations */
           const int my_id, /**< INPUT: the process ID */
           const int num_proc, /**< INPUT: the number of processes */
           const string convergenceDetailedOutput)
{
  std::ofstream ofs (convergenceDetailedOutput.c_str());
  if (! ofs.is_open()) { 
    cout << "error opening " << convergenceDetailedOutput << endl; 
    exit(1);
  }
  ofs.precision(8);
  ofs << "# GMRES algorithm" << endl;
  ofs << "# output showing the convergence for tol = " << tol << endl;

  blitz::Range all = blitz::Range::all();
  flag = 0;
  iter = 0;
  int ierror;
  bool is_return = false;
  const int N_local = x.extent(0), m = RESTRT; // size of the system
  // dimensions and other checks
  if (RESTRT < 1) {
    std::cout << "Bad restart value. RESTRT = " << RESTRT << std::endl;
    exit(1);
  }
  if (MAXITER < 1) {
    std::cout << "Bad maxiter value. MAXITER = " << MAXITER << std::endl;
    exit(1);
  }

  double bnorm2, local_bnorm2, rnorm2, local_rnorm2, wnorm2, local_wnorm2;
  local_bnorm2 = squareNorm2(b);
  ierror = MPI_Allreduce(&local_bnorm2, &bnorm2, 1, MPI::DOUBLE, MPI::SUM, MPI::COMM_WORLD);
  bnorm2 = sqrt(abs(bnorm2));
  if (bnorm2==0.0) bnorm2 = 1.0;

  T temp;
  blitz::Array<T, 1> rTmp(N_local), r(N_local);

  // workspaces definitions
  blitz::Array<T, 1> cs(blitz::Range(1, m)), sn(blitz::Range(1, m)), s(m+1), wTmp(N_local), w(N_local), y;
  blitz::Array<T, 2> V(N_local, m+1), H(blitz::Range(1, m+1), blitz::Range(1, m+1));

  for (iter=0 ; iter<MAXITER ; iter++) {
    rTmp = b - matvec(x);
    r = psolve(rTmp);
    local_rnorm2 = squareNorm2(r);
    ierror = MPI_Allreduce(&local_rnorm2, &rnorm2, 1, MPI::DOUBLE, MPI::SUM, MPI::COMM_WORLD);
    rnorm2 = sqrt(abs(rnorm2));
    V(all, 0) = r/static_cast<T>(rnorm2);
    s = 0.0;
    s(0) = static_cast<T>(rnorm2);

    H = 0.0;
    H(1, m+1) = rnorm2;

    for (int jH=1 ; jH<m+1 ; ++jH){
      // if preconditioning: w = M^(-1) * (A*V(all, jH))
      wTmp = matvec(V(all, jH-1));
      w = psolve(wTmp);
      for (int j=1 ; j<jH+1 ; ++j) H(j, jH) = 0.0;

      // construct orthonormal basis using Modified Gram-Schmidt
      T dloo = 0.0;
      for (int j=1 ; j<jH+1 ; ++j) {
        complex<double> H_local = sum(w * conjArray(V(all, j-1))), H_global;
        ierror = MPI_Allreduce(&H_local, &H_global, 1, MPI::DOUBLE_COMPLEX, MPI::SUM, MPI::COMM_WORLD);
        H(j, jH) = H_global;
        w -= H(j, jH) * V(all, j-1);
        dloo += abs(H_global)*abs(H_global);
      }
      dloo = sqrt(dloo);
      local_wnorm2 = squareNorm2(w);
      ierror = MPI_Allreduce(&local_wnorm2, &wnorm2, 1, MPI::DOUBLE, MPI::SUM, MPI::COMM_WORLD);
      wnorm2 = sqrt(abs(wnorm2));

      H(jH+1, jH) = wnorm2;
      if (jH<m) V(all, jH) = w / static_cast<T>(wnorm2);

      // now we apply the GIVENS rotations
      for (int j=1 ; j<jH ; ++j) {
        temp = cs(j) * H(j,jH) + conjScalar(sn(j)) * H(j+1,jH);
        H(j+1,jH) = cs(j) * H(j+1,jH) - conjScalar(sn(j)) * H(j,jH);
        H(j,jH) = temp;
      }

      // form jH-th rotation matrix
      T auxHjj = H(jH,jH), auxHjp1j= H(jH+1,jH), temp;
      rotmat(temp, sn(jH), auxHjj, auxHjp1j);
      cs(jH)= temp;

      temp = cs(jH) * H(jH, m+1) + conjScalar(sn(jH)) * H(jH+1, m+1);
      H(jH+1, m+1) = cs(jH) * H(jH+1, m+1) - conjScalar(sn(jH)) * H(jH, m+1);
      H(jH, m+1) = temp;

      temp = cs(jH) * H(jH, jH) + conjScalar(sn(jH)) * H(jH+1, jH);
      H(jH+1, jH) = 0.0;
      H(jH, jH) = temp;

      // approximate residual norm
      s(jH) = H(jH+1, m+1);
      s(jH-1) = H(jH, m+1);
      error = abs(s(jH)) / bnorm2;
      /*error = abs(H(jH+1, m+1));*/
      ofs << error << endl;

      // update approximation x
      if ( error<=tol ) {
        y.resize(jH); // Range: 0..jH-1
        blitz::Array<T, 2> H2(jH, jH);
        H2 = H(blitz::Range(1, jH), blitz::Range(1, jH));
        triangleUpSolve( y, H2, s(blitz::Range(0, jH-1)) );
        x += matrixMultiply(V(all, blitz::Range(0, jH-1)), y);
        ofs.close();
        return;
      }
    } // end for (jH =...)
    if ( error<=tol ) {
      ofs.close();
      return;
    }
    y.resize(m); // Range: 0..m-1
    blitz::Array<T, 2> H2(m, m);
    H2 = H(blitz::Range(1, m), blitz::Range(1, m));
    triangleUpSolve( y, H2, s(blitz::Range(0, m-1)) );
    x += matrixMultiply(V(all, blitz::Range(0, m-1)), y);
    rTmp = b - matvec(x);
    r = psolve(rTmp);
    local_rnorm2 = squareNorm2(r);
    ierror = MPI_Allreduce(&local_rnorm2, &rnorm2, 1, MPI::DOUBLE, MPI::SUM, MPI::COMM_WORLD);
    rnorm2 = sqrt(abs(rnorm2));
    s(m) = abs(rnorm2);
    // check convergence
    error = abs(s(m)) / bnorm2;
    ofs << "intermediate error " << error << endl;
    if ( error<=tol ) {
      ofs.close();
      return;
    }
  } // end for (iter =...)

  // bad ending...
  if (error>tol) {
    ofs.close();
    flag = 1;
  }
}

// right preconditioned or flexible GMRES
template <typename T, typename TClassA, typename TClassB>
void fgmres(blitz::Array<T, 1>& x, /**< OUTPUT: converged solution */
            double & error, /**< OUTPUT: the error */
            int & iter, /**< OUTPUT: number of iterations needed */
            int & flag, /**< OUTPUT: success flag: 0 if OK */
            MatvecFunctor<T, TClassA> matvec, /**< INPUT: matvec functor */
            PrecondFunctor<T, TClassB> psolve, /**< INPUT: precond functor */
            const blitz::Array<T, 1>& b, /**< INPUT: right-hand side */
            const double tol, /**< INPUT: tolerance on solution */
            const int RESTRT, /**< INPUT: restart number */
            const int MAXITER, /**< INPUT: max number of iterations */
            const int my_id, /**< INPUT: the process ID */
            const int num_proc, /**< INPUT: the number of processes */
            const string convergenceDetailedOutput)
{
  std::ofstream ofs (convergenceDetailedOutput.c_str());
  if (! ofs.is_open()) { 
    cout << "error opening " << convergenceDetailedOutput << endl; 
    exit(1);
  }
  ofs.precision(8);
  ofs << "# FGMRES algorithm" << endl;
  ofs << "# output showing the convergence for tol = " << tol << endl;

  blitz::Range all = blitz::Range::all();
  flag = 0;
  iter = 0;
  int ierror;
  bool is_return = false;
  const int N_local = x.extent(0), m = RESTRT; // size of the system
  // dimensions and other checks
  if (RESTRT < 1) {
    std::cout << "Bad restart value. RESTRT = " << RESTRT << std::endl;
    exit(1);
  }
  if (MAXITER < 1) {
    std::cout << "Bad maxiter value. MAXITER = " << MAXITER << std::endl;
    exit(1);
  }

  double bnorm2, local_bnorm2, rnorm2, local_rnorm2, wnorm2, local_wnorm2;
  local_bnorm2 = squareNorm2(b);
  ierror = MPI_Allreduce(&local_bnorm2, &bnorm2, 1, MPI::DOUBLE, MPI::SUM, MPI::COMM_WORLD);
  bnorm2 = sqrt(abs(bnorm2));
  if (bnorm2==0.0) bnorm2 = 1.0;

  T temp;
  blitz::Array<T, 1> rTmp(N_local), r(N_local);

  // workspaces definitions
  blitz::Array<T, 1> cs(blitz::Range(1, m)), sn(blitz::Range(1, m)), s(m+1), wTmp(N_local), w(N_local), y;
  blitz::Array<T, 2> V(N_local, m+1), Z(N_local, m+1), H(blitz::Range(1, m+1), blitz::Range(1, m+1));

  for (iter=0 ; iter<MAXITER ; iter++) {
    r = b - matvec(x);
    local_rnorm2 = squareNorm2(r);
    ierror = MPI_Allreduce(&local_rnorm2, &rnorm2, 1, MPI::DOUBLE, MPI::SUM, MPI::COMM_WORLD);
    rnorm2 = sqrt(abs(rnorm2));
    V(all, 0) = r/static_cast<T>(rnorm2);
    s = 0.0;
    s(0) = static_cast<T>(rnorm2);

    H = 0.0;
    H(1, m+1) = rnorm2;

    for (int jH=1 ; jH<m+1 ; ++jH){
      // if right preconditioning: z = M^-1 * V(all, jH)
      Z(all, jH-1) = psolve(V(all, jH-1));
      w = matvec(Z(all, jH-1));

      for (int j=1 ; j<jH+1 ; ++j) H(j, jH) = 0.0;

      // construct orthonormal basis using Modified Gram-Schmidt
      T dloo = 0.0;
      for (int j=1 ; j<jH+1 ; ++j) {
        complex<double> H_local = sum(w * conjArray(V(all, j-1))), H_global;
        ierror = MPI_Allreduce(&H_local, &H_global, 1, MPI::DOUBLE_COMPLEX, MPI::SUM, MPI::COMM_WORLD);
        H(j, jH) = H_global;
        w -= H(j, jH) * V(all, j-1);
        dloo += abs(H_global)*abs(H_global);
      }
      dloo = sqrt(dloo);
      local_wnorm2 = squareNorm2(w);
      ierror = MPI_Allreduce(&local_wnorm2, &wnorm2, 1, MPI::DOUBLE, MPI::SUM, MPI::COMM_WORLD);
      wnorm2 = sqrt(abs(wnorm2));

      H(jH+1, jH) = wnorm2;
      if (jH<m) V(all, jH) = w / static_cast<T>(wnorm2);

      // now we apply the GIVENS rotations
      for (int j=1 ; j<jH ; ++j) {
        temp = cs(j) * H(j,jH) + conjScalar(sn(j)) * H(j+1,jH);
        H(j+1,jH) = cs(j) * H(j+1,jH) - conjScalar(sn(j)) * H(j,jH);
        H(j,jH) = temp;
      }

      // form jH-th rotation matrix
      T auxHjj = H(jH,jH), auxHjp1j= H(jH+1,jH), temp;
      rotmat(temp, sn(jH), auxHjj, auxHjp1j);
      cs(jH)= temp;

      temp = cs(jH) * H(jH, m+1) + conjScalar(sn(jH)) * H(jH+1, m+1);
      H(jH+1, m+1) = cs(jH) * H(jH+1, m+1) - conjScalar(sn(jH)) * H(jH, m+1);
      H(jH, m+1) = temp;

      temp = cs(jH) * H(jH, jH) + conjScalar(sn(jH)) * H(jH+1, jH);
      H(jH+1, jH) = 0.0;
      H(jH, jH) = temp;

      // approximate residual norm
      s(jH) = H(jH+1, m+1);
      s(jH-1) = H(jH, m+1);
      error = abs(s(jH)) / bnorm2;
      /*error = abs(H(jH+1, m+1));*/
      ofs << error << endl;

      // update approximation x
      if ( error<=tol ) {
        y.resize(jH); // Range: 0..jH-1
        blitz::Array<T, 2> H2(jH, jH);
        H2 = H(blitz::Range(1, jH), blitz::Range(1, jH));
        triangleUpSolve( y, H2, s(blitz::Range(0, jH-1)) );
        x += matrixMultiply(Z(all, blitz::Range(0, jH-1)), y);
        ofs.close();
        return;
      }
    } // end for (jH =...)
    if ( error<=tol ) {
      ofs.close();
      return;
    }
    y.resize(m); // Range: 0..m-1
    blitz::Array<T, 2> H2(m, m);
    H2 = H(blitz::Range(1, m), blitz::Range(1, m));
    triangleUpSolve( y, H2, s(blitz::Range(0, m-1)) );
    x += matrixMultiply(Z(all, blitz::Range(0, m-1)), y);
    r = b - matvec(x);
    local_rnorm2 = squareNorm2(r);
    ierror = MPI_Allreduce(&local_rnorm2, &rnorm2, 1, MPI::DOUBLE, MPI::SUM, MPI::COMM_WORLD);
    rnorm2 = sqrt(abs(rnorm2));
    s(m) = abs(rnorm2);
    // check convergence
    error = abs(s(m)) / bnorm2;
    ofs << "intermediate error " << error << endl;
    if ( error<=tol ) {
      ofs.close();
      return;
    }
  } // end for (iter =...)

  // bad ending...
  if (error>tol) {
    ofs.close();
    flag = 1;
  }
}

// BICGSTAB
template <typename T, typename TClassA, typename TClassB>
void bicgstab(blitz::Array<T, 1>& x, /**< OUTPUT: converged solution */
              double & error, /**< OUTPUT: the error */
              int & iter, /**< OUTPUT: number of iterations needed */
              int & flag, /**< OUTPUT: success flag: 0 if OK */
              MatvecFunctor<T, TClassA> matvec, /**< INPUT: matvec functor */
              PrecondFunctor<T, TClassB> psolve, /**< INPUT: precond functor */
              const blitz::Array<T, 1>& b, /**< INPUT: right-hand side */
              const double tol, /**< INPUT: tolerance on solution */
              const int MAXITER,  /**< INPUT: max number of iterations */
              const int my_id, /**< INPUT: the process ID */
              const int num_proc, /**< INPUT: the number of processes */
              const string convergenceDetailedOutput)
{
  std::ofstream ofs (convergenceDetailedOutput.c_str());
  if (! ofs.is_open()) { 
    cout << "error opening " << convergenceDetailedOutput << endl; 
    exit(1);
  }
  ofs.precision(8);
  ofs << "# BiCGSTAB algorithm" << endl;
  ofs << "# output showing the convergence for tol = " << tol << endl;

  blitz::Range all = blitz::Range::all();
  flag = 0;
  iter = 0;
  int ierror;
  bool is_return = false;
  const int N_local = x.extent(0); // size of the system
  // dimensions and other checks
  if (MAXITER < 1) {
    std::cout << "BiCGSTAB: Bad maxiter value. MAXITER = " << MAXITER << std::endl;
    exit(1);
  }
  T alpha, omega, rho, rho_local, rho_1;
  double bnorm2, snorm2, local_snorm2, local_bnorm2, local_error;
  // local arrays
  blitz::Array<T, 1> p(N_local), p_hat(N_local), r_tld(N_local), r(N_local), v(N_local), s(N_local), s_hat(N_local), t(N_local);
  local_bnorm2 = squareNorm2(b);
  ierror = MPI_Allreduce(&local_bnorm2, &bnorm2, 1, MPI::DOUBLE, MPI::SUM, MPI::COMM_WORLD);
  bnorm2 = sqrt(abs(bnorm2));
  if (bnorm2==0.0) bnorm2 = 1.0;
  v = matvec(x); // v is used as a temporary vector here
  r = b - v;
  local_error = squareNorm2(r);
  ierror = MPI_Allreduce(&local_error, &error, 1, MPI::DOUBLE, MPI::SUM, MPI::COMM_WORLD);
  error = sqrt( abs(error) ) / bnorm2;
  ofs << error << endl;

  if ( error <= tol ) {ofs.close(); return;}
  omega  = 1.0;
  r_tld = r;
  // beginning of the iterations
  for (iter=1 ; iter<=MAXITER ; iter++) {
    rho_local = sum(conjArray(r_tld) * r);
    ierror = MPI_Allreduce(&rho_local, &rho, 1, MPI::COMPLEX, MPI::SUM, MPI::COMM_WORLD);
    
    if (rho==ZZERO) {ofs.close(); return;}
    if (iter==1) p = r;
    else{
      T beta  = ( rho/rho_1 )*( alpha/omega );
      p = r + beta*( p - omega*v );
    }
    ierror = MPI_Barrier(MPI::COMM_WORLD);
    p_hat = psolve(p);
    v = matvec(p_hat);
    T alpha_denom_local, alpha_denom;
    alpha_denom_local = sum( conjArray(r_tld) * v );
    ierror = MPI_Allreduce(&alpha_denom_local, &alpha_denom, 1, MPI::COMPLEX, MPI::SUM, MPI::COMM_WORLD);
    alpha = rho / alpha_denom;
    s = r - alpha*v;
    local_snorm2 = squareNorm2(s);
    ierror = MPI_Allreduce(&local_snorm2, &snorm2, 1, MPI::DOUBLE, MPI::SUM, MPI::COMM_WORLD);
    snorm2 = sqrt(abs(snorm2));
    if ( snorm2 / bnorm2 < tol ) {
      x += alpha*p_hat;
      error = snorm2 / bnorm2;
      ofs << error << endl;
      ofs.close();
      return;
    }
    // stabilizer
    s_hat = psolve(s);
    t = matvec(s_hat);
    T local_omega_num, omega_num, local_omega_denom, omega_denom;
    local_omega_num = sum(conjArray(t) * s);
    local_omega_denom = sum(conjArray(t) * t);
    ierror = MPI_Allreduce(&local_omega_num, &omega_num, 1, MPI::COMPLEX, MPI::SUM, MPI::COMM_WORLD);
    ierror = MPI_Allreduce(&local_omega_denom, &omega_denom, 1, MPI::COMPLEX, MPI::SUM, MPI::COMM_WORLD);
    omega = omega_num/omega_denom;
    x += alpha*p_hat + omega*s_hat;
    r = s - omega*t;
    local_error = squareNorm2(r);
    ierror = MPI_Allreduce(&local_error, &error, 1, MPI::DOUBLE, MPI::SUM, MPI::COMM_WORLD);
    error = sqrt( abs(error) ) / bnorm2;
    ofs << error << endl;

    if ( error <= tol ) {ofs.close(); return;}
    if ( omega == ZZERO ) {ofs.close(); return;}
    rho_1 = rho;
  } // end for
  local_snorm2 = squareNorm2(s);
  ierror = MPI_Allreduce(&local_snorm2, &snorm2, 1, MPI::DOUBLE, MPI::SUM, MPI::COMM_WORLD);
  snorm2 = sqrt(abs(snorm2));
  if ( ( error <= tol ) || ( norm2( s ) <= tol ) ) { // converged
    if ( snorm2 / bnorm2 <= tol ) {
      error = snorm2 / bnorm2;
      ofs << error << endl;
    }
    flag = 0;
  }
  else if ( omega == ZZERO ) flag = -2; // breakdown
  else if ( rho == ZZERO ) flag = -1; // breakdown
  else flag = 1; // no convergence
  ofs.close();
  return;
}

#endif
