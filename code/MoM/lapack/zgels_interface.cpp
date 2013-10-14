#include <iostream>
#include <complex>
#include <blitz/array.h>

/*%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*% interface to the actual (fortran) ZBESH routine
%*/

using namespace std;

extern "C" {
  //SUBROUTINE ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
  void zgels_(char *, int *, int *, int *, std::complex<double> *, int *, std::complex<double> *, int *, std::complex<double> *, int *, int *);
}

void zgels(char & trans, int & m, int & n, int & nrhs, blitz::Array<std::complex<double>, 2> & A, int & lda, blitz::Array<std::complex<double>, 2> & B, int & ldb, blitz::Array<std::complex<double>, 1> & work, int & lwork, int & info) {
  zgels_(&trans, &m, &n, &nrhs, A.data(), &lda, B.data(), &ldb, work.data(), &lwork, &info);
}

void zgels2(char & trans, int & m, int & n, int & nrhs, blitz::Array<std::complex<double>, 1> & A, int & lda, blitz::Array<std::complex<double>, 1> & B, int & ldb, blitz::Array<std::complex<double>, 1> & work, int & lwork, int & info) {
  zgels_(&trans, &m, &n, &nrhs, A.data(), &lda, B.data(), &ldb, work.data(), &lwork, &info);
}
