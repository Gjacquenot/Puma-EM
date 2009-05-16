#include <iostream>
#include <complex>
#include <blitz/array.h>

/*%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*% interface to the actual (fortran) ZTRSM routine
%*/

using namespace std;

const complex<double> I (0.0, 1.0);

extern "C" {
  //SUBROUTINE ZTRSM( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB )
  void ztrsm_(char *, char *, char *, char *, int *, int *, complex<double> *, complex<double> *, int *, complex<double> *, int *);
}

void ztrsm(char & side, char & uplo, char & transa, char & diag, int & m, int & n, complex<double> & alpha, blitz::Array<std::complex<double>, 2> & A, int & lda, blitz::Array<std::complex<double>, 2> & B, int & ldb) {
  ztrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A.data(), &lda, B.data(), &ldb);
}
