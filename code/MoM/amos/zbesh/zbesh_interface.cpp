#include <iostream>
#include <complex>
//#include <blitz/array.h>

/*%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*% interface to the actual (fortran) ZBESH routine
%*/

using namespace std;

const complex<double> I (0.0, 1.0);

extern "C" {
  //SUBROUTINE ZBESH(ZR, ZI, FNU, KODE, M, N, CYR, CYI, NZ, IERR)
  void zbesh_(double *, double *, double *, int *, int *, int *, double *, double *, int *, int *);
}

void zbesh(complex<double>& z_cmplx, double & fnu, int & kode, int & M, int & N, complex<double>& h, int & nz, int & ierr) {

  double zr = real(z_cmplx), zi = imag(z_cmplx), cyr, cyi;
  zbesh_(&zr, &zi, &fnu, &kode, &M, &N, &cyr, &cyi, &nz, &ierr);
  if (ierr>0) {
    if ((zr>0.0) || (zi>0.0)) cout << "Error in computation of zbesh. Error code: " << ierr << endl;
  }
  h = cyr + I * cyi;
}
