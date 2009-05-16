#include <fstream>
#include <iostream>
#include <complex>
#include <cmath>
#include <blitz/array.h>

using namespace blitz;

const complex<double> I (0.0, 1.0);

#include "Zmain.h"
#include "Vmain.h"

void write_Z_matrix (char * filename, Array<complex<double>, 2>& A) {

  int i, j, L = A.rows();
  ofstream ofs (filename);
  if (! ofs.is_open()) { 
    cout << "main : Error opening file : " << filename << endl; 
    exit (1);
  }
  else
  {
    ofs.precision(12);
    for (i=0 ; i<L ; i++) {
      for (j=0 ; j<L ; j++) ofs << real(A (i, j)) << "\t" << imag(A (i, j)) << "\t";
      ofs << "\n";
    }
  }
  ofs.close ();
}

void write_V_matrix (char * filename, Array<complex<double>, 2>& V) {

  int i, j, L = V.rows(), C = V.columns();
  ofstream ofs (filename);
  if (! ofs.is_open()) { 
    cout << "main : Error opening file : " << filename << endl; 
    exit (1);
  }
  else
  {
    ofs.precision(12);
    for (i=0 ; i<L ; i++) {
      for (j=0 ; j<C ; j++) ofs << real(V (i, j)) << "\t" << imag(V (i, j)) << "\t";
      ofs << "\n";
    }
  }
  ofs.close ();
}

int main (void) {
  Array<complex<double>, 2> Z_TE_J, Z_NE_J, Z_TH_J, Z_NH_J, Z_TE_M, Z_NE_M, Z_TH_M, Z_NH_M;
  Zmain (Z_TE_J, Z_NE_J, Z_TH_J, Z_NH_J, Z_TE_M, Z_NE_M, Z_TH_M, Z_NH_M);

  write_Z_matrix ("Z_TE_J.txt", Z_TE_J);
  write_Z_matrix ("Z_NE_J.txt", Z_NE_J);
  write_Z_matrix ("Z_TH_J.txt", Z_TH_J);
  write_Z_matrix ("Z_NH_J.txt", Z_NH_J);
  write_Z_matrix ("Z_TE_M.txt", Z_TE_M);
  write_Z_matrix ("Z_NE_M.txt", Z_NE_M);
  write_Z_matrix ("Z_TH_M.txt", Z_TH_M);
  write_Z_matrix ("Z_NH_M.txt", Z_NH_M);

  complex<double> S11_dip_soil;
  Array<complex<double>, 2> V_EJ, V_HJ, V_EM, V_HM, V_EJ_dip, V_HJ_dip, V_EM_dip, V_HM_dip;
  Vmain (V_EJ, V_HJ, V_EM, V_HM, V_EJ_dip, V_HJ_dip, V_EM_dip, V_HM_dip, S11_dip_soil);
  write_V_matrix ("V_EJ.txt", V_EJ);
  write_V_matrix ("V_HJ.txt", V_HJ);
  write_V_matrix ("V_EM.txt", V_EM);
  write_V_matrix ("V_HM.txt", V_HM);
  write_V_matrix ("V_EJ_dip.txt", V_EJ_dip);
  write_V_matrix ("V_HJ_dip.txt", V_HJ_dip);
  write_V_matrix ("V_EM_dip.txt", V_EM_dip);
  write_V_matrix ("V_HM_dip.txt", V_HM_dip);

  {
    FILE * f = fopen("S11_dip_soil.txt", "w");
    fprintf (f, "%.10g\t%.10g\t", real(S11_dip_soil), imag(S11_dip_soil));
    fclose (f);
  }

  return 0;
}
