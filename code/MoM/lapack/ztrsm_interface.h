

using namespace std;

void ztrsm(char & side, char & uplo, char & transa, char & diag, int & m, int & n, complex<double> & alpha, blitz::Array<std::complex<double>, 2> & A, int & lda, blitz::Array<std::complex<double>, 2> & B, int & ldb);
