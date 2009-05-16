

using namespace std;

void zgels(char & trans, int & m, int & n, int & nrhs, blitz::Array<std::complex<double>, 2> & A, int & lda, blitz::Array<std::complex<double>, 2> & B, int & ldb, blitz::Array<std::complex<double>, 1> & work, int & lwork, int & info);
