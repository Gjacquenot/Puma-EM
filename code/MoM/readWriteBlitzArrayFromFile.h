#ifndef READWRITEBLITZARRAYFROMFILE_H
#define READWRITEBLITZARRAYFROMFILE_H
#include <iostream>
#include <fstream>
#include <string>
#include <blitz/array.h>

using namespace std;

std::string intToString(const int & x);
double stringToDouble(const std::string& s);
// write scalar to text file
template <typename T>
void writeScalarToASCIIFile(const string filename, const T & x);
void writeIntToASCIIFile(const string filename, const int & x);
void writeFloatToASCIIFile(const string filename, const float & x);
void writeDoubleToASCIIFile(const string filename, const double & x);
void writeComplexFloatToASCIIFile(const string filename, const complex<float> & x);
void writeComplexDoubleToASCIIFile(const string filename, const complex<double> & x);

// read string from file
void readStringFromASCIIFile(const string filename, string & x);
// read scalar from text file
template <typename T>
void readScalarFromASCIIFile(const string filename, T & x);
void readIntFromASCIIFile(const string filename, int & x);
void readFloatFromASCIIFile(const string filename, float & x);
void readDoubleFromASCIIFile(const string filename, double & x);
void readComplexFloatFromASCIIFile(const string filename, complex<float> & x);
void readComplexDoubleFromASCIIFile(const string filename, complex<double> & x);

// read ASCII 1-D arrays
template<typename T>
void readBlitzArray1DFromASCIIFile(const string filename, blitz::Array<T, 1>& A);
void readIntBlitzArray1DFromASCIIFile(const string filename, blitz::Array<int, 1>& A);
void readFloatBlitzArray1DFromASCIIFile(const string filename, blitz::Array<float, 1>& A);
void readDoubleBlitzArray1DFromASCIIFile(const string filename, blitz::Array<double, 1>& A);
void readComplexFloatBlitzArray1DFromASCIIFile(const string filename, blitz::Array<std::complex<float>, 1>& A);
void readComplexDoubleBlitzArray1DFromASCIIFile(const string filename, blitz::Array<std::complex<double>, 1>& A);

// read binary 1-D arrays
void readIntBlitzArray1DFromBinaryFile(const string filename, blitz::Array<int, 1>& A);
void readFloatBlitzArray1DFromBinaryFile(const string filename, blitz::Array<float, 1>& A);
void readDoubleBlitzArray1DFromBinaryFile(const string filename, blitz::Array<double, 1>& A);
void readComplexFloatBlitzArray1DFromBinaryFile(const string filename, blitz::Array<std::complex<float>, 1>& A);
void readComplexDoubleBlitzArray1DFromBinaryFile(const string filename, blitz::Array<std::complex<double>, 1>& A);

// write ASCII and Binary 1-D arrays
void writeIntBlitzArray1DToASCIIFile(const string filename, blitz::Array<int, 1>& A);
void writeComplexFloatBlitzArray1DToASCIIFile(const string filename, const blitz::Array<complex<float>, 1>& A);
void writeFloatBlitzArray1DToASCIIFile(const string filename, const blitz::Array<float, 1>& A);
void writeComplexFloatBlitzArray1DToBinaryFile(const string filename, const blitz::Array<complex<float>, 1>& A);
void writeDoubleBlitzArray1DToBinaryFile(const string filename, blitz::Array<double, 1>& A);
void writeIntBlitzArray1DToBinaryFile(const string filename, blitz::Array<int, 1>& A);


// read ASCII 2-D arrays
template<typename T>
void readBlitzArray2DFromASCIIFile(const string filename, blitz::Array<T, 2>& A);
void readIntBlitzArray2DFromASCIIFile(const string filename, blitz::Array<int, 2>& A);
void readFloatBlitzArray2DFromASCIIFile(const string filename, blitz::Array<float, 2>& A);
void readDoubleBlitzArray2DFromASCIIFile(const string filename, blitz::Array<double, 2>& A);
void readComplexFloatBlitzArray2DFromASCIIFile(const string filename, blitz::Array<std::complex<float>, 2>& A);
void readComplexDoubleBlitzArray2DFromASCIIFile(const string filename, blitz::Array<std::complex<double>, 2>& A);

// read binary 2-D arrays
void readIntBlitzArray2DFromBinaryFile(const string filename, blitz::Array<int, 2>& A);
void readFloatBlitzArray2DFromBinaryFile(const string filename, blitz::Array<float, 2>& A);
void readDoubleBlitzArray2DFromBinaryFile(const string filename, blitz::Array<double, 2>& A);
void readComplexFloatBlitzArray2DFromBinaryFile(const string filename, blitz::Array<std::complex<float>, 2>& A);
void readComplexDoubleBlitzArray2DFromBinaryFile(const string filename, blitz::Array<std::complex<double>, 2>& A);

// write ASCII 2-D array
void writeIntBlitzArray2DToASCIIFile(const string filename, blitz::Array<int, 2>& A);
void writeComplexFloatBlitzArray2DToASCIIFile(const string filename, blitz::Array<complex<float>, 2>& A);
void writeComplexDoubleBlitzArray2DToASCIIFile(const string filename, blitz::Array<complex<double>, 2>& A);
void writeDoubleBlitzArray2DToASCIIFile(const string filename, blitz::Array<double, 2>& A);
void writeFloatBlitzArray2DToASCIIFile(const string filename, blitz::Array<float, 2>& A);
// write 2D binary arrays
void writeIntBlitzArray2DToBinaryFile(const string filename, blitz::Array<int, 2>& A);
void writeFloatBlitzArray2DToBinaryFile(const string filename, blitz::Array<float, 2>& A);
void writeDoubleBlitzArray2DToBinaryFile(const string filename, blitz::Array<double, 2>& A);
void writeComplexFloatBlitzArray2DToBinaryFile(const string filename, blitz::Array<complex<float>, 2>& A);

#endif //READWRITEBLITZARRAYFROMFILE_H
