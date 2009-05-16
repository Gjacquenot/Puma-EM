#include "readWriteBlitzArrayFromFile.h"

template <class T>
bool stringTo(T& t, 
              const std::string& s, 
              std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

double stringToDouble(const std::string& s)
{ 
  double result;
  stringTo(result, s, std::dec);
  return result;
}

template<typename T> 
std::string toString(const T& x)
{
  std::ostringstream oss;
  oss << x;
  return oss.str();
}
std::string intToString(const int & x)
{
  return toString(x);
}

// write scalar values
template <typename T>
void writeScalarToASCIIFile(const string filename, const T & x)
{
  ofstream ofs (filename.c_str());
  if (! ofs.is_open()) { 
    cout << "readWriteBlitzArrayFromFile.cpp::writeScalarToASCIIFile: error opening " << filename << endl; 
    exit (1);
  }
  ofs.precision(16);
  ofs << x << endl;
  ofs.close();
}
// Since there are problems with Scipy.weave when using templates, I define the following functions
void writeIntToASCIIFile(const string filename, const int & x) {
  writeScalarToASCIIFile(filename, x);
}
void writeFloatToASCIIFile(const string filename, const float & x) {
  writeScalarToASCIIFile(filename, x);
}
void writeDoubleToASCIIFile(const string filename, const double & x) {
  writeScalarToASCIIFile(filename, x);
}
void writeComplexFloatToASCIIFile(const string filename, const complex<float> & x) {
  writeScalarToASCIIFile(filename, x);
}
void writeComplexDoubleToASCIIFile(const string filename, const complex<double> & x) {
  writeScalarToASCIIFile(filename, x);
}


// read scalar values
template <typename T>
void readScalarFromASCIIFile(const string filename, T & x)
{
  ifstream ifs (filename.c_str());
  if (! ifs.is_open()) { 
    cout << "readWriteBlitzArrayFromFile.cpp::readScalarFromASCIIFile : error opening " << filename << endl; 
    exit (1);
  }
  else ifs >> x;
  ifs.close();
}
void readStringFromASCIIFile(const string filename, string & x)
{
  ifstream ifs (filename.c_str());
  if (! ifs.is_open()) { 
    cout << "readWriteBlitzArrayFromFile.cpp::readStringFromASCIIFile : error opening " << filename << endl; 
    exit (1);
  }
  else ifs >> x;
  ifs.close();
}
// Since there are problems with Scipy.weave when using templates, I define the following functions
void readIntFromASCIIFile(const string filename, int & x) {
  readScalarFromASCIIFile(filename, x);
}
void readFloatFromASCIIFile(const string filename, float & x) {
  readScalarFromASCIIFile(filename, x);
}
void readDoubleFromASCIIFile(const string filename, double & x) {
  readScalarFromASCIIFile(filename, x);
}
void readComplexFloatFromASCIIFile(const string filename, complex<float> & x) {
  readScalarFromASCIIFile(filename, x);
}
void readComplexDoubleFromASCIIFile(const string filename, complex<double> & x) {
  readScalarFromASCIIFile(filename, x);
}

/*********************************************
/*
/*               1D Arrays
/*
/*********************************************/
template<typename T>
void readBlitzArray1DFromASCIIFile(const string filename, blitz::Array<T, 1>& A)
{
  blitz::ifstream ifs(filename.c_str());
  if (! ifs.is_open()) { 
    cout << "readWriteBlitzArrayFromFile.cpp::readBlitzArray1DFromASCIIFile: error opening " << filename << endl;
    exit (1);
  }
  ifs >> A;
  ifs.close();
}

template<typename T>
void writeBlitzArray1DToASCIIFile(const string filename, const blitz::Array<T, 1>& A)
{
  blitz::ofstream ofs(filename.c_str());
  if (! ofs.is_open()) { 
    cout << "readWriteBlitzArrayFromFile.cpp::writeBlitzArray1DToASCIIFile: error opening " << filename << endl;
    exit (1);
  }
  ofs << A << endl;
  ofs.close();
}
template<typename T>
void writeBlitzArray1DToBinaryFile(const string filename, const blitz::Array<T, 1>& A, const int sizeOfItem)
{
  blitz::ofstream ofs(filename.c_str(), blitz::ios::binary);
  if (! ofs.is_open()) { 
    cout << "readWriteBlitzArrayFromFile.cpp::writeBlitzArray1DToBinaryFile: error opening " << filename << endl;
    exit (1);
  }
  ofs.write((char *)(A.data()), A.size()*sizeOfItem);
  ofs.close();
}

void readIntBlitzArray1DFromASCIIFile(const string filename, blitz::Array<int, 1>& A)
{
  readBlitzArray1DFromASCIIFile(filename, A);
}

void readFloatBlitzArray1DFromASCIIFile(const string filename, blitz::Array<float, 1>& A)
{
  readBlitzArray1DFromASCIIFile(filename, A);
}

void readDoubleBlitzArray1DFromASCIIFile(const string filename, blitz::Array<double, 1>& A)
{
  readBlitzArray1DFromASCIIFile(filename, A);
}

void readComplexFloatBlitzArray1DFromASCIIFile(const string filename, blitz::Array<std::complex<float>, 1>& A)
{
  readBlitzArray1DFromASCIIFile(filename, A);
}

void readComplexDoubleBlitzArray1DFromASCIIFile(const string filename, blitz::Array<std::complex<double>, 1>& A)
{
  readBlitzArray1DFromASCIIFile(filename, A);
}

void writeIntBlitzArray1DToASCIIFile(const string filename, blitz::Array<int, 1>& A)
{
  writeBlitzArray1DToASCIIFile(filename, A);
}
void writeFloatBlitzArray1DToASCIIFile(const string filename, const blitz::Array<float, 1>& A)
{
  writeBlitzArray1DToASCIIFile(filename, A);
}
void writeComplexFloatBlitzArray1DToASCIIFile(const string filename, const blitz::Array<complex<float>, 1>& A)
{
  writeBlitzArray1DToASCIIFile(filename, A);
}
void writeComplexFloatBlitzArray1DToBinaryFile(const string filename, const blitz::Array<complex<float>, 1>& A)
{
  const int itemsize = 8;
  writeBlitzArray1DToBinaryFile(filename, A, itemsize);
}
void writeDoubleBlitzArray1DToBinaryFile(const string filename, blitz::Array<double, 1>& A) {
  const int itemsize = 8;
  writeBlitzArray1DToBinaryFile(filename, A, itemsize);
}
void writeIntBlitzArray1DToBinaryFile(const string filename, blitz::Array<int, 1>& A) {
  const int itemsize = 4;
  writeBlitzArray1DToBinaryFile(filename, A, itemsize);
}



////////////////////////////////////////////////////////////////////////////////
// read 1-D binary blitz Arrays
template<typename T>
void readBlitzArray1DFromBinaryFile(const string filename, blitz::Array<T, 1>& A, const int itemsize)
{
  blitz::ifstream ifs(filename.c_str(), blitz::ios::binary);
  if (! ifs.is_open()) { 
    cout << "readWriteBlitzArrayFromFile.cpp::readBlitzArray1DFromBinaryFile: error opening " << filename << endl;
    exit (1);
  }
  ifs.read((char *)(A.data()), A.size()*itemsize);
  ifs.close();
}

void readIntBlitzArray1DFromBinaryFile(const string filename, blitz::Array<int, 1>& A)
{
  const int itemsize = 4;
  readBlitzArray1DFromBinaryFile(filename, A, itemsize);
}
void readFloatBlitzArray1DFromBinaryFile(const string filename, blitz::Array<float, 1>& A)
{
  const int itemsize = 4;
  readBlitzArray1DFromBinaryFile(filename, A, itemsize);
}
void readDoubleBlitzArray1DFromBinaryFile(const string filename, blitz::Array<double, 1>& A)
{
  const int itemsize = 8;
  readBlitzArray1DFromBinaryFile(filename, A, itemsize);
}
void readComplexFloatBlitzArray1DFromBinaryFile(const string filename, blitz::Array<std::complex<float>, 1>& A)
{
  const int itemsize = 8;
  readBlitzArray1DFromBinaryFile(filename, A, itemsize);
}
void readComplexDoubleBlitzArray1DFromBinaryFile(const string filename, blitz::Array<std::complex<double>, 1>& A)
{
  const int itemsize = 16;
  readBlitzArray1DFromBinaryFile(filename, A, itemsize);
}

/*********************************************
/*
/*               2D Arrays
/*
/*********************************************/

// read and write 2-D ASCII blitz Arrays: generic functions
template<typename T>
void readBlitzArray2DFromASCIIFile(const string filename, blitz::Array<T, 2>& A)
{
  blitz::ifstream ifs(filename.c_str());
  if (! ifs.is_open()) { 
    cout << "readWriteBlitzArrayFromFile.cpp::readBlitzArray2DFromASCIIFile: error opening " << filename << endl;
    exit (1);
  }
  ifs >> A;
  ifs.close();
}
template<typename T>
void readBlitzArray2DFromBinaryFile(const string filename, blitz::Array<T, 2>& A, const int itemsize)
{
  blitz::ifstream ifs(filename.c_str(), blitz::ios::binary);
  if (! ifs.is_open()) {
    cout << "readWriteBlitzArrayFromFile.cpp::readBlitzArray2DFromBinaryFile: error opening " << filename << endl;
    exit (1);
  }
  ifs.read((char *)(A.data()), A.size()*itemsize);
  ifs.close();
}
template<typename T>
void writeBlitzArray2DToASCIIFile(const string filename, blitz::Array<T, 2>& A)
{
  blitz::ofstream ofs(filename.c_str());
  if (! ofs.is_open()) { 
    cout << "readWriteBlitzArrayFromFile.cpp::writeBlitzArray2DToASCIIFile: error opening " << filename << endl;
    exit (1);
  }
  ofs << A << endl;
  ofs.close();
}

template<typename T>
void writeBlitzArray2DToBinaryFile(const string filename, blitz::Array<T, 2>& A, const int sizeOfItem)
{
  blitz::ofstream ofs(filename.c_str(), blitz::ios::binary);
  if (! ofs.is_open()) { 
    cout << "readWriteBlitzArrayFromFile.cpp::writeBlitzArray2DToBinaryFile: error opening " << filename << endl;
    exit (1);
  }
  ofs.write((char *)(A.data()), A.size()*sizeOfItem);
  ofs.close();
}

// read 2D arrays from ASCII
void readIntBlitzArray2DFromASCIIFile(const string filename, blitz::Array<int, 2>& A)
{
  readBlitzArray2DFromASCIIFile(filename, A);
}

void readFloatBlitzArray2DFromASCIIFile(const string filename, blitz::Array<float, 2>& A)
{
  readBlitzArray2DFromASCIIFile(filename, A);
}

void readDoubleBlitzArray2DFromASCIIFile(const string filename, blitz::Array<double, 2>& A)
{
  readBlitzArray2DFromASCIIFile(filename, A);
}

void readComplexFloatBlitzArray2DFromASCIIFile(const string filename, blitz::Array<std::complex<float>, 2>& A)
{
  readBlitzArray2DFromASCIIFile(filename, A);
}

void readComplexDoubleBlitzArray2DFromACSIIFile(const string filename, blitz::Array<std::complex<double>, 2>& A)
{
  readBlitzArray2DFromASCIIFile(filename, A);
}

// read 2-D binary blitz Arrays
void readIntBlitzArray2DFromBinaryFile(const string filename, blitz::Array<int, 2>& A)
{
  const int itemsize = 4;
  readBlitzArray2DFromBinaryFile(filename, A, itemsize);
}

void readFloatBlitzArray2DFromBinaryFile(const string filename, blitz::Array<float, 2>& A)
{
  const int itemsize = 4;
  readBlitzArray2DFromBinaryFile(filename, A, itemsize);
}

void readDoubleBlitzArray2DFromBinaryFile(const string filename, blitz::Array<double, 2>& A)
{
  const int itemsize = 8;
  readBlitzArray2DFromBinaryFile(filename, A, itemsize);
}

void readComplexFloatBlitzArray2DFromBinaryFile(const string filename, blitz::Array<std::complex<float>, 2>& A)
{
  const int itemsize = 8;
  readBlitzArray2DFromBinaryFile(filename, A, itemsize);
}

void readComplexDoubleBlitzArray2DFromBinaryFile(const string filename, blitz::Array<std::complex<double>, 2>& A)
{
  const int itemsize = 16;
  readBlitzArray2DFromBinaryFile(filename, A, itemsize);
}

// writing 2D arrays: ASCII
void writeIntBlitzArray2DToASCIIFile(const string filename, blitz::Array<int, 2>& A) {
  writeBlitzArray2DToASCIIFile(filename, A);
}
void writeComplexFloatBlitzArray2DToASCIIFile(const string filename, blitz::Array<std::complex<float>, 2>& A) {
  writeBlitzArray2DToASCIIFile(filename, A);
}
void writeComplexDoubleBlitzArray2DToASCIIFile(const string filename, blitz::Array<std::complex<double>, 2>& A) {
  writeBlitzArray2DToASCIIFile(filename, A);
}
void writeDoubleBlitzArray2DToASCIIFile(const string filename, blitz::Array<double, 2>& A) {
  writeBlitzArray2DToASCIIFile(filename, A);
}
void writeFloatBlitzArray2DToASCIIFile(const string filename, blitz::Array<float, 2>& A) {
  writeBlitzArray2DToASCIIFile(filename, A);
}
// writing 2D arrays: BINARY
void writeIntBlitzArray2DToBinaryFile(const string filename, blitz::Array<int, 2>& A) {
  const int itemsize = 4;
  writeBlitzArray2DToBinaryFile(filename, A, itemsize);
}
void writeDoubleBlitzArray2DToBinaryFile(const string filename, blitz::Array<double, 2>& A) {
  const int itemsize = 8;
  writeBlitzArray2DToBinaryFile(filename, A, itemsize);
}
void writeComplexFloatBlitzArray2DToBinaryFile(const string filename, blitz::Array<complex<float>, 2>& A) {
  const int itemsize = 8;
  writeBlitzArray2DToBinaryFile(filename, A, itemsize);
}


