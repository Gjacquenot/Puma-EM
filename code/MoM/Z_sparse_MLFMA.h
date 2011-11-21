#ifndef Z_SPARSE_MLFMA_H
#define Z_SPARSE_MLFMA_H

#include <complex>

using namespace std;

#include "readWriteBlitzArrayFromFile.h"
 
class Z_sparse_MLFMA {
    int N_near;
    int N_test_RWG;
    int N_src_RWG;
    blitz::Array<int, 1> test_RWG_numbers;
    blitz::Array<int, 1> src_RWG_numbers;
    blitz::Array<int, 2> rowIndexToColumnIndexes;
    blitz::Array<std::complex<float>, 1> Z_CFIE_near;
  public:
    // constructors
    Z_sparse_MLFMA(void);
    Z_sparse_MLFMA(const blitz::Array<int, 1> /*test_RWG_numbers*/,
                   const blitz::Array<int, 1> /*src_RWG_numbers*/,
                   const blitz::Array<int, 2> /*rowIndexToColumnIndexes*/,
                   const blitz::Array<std::complex<float>, 1> /*Z_CFIE_near*/);
    // not really a constructor, but oh well....
    void setZ_sparse_MLFMAFromFile(const string /*path*/,
                                   const string /*Z_name*/,
                                   const int /*chunkNumber*/);
    void copyZ_sparse_MLFMA (const Z_sparse_MLFMA&);
    Z_sparse_MLFMA(const Z_sparse_MLFMA&); // copy constructor
    Z_sparse_MLFMA& operator=(const Z_sparse_MLFMA&); // copy assignment operator
    ~Z_sparse_MLFMA();

    // specific functions
    const int getN_near(void) const {return N_near;};
    void setN_near(const int N) {N_near = N;};
    const int getN_test_RWG(void) const {return N_test_RWG;};
    void setN_test_RWG(const int N) {N_test_RWG = N;};
    const int getN_src_RWG(void) const {return N_src_RWG;};
    void setN_src_RWG(const int N) {N_src_RWG = N;};

    const blitz::Array<int, 1> get_test_RWG_numbers(void) const {return test_RWG_numbers;};
    void set_test_RWG_numbers(const blitz::Array<int, 1>& V) {test_RWG_numbers.resize(V.size()); test_RWG_numbers = V;};
    const blitz::Array<int, 1> get_src_RWG_numbers(void) const {return src_RWG_numbers;};
    void set_src_RWG_numbers(const blitz::Array<int, 1>& V) {src_RWG_numbers.resize(V.size()); src_RWG_numbers = V;};
    const blitz::Array<int, 2> getRowIndexToColumnIndexes(void) const {return rowIndexToColumnIndexes;};
    void setRowIndexToColumnIndexes(const blitz::Array<int, 2>& V) {rowIndexToColumnIndexes.resize(V.size()); rowIndexToColumnIndexes = V;};
    const blitz::Array<std::complex<float>, 1> getZ_CFIE_near(void) const {return Z_CFIE_near;};
    void setZ_CFIE_near(const blitz::Array<std::complex<float>, 1>& V) {Z_CFIE_near.resize(V.size()); Z_CFIE_near = V;};
    void printZ_CFIE_near(void) {blitz::cout << "Z_CFIE_near = " << Z_CFIE_near << endl;}
    void matvec_Z_PQ_near(blitz::Array<std::complex<float>, 1>& /*ZI_PQ*/,
                          const blitz::Array<std::complex<float>, 1>& /*I_PQ*/);
};

Z_sparse_MLFMA::Z_sparse_MLFMA(void){};

Z_sparse_MLFMA::Z_sparse_MLFMA(const blitz::Array<int, 1> test_RWG_numbers,
                   const blitz::Array<int, 1> src_RWG_numbers,
                   const blitz::Array<int, 2> rowIndexToColumnIndexes,
                   const blitz::Array<std::complex<float>, 1> Z_CFIE_near) 
{
  set_test_RWG_numbers(test_RWG_numbers);
  set_src_RWG_numbers(src_RWG_numbers);
  setRowIndexToColumnIndexes(rowIndexToColumnIndexes);
  setZ_CFIE_near(Z_CFIE_near);
  setN_test_RWG(test_RWG_numbers.size());
  setN_near(Z_CFIE_near.size());
  setN_src_RWG(src_RWG_numbers.size());
};

void Z_sparse_MLFMA::copyZ_sparse_MLFMA(const Z_sparse_MLFMA& Z_sparse_MLFMAToCopy) /// copy member function
{
  N_near = Z_sparse_MLFMAToCopy.getN_near();
  N_test_RWG = Z_sparse_MLFMAToCopy.getN_test_RWG();
  N_src_RWG = Z_sparse_MLFMAToCopy.getN_src_RWG();
  set_test_RWG_numbers(Z_sparse_MLFMAToCopy.get_test_RWG_numbers());
  set_src_RWG_numbers(Z_sparse_MLFMAToCopy.get_src_RWG_numbers());
  setRowIndexToColumnIndexes(Z_sparse_MLFMAToCopy.getRowIndexToColumnIndexes());
  setZ_CFIE_near(Z_sparse_MLFMAToCopy.getZ_CFIE_near());
};

Z_sparse_MLFMA::Z_sparse_MLFMA(const Z_sparse_MLFMA& Z_sparse_MLFMAToCopy) /// copy constructor
{
  copyZ_sparse_MLFMA(Z_sparse_MLFMAToCopy);
};

Z_sparse_MLFMA& Z_sparse_MLFMA::operator=(const Z_sparse_MLFMA& Z_sparse_MLFMAToCopy) { /// copy assignment
  copyZ_sparse_MLFMA(Z_sparse_MLFMAToCopy);
  return *this;
};

Z_sparse_MLFMA::~Z_sparse_MLFMA() {
  test_RWG_numbers.free();
  src_RWG_numbers.free();
  rowIndexToColumnIndexes.free();
  Z_CFIE_near.free();
};

void Z_sparse_MLFMA::setZ_sparse_MLFMAFromFile(const string path, const string Z_name, const int chunkNumber)
{
  int N_RWG, N_near, N_q_array;
  const string chunkNumberString = intToString(chunkNumber);
  {
    string filename = path + "N_test_RWG" + chunkNumberString + ".txt";
    readIntFromASCIIFile(filename, N_RWG);
  }
  {
    string filename = path + "N_src_RWG" + chunkNumberString + ".txt";
    readIntFromASCIIFile(filename, N_q_array);
  }
  {
    string filename = path + "N_near" + chunkNumberString + ".txt";
    readIntFromASCIIFile(filename, N_near);
  }
  setN_near(N_near);
  setN_test_RWG(N_RWG);
  test_RWG_numbers.resize(N_RWG);
  rowIndexToColumnIndexes.resize(N_RWG, 2);
  src_RWG_numbers.resize(N_q_array);
  Z_CFIE_near.resize(N_near);
  // reading test_RWG_numbers
  {
    string filename = path + "test_RWG_numbers" + chunkNumberString + ".txt";
    readIntBlitzArray1DFromBinaryFile(filename, test_RWG_numbers);
  }
  // reading rowIndexToColumnIndexes
  {
    string filename = path + "rowIndexToColumnIndexes" + chunkNumberString + ".txt";
    readIntBlitzArray2DFromBinaryFile(filename, rowIndexToColumnIndexes);
  }
  // reading src_RWG_numbers
  {
    string filename = path + "src_RWG_numbers" + chunkNumberString + ".txt";
    readIntBlitzArray1DFromBinaryFile(filename, src_RWG_numbers);
  }
  // reading Z_CFIE_near
  {
    string filename = path + Z_name + chunkNumberString + ".txt";
    readComplexFloatBlitzArray1DFromBinaryFile(filename, Z_CFIE_near);
  }
};

void Z_sparse_MLFMA::matvec_Z_PQ_near(blitz::Array<std::complex<float>, 1>& ZI_PQ,
                                      const blitz::Array<std::complex<float>, 1>& I_PQ)
/**
 * matrix-vector multiplication for a sparse matrix data structure
 * such as the compressed row storage scheme
 */
{
  int indexInZ_CFIE = 0;
  for (int i=0 ; i<N_test_RWG ; i++) {
    const int test_RWG_number = test_RWG_numbers(i);
    const int startIndexInSrcRWG_numbers = rowIndexToColumnIndexes(i, 0);
    const int stopIndexInSrc_RWG_numbers = rowIndexToColumnIndexes(i, 1);
    for (int j=startIndexInSrcRWG_numbers ; j<stopIndexInSrc_RWG_numbers ; j++) {
      ZI_PQ(test_RWG_number) += Z_CFIE_near(indexInZ_CFIE) * I_PQ(src_RWG_numbers(j));
      indexInZ_CFIE++;
    }
  }
};

#endif
