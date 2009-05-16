#ifndef Z_SPARSE_MLFMA_H
#define Z_SPARSE_MLFMA_H

#include <complex>

using namespace std;

#include "readWriteBlitzArrayFromFile.h"
 
class Z_sparse_MLFMA {
    int N_near;
    int N_RWG;
    int N_q_array;
    blitz::Array<int, 1> RWG_numbers;
    blitz::Array<int, 1> q_array;
    blitz::Array<int, 2> rowIndexToColumnIndexes;
    blitz::Array<std::complex<float>, 1> Z_CFIE_near;
  public:
    // constructors
    Z_sparse_MLFMA(void);
    Z_sparse_MLFMA(const blitz::Array<int, 1> /*RWG_numbers*/,
                   const blitz::Array<int, 1> /*q_array*/,
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
    const int getN_RWG(void) const {return N_RWG;};
    void setN_RWG(const int N) {N_RWG = N;};
    const int getN_q_array(void) const {return N_q_array;};
    void setN_q_Array(const int N) {N_q_array = N;};

    const blitz::Array<int, 1> getRWG_numbers(void) const {return RWG_numbers;};
    void setRWG_numbers(const blitz::Array<int, 1>& V) {RWG_numbers.resize(V.size()); RWG_numbers = V;};
    const blitz::Array<int, 1> getQ_array(void) const {return q_array;};
    void setQ_array(const blitz::Array<int, 1>& V) {q_array.resize(V.size()); q_array = V;};
    const blitz::Array<int, 2> getRowIndexToColumnIndexes(void) const {return rowIndexToColumnIndexes;};
    void setRowIndexToColumnIndexes(const blitz::Array<int, 2>& V) {rowIndexToColumnIndexes.resize(V.size()); rowIndexToColumnIndexes = V;};
    const blitz::Array<std::complex<float>, 1> getZ_CFIE_near(void) const {return Z_CFIE_near;};
    void setZ_CFIE_near(const blitz::Array<std::complex<float>, 1>& V) {Z_CFIE_near.resize(V.size()); Z_CFIE_near = V;};
    void printZ_CFIE_near(void) {blitz::cout << "Z_CFIE_near = " << Z_CFIE_near << endl;}
    void matvec_Z_PQ_near(blitz::Array<std::complex<float>, 1>& /*ZI_PQ*/,
                          const blitz::Array<std::complex<float>, 1>& /*I_PQ*/);
};

Z_sparse_MLFMA::Z_sparse_MLFMA(void){};

Z_sparse_MLFMA::Z_sparse_MLFMA(const blitz::Array<int, 1> RWG_numbers,
                   const blitz::Array<int, 1> q_array,
                   const blitz::Array<int, 2> rowIndexToColumnIndexes,
                   const blitz::Array<std::complex<float>, 1> Z_CFIE_near) 
{
  setRWG_numbers(RWG_numbers);
  setQ_array(q_array);
  setRowIndexToColumnIndexes(rowIndexToColumnIndexes);
  setZ_CFIE_near(Z_CFIE_near);
  setN_RWG(RWG_numbers.size());
  setN_near(Z_CFIE_near.size());
  setN_q_Array(q_array.size());
};

void Z_sparse_MLFMA::copyZ_sparse_MLFMA(const Z_sparse_MLFMA& Z_sparse_MLFMAToCopy) /// copy member function
{
  N_near = Z_sparse_MLFMAToCopy.getN_near();
  N_RWG = Z_sparse_MLFMAToCopy.getN_RWG();
  N_q_array = Z_sparse_MLFMAToCopy.getN_q_array();
  setRWG_numbers(Z_sparse_MLFMAToCopy.getRWG_numbers());
  setQ_array(Z_sparse_MLFMAToCopy.getQ_array());
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
  RWG_numbers.free();
  q_array.free();
  rowIndexToColumnIndexes.free();
  Z_CFIE_near.free();
};

void Z_sparse_MLFMA::setZ_sparse_MLFMAFromFile(const string path, const string Z_name, const int chunkNumber)
{
  int N_RWG, N_near, N_q_array;
  const string chunkNumberString = intToString(chunkNumber);
  {
    string filename = path + "N_RWG" + chunkNumberString + ".txt";
    readIntFromASCIIFile(filename, N_RWG);
  }
  {
    string filename = path + "N_q_array" + chunkNumberString + ".txt";
    readIntFromASCIIFile(filename, N_q_array);
  }
  {
    string filename = path + "N_near" + chunkNumberString + ".txt";
    readIntFromASCIIFile(filename, N_near);
  }
  setN_near(N_near);
  setN_RWG(N_RWG);
  RWG_numbers.resize(N_RWG);
  rowIndexToColumnIndexes.resize(N_RWG, 2);
  q_array.resize(N_q_array);
  Z_CFIE_near.resize(N_near);
  // reading RWG_numbers
  {
    string filename = path + "RWG_numbers" + chunkNumberString + ".txt";
    readIntBlitzArray1DFromBinaryFile(filename, RWG_numbers);
  }
  // reading rowIndexToColumnIndexes
  {
    string filename = path + "rowIndexToColumnIndexes" + chunkNumberString + ".txt";
    readIntBlitzArray2DFromBinaryFile(filename, rowIndexToColumnIndexes);
  }
  // reading q_array
  {
    string filename = path + "q_array" + chunkNumberString + ".txt";
    readIntBlitzArray1DFromBinaryFile(filename, q_array);
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
  for (int i=0 ; i<N_RWG ; i++) {
    const int RWG_number = RWG_numbers(i);
    const int startIndexInQArray = rowIndexToColumnIndexes(i, 0);
    const int stopIndexInQArray = rowIndexToColumnIndexes(i, 1);
    for (int j=startIndexInQArray ; j<stopIndexInQArray ; j++) {
      ZI_PQ(RWG_number) += Z_CFIE_near(indexInZ_CFIE) * I_PQ(q_array(j));
      indexInZ_CFIE++;
    }
  }
};

#endif
