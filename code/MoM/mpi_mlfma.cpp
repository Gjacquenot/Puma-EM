#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <blitz/array.h>
#include <vector>
#include <algorithm>
#include <mpi.h>

using namespace std;

#include "GetMemUsage.h"
#include "readWriteBlitzArrayFromFile.h"
#include "Z_sparse_MLFMA.h"
#include "octtree.h"
#include "./iterative/iterative.h"
#include "V_E_V_H.h"
#include "mesh.h"
#include "EMConstants.h"

/****************************************************************************/
/******************************* MLFMA matvec *******************************/
/****************************************************************************/
void gatherAndRedistribute(blitz::Array<std::complex<float>, 1>& x, const int my_id, const int num_procs)
{
  int ierror, N_RWG = x.size();
  blitz::Array<std::complex<float>, 1> recvBuf;
  if (my_id==0) recvBuf.resize(N_RWG);
  ierror = MPI_Reduce(x.data(), recvBuf.data(), N_RWG, MPI_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
  if (my_id==0) x = recvBuf;
  ierror = MPI_Bcast(x.data(), N_RWG, MPI_COMPLEX, 0, MPI_COMM_WORLD);
}

class MatvecMLFMA {

  public:
    Octtree * pOcttree;
    int procNumber;
    int totalProcNumber;
    int N_RWG;
    blitz::Array<int, 1> localRWGnumbers;
    string simuDir;

    // constructors
    MatvecMLFMA(void){;};
    MatvecMLFMA(Octtree & /*octtree*/,
                const int /*numberOfRWG*/,
                const blitz::Array<int, 1>& /*localRWGindexes*/,
		const string & /*simuDir*/);
    // destructor
    ~MatvecMLFMA(void){localRWGnumbers.free();};
    // copy operators
    void copyMatvecMLFMA (const MatvecMLFMA&);
    MatvecMLFMA(const MatvecMLFMA&); // copy constructor
    MatvecMLFMA& operator=(const MatvecMLFMA&); // copy assignment operator

    // function
    void setProcNumber(const int n) {procNumber = n;};
    const int getProcNumber(void) const {return procNumber;};
    void setTotalProcNumber(const int n) {totalProcNumber = n;};
    const int getTotalProcNumber(void) const {return totalProcNumber;};
    const int getN_RWG(void) const {return N_RWG;};
    void matvecZnear(blitz::Array<std::complex<float>, 1> & /*y*/,
                     const blitz::Array<std::complex<float>, 1> & /*x*/);
    blitz::Array<std::complex<float>, 1> matvec(const blitz::Array<std::complex<float>, 1> & /*x*/);
};

MatvecMLFMA::MatvecMLFMA(Octtree & octtree,
                         const int numberOfRWG,
                         const blitz::Array<int, 1>& localRWGindexes,
                         const string & simu_dir)

{
  pOcttree = &octtree;
  procNumber = MPI::COMM_WORLD.Get_rank();
  totalProcNumber = MPI::COMM_WORLD.Get_size();
  N_RWG = numberOfRWG;
  localRWGnumbers.resize(localRWGindexes.size());
  localRWGnumbers = localRWGindexes;
  simuDir = simu_dir;
}

void MatvecMLFMA::copyMatvecMLFMA(const MatvecMLFMA& matvecMLFMAtoCopy) // copy member function
{
  pOcttree = matvecMLFMAtoCopy.pOcttree;
  procNumber = matvecMLFMAtoCopy.procNumber;
  totalProcNumber = matvecMLFMAtoCopy.totalProcNumber;
  N_RWG = matvecMLFMAtoCopy.N_RWG;
  localRWGnumbers.resize(matvecMLFMAtoCopy.localRWGnumbers.size());
  localRWGnumbers = matvecMLFMAtoCopy.localRWGnumbers;
  simuDir = matvecMLFMAtoCopy.simuDir;
}

MatvecMLFMA::MatvecMLFMA(const MatvecMLFMA& matvecMLFMAtoCopy) // copy constructor
{
  copyMatvecMLFMA(matvecMLFMAtoCopy);
}

MatvecMLFMA& MatvecMLFMA::operator=(const MatvecMLFMA& matvecMLFMAtoCopy) { // copy assignment
  copyMatvecMLFMA(matvecMLFMAtoCopy);
  return *this;
}

void MatvecMLFMA::matvecZnear(blitz::Array<std::complex<float>, 1> & y, const blitz::Array<std::complex<float>, 1> & x)
{
  const int my_id = MPI::COMM_WORLD.Get_rank();
  const string pathToReadFrom = simuDir + "/tmp" + intToString(my_id) + "/Z_near/", Z_name = "Z_CFIE_near";
  Array<int, 1> chunkNumbers;
  readIntBlitzArray1DFromASCIIFile(pathToReadFrom + "chunkNumbers.txt", chunkNumbers);
  Z_sparse_MLFMA Z_near;
  for (int i=0 ; i<chunkNumbers.size() ; i++) {
    int number = chunkNumbers(i);
    Z_near.setZ_sparse_MLFMAFromFile(pathToReadFrom, Z_name, number);
    Z_near.matvec_Z_PQ_near(y, x);
  }
}

blitz::Array<std::complex<float>, 1> MatvecMLFMA::matvec(const blitz::Array<std::complex<float>, 1> & x)
{
  // creation of a local solution vector for MLFMA
  blitz::Array<std::complex<float>, 1> y_local_MLFMA(this->localRWGnumbers.size());
  y_local_MLFMA = 0.0;

  // far-field multiplication
  pOcttree->ZIFarComputation(y_local_MLFMA, x);
  
  // distribution of x among processes
  blitz::Array<std::complex<float>, 1> x_global(this->N_RWG);
  x_global = 0.0;
  for (int i=0 ; i<this->localRWGnumbers.size() ; ++i) x_global(this->localRWGnumbers(i)) = x(i);
  // we should now gather and redistribute the result among the processes
  gatherAndRedistribute(x_global, getProcNumber(), getTotalProcNumber());

  // near-field multiplication
  const int my_id = MPI::COMM_WORLD.Get_rank();
  const string pathToReadFrom = simuDir + "/tmp" + intToString(my_id) + "/Z_near/";
  int local_N_src_RWG;
  string filename = pathToReadFrom + "local_N_src_RWG.txt";
  readIntFromASCIIFile(filename, local_N_src_RWG);
  blitz::Array<int, 1> local_src_RWG_numbers(local_N_src_RWG);
  filename = pathToReadFrom + "local_src_RWG_numbers.txt";
  readIntBlitzArray1DFromBinaryFile(filename, local_src_RWG_numbers);
  blitz::Array<std::complex<float>, 1> x_local_Z(local_N_src_RWG);
  for (int i=0; i<local_N_src_RWG; i++) x_local_Z(i) = x_global(local_src_RWG_numbers(i));
  x_global.free();
  
  int local_N_test_RWG;
  filename = pathToReadFrom + "local_N_test_RWG.txt";
  readIntFromASCIIFile(filename, local_N_test_RWG);
  blitz::Array<int, 1> local_test_RWG_numbers(local_N_test_RWG);
  filename = pathToReadFrom + "local_test_RWG_numbers.txt";
  readIntBlitzArray1DFromBinaryFile(filename, local_test_RWG_numbers);
  blitz::Array<std::complex<float>, 1> y_local_Z(local_N_test_RWG);
  y_local_Z = 0.0;
  matvecZnear(y_local_Z, x_local_Z);
  x_local_Z.free();
  // we should now gather and redistribute the result among the processes
  // but first we add the results from the multiplications
  blitz::Array<std::complex<float>, 1> y_global(this->N_RWG);
  y_global = 0.0;
  for (int i=0; i<local_N_test_RWG; i++) y_global(local_test_RWG_numbers(i)) += y_local_Z(i);
  y_local_Z.free();
  for (int i=0 ; i<this->localRWGnumbers.size() ; ++i) y_global(this->localRWGnumbers(i)) += y_local_MLFMA(i);
  y_local_MLFMA.free();
  gatherAndRedistribute(y_global, getProcNumber(), getTotalProcNumber());
  // we now select only the elements to return
  blitz::Array<std::complex<float>, 1> yTmp(this->localRWGnumbers.size());
  for (int i=0 ; i<this->localRWGnumbers.size() ; ++i) yTmp(i) = y_global(this->localRWGnumbers(i));
  return yTmp;
}

/****************************************************************************/
/******************************* Left Frob precond **************************/
/****************************************************************************/
class LeftFrobPsolveMLFMA {

  public:
    int procNumber;
    int totalProcNumber;
    int N_RWG;
    blitz::Array<int, 1> localRWGnumbers;
    string simuDir;

    // constructors
    LeftFrobPsolveMLFMA(void){;};
    LeftFrobPsolveMLFMA(const int /*numberOfRWG*/,
                        const blitz::Array<int, 1>& /*localRWGindexes*/,
                        const string & /*simuDir*/);
    // destructor
    ~LeftFrobPsolveMLFMA(void){localRWGnumbers.free();};
    // copy operators
    void copyLeftFrobPsolveMLFMA (const LeftFrobPsolveMLFMA&);
    LeftFrobPsolveMLFMA(const LeftFrobPsolveMLFMA&); // copy constructor
    LeftFrobPsolveMLFMA& operator=(const LeftFrobPsolveMLFMA&); // copy assignment operator

    // function
    blitz::Array<std::complex<float>, 1> psolve(const blitz::Array<std::complex<float>, 1> & /*x*/);
};

LeftFrobPsolveMLFMA::LeftFrobPsolveMLFMA(const int numberOfRWG,
                                         const blitz::Array<int, 1>& localRWGindexes,
                                         const string & simu_dir)
{
  procNumber = MPI::COMM_WORLD.Get_rank();
  totalProcNumber = MPI::COMM_WORLD.Get_size();
  N_RWG = numberOfRWG;
  localRWGnumbers.resize(localRWGindexes.size());
  localRWGnumbers = localRWGindexes;
  simuDir = simu_dir;
}

void LeftFrobPsolveMLFMA::copyLeftFrobPsolveMLFMA(const LeftFrobPsolveMLFMA& leftFrobPsolveMLFMAtoCopy) // copy member function
{
  procNumber = leftFrobPsolveMLFMAtoCopy.procNumber;
  totalProcNumber = leftFrobPsolveMLFMAtoCopy.totalProcNumber;
  N_RWG = leftFrobPsolveMLFMAtoCopy.N_RWG;
  localRWGnumbers.resize(leftFrobPsolveMLFMAtoCopy.localRWGnumbers.size());
  localRWGnumbers = leftFrobPsolveMLFMAtoCopy.localRWGnumbers;
  simuDir = leftFrobPsolveMLFMAtoCopy.simuDir;
}

LeftFrobPsolveMLFMA::LeftFrobPsolveMLFMA(const LeftFrobPsolveMLFMA& leftFrobPsolveMLFMAtoCopy) // copy constructor
{
  copyLeftFrobPsolveMLFMA(leftFrobPsolveMLFMAtoCopy);
}

LeftFrobPsolveMLFMA& LeftFrobPsolveMLFMA::operator=(const LeftFrobPsolveMLFMA& leftFrobPsolveMLFMAtoCopy) { // copy assignment
  copyLeftFrobPsolveMLFMA(leftFrobPsolveMLFMAtoCopy);
  return *this;
}

blitz::Array<std::complex<float>, 1> LeftFrobPsolveMLFMA::psolve(const blitz::Array<std::complex<float>, 1> & x)
{
  const int my_id = MPI::COMM_WORLD.Get_rank();
  blitz::Array<std::complex<float>, 1> x_global(this->N_RWG);
  x_global = 0.0;
  for (int i=0 ; i<localRWGnumbers.size() ; ++i) x_global(localRWGnumbers(i)) = x(i);
  gatherAndRedistribute(x_global, procNumber, totalProcNumber);

  const string pathToReadFrom = simuDir + "/tmp" + intToString(my_id) + "/Mg_LeftFrob/", Z_name = "Mg_LeftFrob";
  int local_N_src_RWG;
  string filename = pathToReadFrom + "local_N_src_RWG.txt";
  readIntFromASCIIFile(filename, local_N_src_RWG);
  blitz::Array<int, 1> local_src_RWG_numbers(local_N_src_RWG);
  filename = pathToReadFrom + "local_src_RWG_numbers.txt";
  readIntBlitzArray1DFromBinaryFile(filename, local_src_RWG_numbers);
  blitz::Array<std::complex<float>, 1> x_local_Y(local_N_src_RWG);
  for (int i=0; i<local_N_src_RWG; i++) x_local_Y(i) = x_global(local_src_RWG_numbers(i));
  x_global.free();  
  
  int local_N_test_RWG;
  filename = pathToReadFrom + "local_N_test_RWG.txt";
  readIntFromASCIIFile(filename, local_N_test_RWG);
  blitz::Array<int, 1> local_test_RWG_numbers(local_N_test_RWG);
  filename = pathToReadFrom + "local_test_RWG_numbers.txt";
  readIntBlitzArray1DFromBinaryFile(filename, local_test_RWG_numbers);
  blitz::Array<std::complex<float>, 1> y_local_Y(local_N_test_RWG);
  y_local_Y = 0.0;
  blitz::Array<int, 1> chunkNumbers;
  readIntBlitzArray1DFromASCIIFile(pathToReadFrom + "chunkNumbers.txt", chunkNumbers);
  Z_sparse_MLFMA Mg_LeftFrob;
  for (int i=0 ; i<chunkNumbers.size() ; i++) {
    int number = chunkNumbers(i);
    Mg_LeftFrob.setZ_sparse_MLFMAFromFile(pathToReadFrom, Z_name, number);
    //Mg_LeftFrob.printZ_CFIE_near();
    Mg_LeftFrob.matvec_Z_PQ_near(y_local_Y, x_local_Y);
  }
  x_local_Y.free();
  // we should now gather and redistribute the result among the processes
  blitz::Array<std::complex<float>, 1> y_global(this->N_RWG);
  y_global = 0.0;
  for (int i=0; i<local_N_test_RWG; i++) y_global(local_test_RWG_numbers(i)) += y_local_Y(i);
  y_local_Y.free();
  gatherAndRedistribute(y_global, procNumber, totalProcNumber);
  // we now select only the elements to return
  blitz::Array<std::complex<float>, 1> yTmp(localRWGnumbers.size());
  for (int i=0 ; i<localRWGnumbers.size() ; ++i) yTmp(i) = y_global(localRWGnumbers(i));
  return yTmp;
}

/****************************************************************************/
/*************************** AMLFMA Preconditioner **************************/
/****************************************************************************/
class PsolveAMLFMA { // should use FGMRES and _NOT_ GMRES or BiCGSTAB
    MatvecMLFMA *pMatvecMLFMA;
    LeftFrobPsolveMLFMA *pLeftFrobPsolveMLFMA;
    int inner_maxiter;
    int inner_restart;
    int my_id;
    int num_procs;
    int N_RWG;
    int NlocalRWG;
    float inner_tol;
    string inner_solver;
    string simuDir;

  public:
    // constructors
    PsolveAMLFMA(MatvecMLFMA & /*matvecMLFMA*/,
                 LeftFrobPsolveMLFMA & /*leftFrobPsolveMLFMA*/,
                 const float /*INNER_TOL*/,
                 const int /*INNER_MAXITER*/,
                 const int /*INNER_RESTART*/,
                 const int /*numberOfRWG*/,
                 const string & /*INNER_SOLVER*/,
                 const string & /*simuDir*/);
    ~PsolveAMLFMA(void){};
    // function
    blitz::Array<std::complex<float>, 1> psolve(const blitz::Array<std::complex<float>, 1> & /*x*/);
};

PsolveAMLFMA::PsolveAMLFMA (MatvecMLFMA & matvecMLFMA,
                            LeftFrobPsolveMLFMA & leftFrobPsolveMLFMA,
                            const float INNER_TOL,
                            const int INNER_MAXITER,
                            const int INNER_RESTART,
                            const int numberOfRWG,
                            const string & INNER_SOLVER,
                            const string & simu_dir)
{
  pMatvecMLFMA = &matvecMLFMA;
  pLeftFrobPsolveMLFMA = &leftFrobPsolveMLFMA;
  inner_maxiter = INNER_MAXITER;
  inner_restart = INNER_RESTART;
  inner_tol = INNER_TOL;
  my_id = MPI::COMM_WORLD.Get_rank();
  num_procs = MPI::COMM_WORLD.Get_size();
  N_RWG = numberOfRWG;
  NlocalRWG = matvecMLFMA.localRWGnumbers.size();
  inner_solver = INNER_SOLVER;
  simuDir = simu_dir;
}

blitz::Array<std::complex<float>, 1> PsolveAMLFMA::psolve(const blitz::Array<std::complex<float>, 1>& x)
{
  MatvecFunctor< std::complex<float>, MatvecMLFMA > matvec(pMatvecMLFMA, &MatvecMLFMA::matvec);
  PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > innerPsolve(pLeftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);

  const string TMP = simuDir + "/tmp" + intToString(my_id), ITERATIVE_DATA_PATH = TMP + "/iterative_data/";
  int iter, flag;
  double error;
  blitz::Array<std::complex<float>, 1> y(this->NlocalRWG);
  y = 0.0;
  if (inner_solver == "BICGSTAB") bicgstab(y, error, iter, flag, matvec, innerPsolve, x, inner_tol, 100, my_id, num_procs, ITERATIVE_DATA_PATH + "inner_convergence.txt");
  else if (inner_solver == "GMRES") gmres(y, error, iter, flag, matvec, innerPsolve, x, inner_tol, inner_restart, inner_maxiter, my_id, num_procs, ITERATIVE_DATA_PATH + "inner_convergence.txt");
  else if (inner_solver == "RGMRES") fgmres(y, error, iter, flag, matvec, innerPsolve, x, inner_tol, inner_restart, inner_maxiter, my_id, num_procs, ITERATIVE_DATA_PATH + "inner_convergence.txt");
  else if (inner_solver == "FGMRES") fgmres(y, error, iter, flag, matvec, innerPsolve, x, inner_tol, inner_restart, inner_maxiter, my_id, num_procs, ITERATIVE_DATA_PATH + "inner_convergence.txt");
  return y;
}

/****************************************************************************/
/************************* some useful functions ****************************/
/****************************************************************************/

void computeE_obs(blitz::Array<std::complex<double>, 2>& E_obs,
                  const blitz::Array<double, 2>& r_obs,
                  const LocalMesh &  local_target_mesh,
                  const blitz::Array<std::complex<float>, 1>& ZI,
                  const std::complex<float> eps_r,
                  const std::complex<float> mu_r,
                  const float w)
{
  int ierror, FULL_PRECISION = 1;
  blitz::Range all = blitz::Range::all();
  blitz::Array<std::complex<float>, 1> CFIE_for_tE(4);
  CFIE_for_tE = 1.0, 0.0, 0.0, 0.0;
  E_obs.resize(r_obs.extent(0), r_obs.extent(1));
  for (int j=0 ; j<r_obs.extent(0) ; ++j) {
    for (int i=0 ; i<3 ; ++i) {
      blitz::Array<std::complex<float>, 1> V_tE;
      blitz::Array<std::complex<double>, 1> J_obs(3);
      J_obs = 0.0;
      J_obs(i) = 1.0;
      local_V_CFIE_dipole (V_tE, J_obs, r_obs(j, all), local_target_mesh, w, eps_r, mu_r, CFIE_for_tE, FULL_PRECISION);
      std::complex<float> local_e_tmp = sum(V_tE * ZI), e_tmp;
      ierror = MPI_Allreduce(&local_e_tmp, &e_tmp, 1, MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
      E_obs(j, i) = e_tmp;
    }
  }
}

void computeForOneExcitation(Octtree & octtree,
                             LocalMesh & local_target_mesh,
                             const string SOLVER,
                             const string SIMU_DIR,
                             const string TMP,
                             const string OCTTREE_DATA_PATH,
                             const string MESH_DATA_PATH,
                             const string V_CFIE_DATA_PATH,
                             const string RESULT_DATA_PATH,
                             const string ITERATIVE_DATA_PATH)
{
  blitz::Range all = blitz::Range::all();
  int num_procs = MPI::COMM_WORLD.Get_size(), my_id = MPI::COMM_WORLD.Get_rank();
  const int master = 0, N_local_RWG = local_target_mesh.N_local_RWG;
  const float w = octtree.w;
  const std::complex<float> eps_r = octtree.eps_r, mu_r = octtree.mu_r;
  blitz::Array<int, 1> localRWGNumbers(local_target_mesh.localRWGNumbers.size());
  localRWGNumbers = local_target_mesh.localRWGNumbers;
  // functors declarations
  const string PRECOND = "FROB";
  int N_RWG, iter, flag, RESTART, MAXITER;
  double error, TOL;
  readIntFromASCIIFile(OCTTREE_DATA_PATH + "N_RWG.txt", N_RWG);
  MatvecMLFMA matvecMLFMA(octtree, N_RWG, localRWGNumbers, SIMU_DIR);
  MatvecFunctor< std::complex<float>, MatvecMLFMA > matvec(&matvecMLFMA, &MatvecMLFMA::matvec);
  LeftFrobPsolveMLFMA leftFrobPsolveMLFMA(N_RWG, localRWGNumbers, SIMU_DIR);

  // excitation : V_CFIE computation
  blitz::Array<std::complex<float>, 1> V_CFIE(N_local_RWG);
  V_CFIE = 0.0;
  int DIPOLES_EXCITATION, PLANE_WAVE_EXCITATION, V_FULL_PRECISION;
  readIntFromASCIIFile(V_CFIE_DATA_PATH + "DIPOLES_EXCITATION.txt", DIPOLES_EXCITATION);
  readIntFromASCIIFile(V_CFIE_DATA_PATH + "PLANE_WAVE_EXCITATION.txt", PLANE_WAVE_EXCITATION);
  readIntFromASCIIFile(V_CFIE_DATA_PATH + "V_FULL_PRECISION.txt", V_FULL_PRECISION);
  if (DIPOLES_EXCITATION==1) {
    blitz::Array<std::complex<double>, 2> J_dip, M_dip;
    blitz::Array<double, 2> r_J_dip, r_M_dip;
    int J_DIPOLES_EXCITATION, M_DIPOLES_EXCITATION;
    // electric dipoles
    readIntFromASCIIFile(V_CFIE_DATA_PATH + "J_DIPOLES_EXCITATION.txt", J_DIPOLES_EXCITATION);
    if (J_DIPOLES_EXCITATION==1) {
      readComplexDoubleBlitzArray2DFromASCIIFile( V_CFIE_DATA_PATH + "J_dip.txt", J_dip);
      readDoubleBlitzArray2DFromASCIIFile( V_CFIE_DATA_PATH + "r_J_dip.txt", r_J_dip);
      if (my_id==0) {
        //cout << "J_dip.txt = " << J_dip << endl;
        //cout << "r_J_dip.txt = " << r_J_dip << endl;
      }
      blitz::Array<std::complex<float>, 1> V_CFIE_tmp;
      const char CURRENT_TYPE = 'J';
      local_V_CFIE_dipole_array (V_CFIE_tmp, J_dip, r_J_dip, local_target_mesh, w, eps_r, mu_r, octtree.CFIE, CURRENT_TYPE, V_FULL_PRECISION);
      V_CFIE += V_CFIE_tmp;
    }
    // magnetic dipoles
    readIntFromASCIIFile(V_CFIE_DATA_PATH + "M_DIPOLES_EXCITATION.txt", M_DIPOLES_EXCITATION);
    if (M_DIPOLES_EXCITATION==1) {
      readComplexDoubleBlitzArray2DFromASCIIFile( V_CFIE_DATA_PATH + "M_dip.txt", M_dip);
      readDoubleBlitzArray2DFromASCIIFile( V_CFIE_DATA_PATH + "r_M_dip.txt", r_M_dip);
      if (my_id==0) {
        //cout << "M_dip.txt = " << M_dip << endl;
        //cout << "r_M_dip.txt = " << r_M_dip << endl;
      }
      blitz::Array<std::complex<float>, 1> V_CFIE_tmp;
      const char CURRENT_TYPE = 'M';
      local_V_CFIE_dipole_array (V_CFIE_tmp, J_dip, r_J_dip, local_target_mesh, w, eps_r, mu_r, octtree.CFIE, CURRENT_TYPE, V_FULL_PRECISION);
      V_CFIE += V_CFIE_tmp;
    }
  }
  if (PLANE_WAVE_EXCITATION==1) {
    double theta_inc, phi_inc;
    readDoubleFromASCIIFile(V_CFIE_DATA_PATH + "theta_inc.txt", theta_inc);
    readDoubleFromASCIIFile(V_CFIE_DATA_PATH + "phi_inc.txt", phi_inc);
    // local coordinate system
    blitz::Array<double, 1> r_hat(3), theta_hat(3), phi_hat(3), r_ref(3);
    r_hat = sin(theta_inc)*cos(phi_inc), sin(theta_inc)*sin(phi_inc), cos(theta_inc);
    theta_hat = cos(theta_inc)*cos(phi_inc), cos(theta_inc)*sin(phi_inc), -sin(theta_inc);
    phi_hat = -sin(phi_inc), cos(phi_inc), 0.0;
    blitz::Array<double, 1> k_hat(-1.0 * r_hat);
    r_ref = 0.0;
    // excitation electric field
    // we now read the incoming field amplitude and phase and polarization
    // only 2 components are needed: E_theta and E_phi
    blitz::Array<std::complex<double>, 1> E_inc_components(2), E_inc(3);
    readComplexDoubleBlitzArray1DFromASCIIFile( V_CFIE_DATA_PATH + "E_inc.txt", E_inc_components );
    E_inc = E_inc_components(0) * theta_hat + E_inc_components(1) * phi_hat;
    blitz::Array<std::complex<float>, 1> V_CFIE_tmp;
    local_V_CFIE_plane (V_CFIE_tmp, E_inc, k_hat, r_ref, local_target_mesh, octtree.w, octtree.eps_r, octtree.mu_r, octtree.CFIE, V_FULL_PRECISION);
    V_CFIE += V_CFIE_tmp;
  }
  local_target_mesh.resizeToZero();
  // iterative solving
  readDoubleFromASCIIFile(ITERATIVE_DATA_PATH + "TOL.txt", TOL);
  readIntFromASCIIFile(ITERATIVE_DATA_PATH + "RESTART.txt", RESTART);
  readIntFromASCIIFile(ITERATIVE_DATA_PATH + "MAXITER.txt", MAXITER);
  blitz::Array<std::complex<float>, 1> ZI(N_local_RWG);
  ZI = 0.0;
  if (SOLVER=="BICGSTAB") {
    PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > psolve(&leftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);
    bicgstab(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
  }
  else if (SOLVER=="GMRES") {
    PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > psolve(&leftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);
    gmres(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, RESTART, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
  } 
  else if (SOLVER=="RGMRES") {
    PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > psolve(&leftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);
    fgmres(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, RESTART, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
  } 
  else if (SOLVER=="FGMRES") {
    string INNER_SOLVER;
    readStringFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_SOLVER.txt", INNER_SOLVER);
    int INNER_MAXITER, INNER_RESTART;
    readIntFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_MAXITER.txt", INNER_MAXITER);
    readIntFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_RESTART.txt", INNER_RESTART);
    double INNER_TOL;
    readDoubleFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_TOL.txt", INNER_TOL);
    PsolveAMLFMA psolveAMLFMA(matvecMLFMA, leftFrobPsolveMLFMA, INNER_TOL, INNER_MAXITER, INNER_RESTART, N_RWG, INNER_SOLVER, SIMU_DIR);
    PrecondFunctor< std::complex<float>, PsolveAMLFMA > psolve(&psolveAMLFMA, &PsolveAMLFMA::psolve);
    fgmres(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, RESTART, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
  }
  else {
    cout << "Bad solver choice!! Solver is BICGSTAB or (F)GMRES, and you chose " << SOLVER << endl;
    exit(1);
  }
  octtree.resizeSdownLevelsToZero();

  // now computing the E_field at the r_obs
  blitz::Array<double, 2> r_obs;
  readDoubleBlitzArray2DFromASCIIFile( V_CFIE_DATA_PATH + "r_obs.txt", r_obs);
  blitz::Array<std::complex<double>, 2> E_obs;
  local_target_mesh.setLocalMeshFromFile(MESH_DATA_PATH);
  computeE_obs(E_obs, r_obs, local_target_mesh, ZI, eps_r, mu_r, w);
  local_target_mesh.resizeToZero();
  if (my_id==master) writeComplexDoubleBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "E_obs.txt", E_obs);
  if (my_id==master) writeDoubleBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "r_obs.txt", r_obs);

  // calculating the far fields
  blitz::Array<float, 1> octtreeXthetas_coarsest, octtreeXphis_coarsest;
  readFloatBlitzArray1DFromASCIIFile(OCTTREE_DATA_PATH + "octtreeXphis_coarsest.txt", octtreeXphis_coarsest);
  readFloatBlitzArray1DFromASCIIFile(OCTTREE_DATA_PATH + "octtreeXthetas_coarsest.txt", octtreeXthetas_coarsest);
  blitz::Array<std::complex<float>, 2> e_theta_far, e_phi_far;
  octtree.computeFarField(e_theta_far, e_phi_far, octtreeXthetas_coarsest, octtreeXphis_coarsest, ZI, OCTTREE_DATA_PATH);
  if (my_id==master) {
    writeComplexFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "e_theta_far_ASCII.txt", e_theta_far);
    writeComplexFloatBlitzArray2DToBinaryFile(RESULT_DATA_PATH + "e_theta_far_Binary.txt", e_theta_far);
    writeComplexFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "e_phi_far_ASCII.txt", e_phi_far);
    writeComplexFloatBlitzArray2DToBinaryFile(RESULT_DATA_PATH + "e_phi_far_Binary.txt", e_phi_far);
    // thetas, phis
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "phis_far_field_ASCII.txt", octtreeXphis_coarsest);
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "thetas_far_field_ASCII.txt", octtreeXthetas_coarsest);
  }
  // we now write the solution to disk...
  blitz::Array<std::complex<float>, 1> ZI_global(N_RWG), recvBuf;
  ZI_global = 0.0;
  if ( my_id == master ) {
    recvBuf.resize(N_RWG);
    recvBuf = 0.0;
  }
  for (int i=0 ; i<localRWGNumbers.size() ; ++i) ZI_global(localRWGNumbers(i)) = ZI(i);
  int ierror = MPI_Reduce(ZI_global.data(), recvBuf.data(), N_RWG, MPI_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
  if ( my_id == master ) {
    string filename = TMP + "/ZI/ZI.txt";
    ofstream ofs(filename.c_str(), blitz::ios::binary);
    ofs.write((char *)(recvBuf.data()), recvBuf.size()*8);
    ofs.close();
    writeIntToASCIIFile(ITERATIVE_DATA_PATH + "numberOfMatvecs.txt", octtree.getNumberOfUpdates());
    writeIntToASCIIFile(ITERATIVE_DATA_PATH + "iter.txt", iter);
    cout << endl;
  }
}

void computeMonostaticRCS(Octtree & octtree,
                          LocalMesh & local_target_mesh,
                          const string SOLVER,
                          const string SIMU_DIR,
                          const string TMP,
                          const string OCTTREE_DATA_PATH,
                          const string MESH_DATA_PATH,
                          const string V_CFIE_DATA_PATH,
                          const string RESULT_DATA_PATH,
                          const string ITERATIVE_DATA_PATH)
{
  blitz::Range all = blitz::Range::all();
  int num_procs = MPI::COMM_WORLD.Get_size(), my_id = MPI::COMM_WORLD.Get_rank();
  const int master = 0, N_local_RWG = local_target_mesh.N_local_RWG;
  
  const std::complex<float> eps_r = octtree.eps_r, mu_r = octtree.mu_r;
  blitz::Array<int, 1> localRWGNumbers(local_target_mesh.localRWGNumbers.size());
  localRWGNumbers = local_target_mesh.localRWGNumbers;

  // functors declarations
  const string PRECOND = "FROB";
  int N_RWG, iter, flag, RESTART, MAXITER, USE_PREVIOUS_SOLUTION, MONOSTATIC_BY_BISTATIC_APPROX, V_FULL_PRECISION;
  double error, TOL, MAX_DELTA_PHASE;
  readIntFromASCIIFile(V_CFIE_DATA_PATH + "V_FULL_PRECISION.txt", V_FULL_PRECISION);
  readIntFromASCIIFile(OCTTREE_DATA_PATH + "N_RWG.txt", N_RWG);
  readDoubleFromASCIIFile(ITERATIVE_DATA_PATH + "TOL.txt", TOL);
  readIntFromASCIIFile(ITERATIVE_DATA_PATH + "RESTART.txt", RESTART);
  readIntFromASCIIFile(ITERATIVE_DATA_PATH + "MAXITER.txt", MAXITER);
  readIntFromASCIIFile(TMP + "/USE_PREVIOUS_SOLUTION.txt", USE_PREVIOUS_SOLUTION);
  readIntFromASCIIFile(TMP + "/MONOSTATIC_BY_BISTATIC_APPROX.txt", MONOSTATIC_BY_BISTATIC_APPROX);
  readDoubleFromASCIIFile(TMP + "/MAXIMUM_DELTA_PHASE.txt", MAX_DELTA_PHASE);
  MAX_DELTA_PHASE *= M_PI/180.0;
  if (MONOSTATIC_BY_BISTATIC_APPROX!=1) MAX_DELTA_PHASE = 0.0;

  MatvecMLFMA matvecMLFMA(octtree, N_RWG, localRWGNumbers, SIMU_DIR);
  MatvecFunctor< std::complex<float>, MatvecMLFMA > matvec(&matvecMLFMA, &MatvecMLFMA::matvec);
  LeftFrobPsolveMLFMA leftFrobPsolveMLFMA(N_RWG, localRWGNumbers, SIMU_DIR);
  PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > psolve(&leftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);
  if (SOLVER=="FGMRES") {
    string INNER_SOLVER;
    readStringFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_SOLVER.txt", INNER_SOLVER);
    int INNER_MAXITER, INNER_RESTART;
    readIntFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_MAXITER.txt", INNER_MAXITER);
    readIntFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_RESTART.txt", INNER_RESTART);
    double INNER_TOL;
    readDoubleFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_TOL.txt", INNER_TOL);
    PsolveAMLFMA psolveAMLFMA(matvecMLFMA, leftFrobPsolveMLFMA, INNER_TOL, INNER_MAXITER, INNER_RESTART, N_RWG, INNER_SOLVER, SIMU_DIR);
    PrecondFunctor< std::complex<float>, PsolveAMLFMA > psolve(&psolveAMLFMA, &PsolveAMLFMA::psolve);
  }
  else PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > psolve(&leftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);
  // getting the angles at which monostatic RCS must be computed
  blitz::Array<float, 1> octtreeXthetas_coarsest, octtreeXphis_coarsest;
  readFloatBlitzArray1DFromASCIIFile(OCTTREE_DATA_PATH + "octtreeXphis_coarsest.txt", octtreeXphis_coarsest);
  readFloatBlitzArray1DFromASCIIFile(OCTTREE_DATA_PATH + "octtreeXthetas_coarsest.txt", octtreeXthetas_coarsest);
  const int N_theta(octtreeXthetas_coarsest.size()), N_phi(octtreeXphis_coarsest.size());
  blitz::Array<float, 2> RCS_VV(N_theta, N_phi), RCS_HH(N_theta, N_phi), RCS_HV(N_theta, N_phi), RCS_VH(N_theta, N_phi);
  RCS_VV = 1.0;
  RCS_HH = 1.0;
  RCS_HV = 1.0;
  RCS_VH = 1.0;
  int COMPUTE_RCS_HH, COMPUTE_RCS_HV, COMPUTE_RCS_VH, COMPUTE_RCS_VV;
  readIntFromASCIIFile(TMP + "/COMPUTE_RCS_HH.txt", COMPUTE_RCS_HH);
  readIntFromASCIIFile(TMP + "/COMPUTE_RCS_HV.txt", COMPUTE_RCS_HV);
  readIntFromASCIIFile(TMP + "/COMPUTE_RCS_VH.txt", COMPUTE_RCS_VH);
  readIntFromASCIIFile(TMP + "/COMPUTE_RCS_VV.txt", COMPUTE_RCS_VV);
  // r_ref
  blitz::Array<double, 1> r_ref(3);
  for (int i=0 ; i<3 ; ++i) r_ref(i) = octtree.big_cube_center_coord(i);
  // we now make the calculations of the parameters for the MONOSTATIC_BY_BISTATIC_APPROX
  // see pages 94-96 of Chew 2001, "Fast and Efficient Methods in Computational EM"
  const double big_cube_side_length = 2.0 * (octtree.big_cube_center_coord(0) - octtree.big_cube_lower_coord(0));
  const double LL = big_cube_side_length * sqrt(2.0);
  const double Delta_Phi = octtreeXphis_coarsest(1) - octtreeXphis_coarsest(0);
  double Beta = sqrt(4.0 * c * MAX_DELTA_PHASE / (LL * octtree.w));
  // number of phi points encompassed by Beta: number of points for which I can approximate
  // monostatic RCS by bistatic RCS
  const int BetaPoints = static_cast<int>(floor(Beta/Delta_Phi)) + 1;
  Beta = (BetaPoints-1) * Delta_Phi;
  if (my_id==master) cout << "number of BetaPoints for monostatic-bistatic approximation = " << BetaPoints << endl;
  // loop for monostatic sigma computation
  for (int excitation=0 ; excitation<2 ; ++excitation) { // 0 for H, 1 for V
    const bool HH = ((excitation==0) && (COMPUTE_RCS_HH==1));
    const bool HV = ((excitation==0) && (COMPUTE_RCS_HV==1));
    const bool VH = ((excitation==1) && (COMPUTE_RCS_VH==1));
    const bool VV = ((excitation==1) && (COMPUTE_RCS_VV==1));
    const bool cond = (HH || HV || VH || VV);
    if (cond) {
      for (int t=0 ; t<N_theta ; ++t) {
        blitz::Array<std::complex<float>, 1> ZI(N_local_RWG);
        ZI = 0.0;
        const float theta = octtreeXthetas_coarsest(t);
        float phi_inc = octtreeXphis_coarsest(0) + Beta/2.0;
        int startIndexPhi = 0;
        while (startIndexPhi<N_phi) {
          octtree.resizeSdownLevelsToZero();
          if (my_id==master) {
            if (HH || HV) cout << "\nHH and HV, theta = "<< theta * 180.0/M_PI << ", phi = " << phi_inc * 180.0/M_PI << endl;
            else cout << "\nVV and VH, theta = "<< theta * 180.0/M_PI << ", phi = " << phi_inc * 180.0/M_PI << endl;
            flush(cout);
          }
          // local coordinate system
          blitz::Array<double, 1> r_hat(3), theta_hat(3), phi_hat(3);
          r_hat = sin(theta)*cos(phi_inc), sin(theta)*sin(phi_inc), cos(theta);
          theta_hat = cos(theta)*cos(phi_inc), cos(theta)*sin(phi_inc), -sin(theta);
          phi_hat = -sin(phi_inc), cos(phi_inc), 0.0;
          blitz::Array<double, 1> k_hat(-1.0 * r_hat);

          // excitation field
          blitz::Array<std::complex<double>, 1> E_0(3);
          if (HH || HV) E_0 = 100.0 * (phi_hat + I * 0.0);
          else E_0 = -100.0 * (theta_hat + I * 0.0);
          blitz::Array<std::complex<float>, 1> V_CFIE;
          local_target_mesh.setLocalMeshFromFile(MESH_DATA_PATH);
          local_V_CFIE_plane (V_CFIE, E_0, k_hat, r_ref, local_target_mesh, octtree.w, octtree.eps_r, octtree.mu_r, octtree.CFIE, V_FULL_PRECISION);
          local_target_mesh.resizeToZero();
          // solving
          octtree.setNumberOfUpdates(0);
          if (USE_PREVIOUS_SOLUTION != 1) ZI = 0.0;
          if (SOLVER=="BICGSTAB") bicgstab(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
          else if (SOLVER=="GMRES") gmres(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, RESTART, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
          else if ((SOLVER=="RGMRES") || (SOLVER=="FGMRES")) fgmres(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, RESTART, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
          else {
            cout << "Bad solver choice!! Solver is BICGSTAB or (F)GMRES, and you chose " << SOLVER << endl;
            exit(1);
          }
          // far field computation
          blitz::Array<std::complex<float>, 2> e_theta_far, e_phi_far;
          blitz::Array<float, 1> thetas(1), phis(BetaPoints);
          thetas(0) = theta;
          if ((BetaPoints%2)!=0) { // if BetaPoints is odd and greater than 2
            const float space = 2.0*Delta_Phi;
            for (int j=0 ; j<BetaPoints ; ++j) phis(j) = phi_inc - Beta + j * space;
          }
          else {
            const float space = Delta_Phi;
            for (int j=0 ; j<BetaPoints/2 ; ++j) {
              phis(j) = phi_inc - Beta + j * space;
              phis(BetaPoints-1-j) = phi_inc + Beta - j * space;
            }
          }
          octtree.computeFarField(e_theta_far, e_phi_far, thetas, phis, ZI, OCTTREE_DATA_PATH);
          // filling of the RCS Arrays
          if (HH || HV) {
            for (int j=0 ; j<BetaPoints ; ++j) {
              RCS_HH(t, startIndexPhi + j) = real(e_phi_far(0, j) * conj(e_phi_far(0, j)))/real(sum(E_0 * conj(E_0)) * 4.0*M_PI);
              RCS_HV(t, startIndexPhi + j) = real(e_theta_far(0, j) * conj(e_theta_far(0, j)))/real(sum(E_0 * conj(E_0)) * 4.0*M_PI);
            }
          }
          else {
            for (int j=0 ; j<BetaPoints ; ++j) {
              RCS_VV(t, startIndexPhi + j) = real(e_theta_far(0, j) * conj(e_theta_far(0, j)))/real(sum(E_0 * conj(E_0)) * 4.0*M_PI);
              RCS_VH(t, startIndexPhi + j) = real(e_phi_far(0, j) * conj(e_phi_far(0, j)))/real(sum(E_0 * conj(E_0)) * 4.0*M_PI);
            }
          }
          // phi update
          startIndexPhi += BetaPoints;
          phi_inc += BetaPoints * Delta_Phi; // += Beta;
          if (my_id==master) {
            writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "RCS_HH_ASCII.txt", RCS_HH);
            writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "RCS_HV_ASCII.txt", RCS_HV);
            writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "RCS_VV_ASCII.txt", RCS_VV);
            writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "RCS_VH_ASCII.txt", RCS_VH);
          }
        }
      }
    }
  }
  if (my_id==master) {
    writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "RCS_HH_ASCII.txt", RCS_HH);
    writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "RCS_HV_ASCII.txt", RCS_HV);
    writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "RCS_VH_ASCII.txt", RCS_VH);
    writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "RCS_VV_ASCII.txt", RCS_VV);
    // thetas, phis
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "phis_far_field_ASCII.txt", octtreeXphis_coarsest);
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "thetas_far_field_ASCII.txt", octtreeXthetas_coarsest);
    // final writes
    writeIntToASCIIFile(ITERATIVE_DATA_PATH + "numberOfMatvecs.txt", octtree.getNumberOfUpdates());
    writeIntToASCIIFile(ITERATIVE_DATA_PATH + "iter.txt", iter);
  }
}

void computeMonostaticSAR(Octtree & octtree,
                          LocalMesh & local_target_mesh,
                          const string SOLVER,
                          const string SIMU_DIR,
                          const string TMP,
                          const string OCTTREE_DATA_PATH,
                          const string MESH_DATA_PATH,
                          const string V_CFIE_DATA_PATH,
                          const string RESULT_DATA_PATH,
                          const string ITERATIVE_DATA_PATH)
{
  blitz::Range all = blitz::Range::all();
  int num_procs = MPI::COMM_WORLD.Get_size(), my_id = MPI::COMM_WORLD.Get_rank();
  const int master = 0, N_local_RWG = local_target_mesh.N_local_RWG;
  const double w = octtree.w;
  blitz::Array<int, 1> localRWGNumbers(local_target_mesh.localRWGNumbers.size());
  localRWGNumbers = local_target_mesh.localRWGNumbers;
  const std::complex<double> eps_r = octtree.eps_r, mu_r = octtree.mu_r;
  const std::complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu);

  // functors declarations
  const string PRECOND = "FROB";
  int N_RWG, iter, flag, RESTART, MAXITER, USE_PREVIOUS_SOLUTION, MONOSTATIC_BY_BISTATIC_APPROX, V_FULL_PRECISION;
  double error, TOL, MAX_DELTA_PHASE;
  readIntFromASCIIFile(V_CFIE_DATA_PATH + "V_FULL_PRECISION.txt", V_FULL_PRECISION);
  readIntFromASCIIFile(OCTTREE_DATA_PATH + "N_RWG.txt", N_RWG);
  readDoubleFromASCIIFile(ITERATIVE_DATA_PATH + "TOL.txt", TOL);
  readIntFromASCIIFile(ITERATIVE_DATA_PATH + "RESTART.txt", RESTART);
  readIntFromASCIIFile(ITERATIVE_DATA_PATH + "MAXITER.txt", MAXITER);
  readIntFromASCIIFile(TMP + "/USE_PREVIOUS_SOLUTION.txt", USE_PREVIOUS_SOLUTION);

  MatvecMLFMA matvecMLFMA(octtree, N_RWG, localRWGNumbers, SIMU_DIR);
  MatvecFunctor< std::complex<float>, MatvecMLFMA > matvec(&matvecMLFMA, &MatvecMLFMA::matvec);
  LeftFrobPsolveMLFMA leftFrobPsolveMLFMA(N_RWG, localRWGNumbers, SIMU_DIR);
  PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > psolve(&leftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);
  if (SOLVER=="FGMRES") {
    string INNER_SOLVER;
    readStringFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_SOLVER.txt", INNER_SOLVER);
    int INNER_MAXITER, INNER_RESTART;
    readIntFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_MAXITER.txt", INNER_MAXITER);
    readIntFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_RESTART.txt", INNER_RESTART);
    double INNER_TOL;
    readDoubleFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_TOL.txt", INNER_TOL);
    PsolveAMLFMA psolveAMLFMA(matvecMLFMA, leftFrobPsolveMLFMA, INNER_TOL, INNER_MAXITER, INNER_RESTART, N_RWG, INNER_SOLVER, SIMU_DIR);
    PrecondFunctor< std::complex<float>, PsolveAMLFMA > psolve(&psolveAMLFMA, &PsolveAMLFMA::psolve);
  }
  else PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > psolve(&leftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);
  // getting the positions at which monostatic SAR must be computed
  blitz::Array<double, 1> SAR_local_x_hat(3), SAR_local_y_hat(3), SAR_plane_origin(3);
  readDoubleBlitzArray1DFromASCIIFile(V_CFIE_DATA_PATH + "SAR_local_x_hat.txt", SAR_local_x_hat);
  readDoubleBlitzArray1DFromASCIIFile(V_CFIE_DATA_PATH + "SAR_local_y_hat.txt", SAR_local_y_hat);
  // normalization
  SAR_local_x_hat = SAR_local_x_hat / sqrt(sum(SAR_local_x_hat * SAR_local_x_hat));
  SAR_local_y_hat = SAR_local_y_hat / sqrt(sum(SAR_local_y_hat * SAR_local_y_hat));
  // other SAR data
  readDoubleBlitzArray1DFromASCIIFile(V_CFIE_DATA_PATH + "SAR_plane_origin.txt", SAR_plane_origin);
  double SAR_x_span, SAR_x_span_offset, SAR_y_span, SAR_y_span_offset;
  readDoubleFromASCIIFile(V_CFIE_DATA_PATH + "SAR_x_span.txt", SAR_x_span);
  readDoubleFromASCIIFile(V_CFIE_DATA_PATH + "SAR_x_span_offset.txt", SAR_x_span_offset);
  readDoubleFromASCIIFile(V_CFIE_DATA_PATH + "SAR_y_span.txt", SAR_y_span);
  readDoubleFromASCIIFile(V_CFIE_DATA_PATH + "SAR_y_span_offset.txt", SAR_y_span_offset);
  int SAR_N_x_points, SAR_N_y_points;
  readIntFromASCIIFile(V_CFIE_DATA_PATH + "SAR_N_x_points.txt", SAR_N_x_points);
  readIntFromASCIIFile(V_CFIE_DATA_PATH + "SAR_N_y_points.txt", SAR_N_y_points);
  // protection against stupidity
  if (SAR_N_x_points<1) SAR_N_x_points = 1;
  if (SAR_N_y_points<1) SAR_N_y_points = 1;
  // step increments
  double Delta_x, Delta_y;
  Delta_x = (SAR_N_x_points>1) ? SAR_x_span / (SAR_N_x_points-1) : SAR_x_span/2.0;
  Delta_y = (SAR_N_y_points>1) ? SAR_y_span / (SAR_N_y_points-1) : SAR_y_span/2.0;

  blitz::Array<float, 2> r_SAR(SAR_N_x_points * SAR_N_y_points, 3);
  blitz::Array<float, 1> RCS_VV(SAR_N_x_points * SAR_N_y_points), RCS_HH(SAR_N_x_points * SAR_N_y_points), RCS_HV(SAR_N_x_points * SAR_N_y_points), RCS_VH(SAR_N_x_points * SAR_N_y_points);
  r_SAR = 0.0;
  RCS_VV = 1.0;
  RCS_HH = 1.0;
  RCS_HV = 1.0;
  RCS_VH = 1.0;
  int COMPUTE_RCS_HH, COMPUTE_RCS_HV, COMPUTE_RCS_VH, COMPUTE_RCS_VV;
  readIntFromASCIIFile(TMP + "/COMPUTE_RCS_HH.txt", COMPUTE_RCS_HH);
  readIntFromASCIIFile(TMP + "/COMPUTE_RCS_HV.txt", COMPUTE_RCS_HV);
  readIntFromASCIIFile(TMP + "/COMPUTE_RCS_VH.txt", COMPUTE_RCS_VH);
  readIntFromASCIIFile(TMP + "/COMPUTE_RCS_VV.txt", COMPUTE_RCS_VV);
  // r_ref
  blitz::Array<double, 1> r_ref(3);
  for (int i=0 ; i<3 ; ++i) r_ref(i) = octtree.big_cube_center_coord(i);
  // loop for monostatic sigma computation
  for (int excitation=0 ; excitation<2 ; ++excitation) { // 0 for H, 1 for V
    const bool HH = ((excitation==0) && (COMPUTE_RCS_HH==1));
    const bool HV = ((excitation==0) && (COMPUTE_RCS_HV==1));
    const bool VH = ((excitation==1) && (COMPUTE_RCS_VH==1));
    const bool VV = ((excitation==1) && (COMPUTE_RCS_VV==1));
    const bool cond = (HH || HV || VV);
    if (cond) {
      for (int t=0 ; t<SAR_N_y_points ; ++t) {
        blitz::Array<std::complex<float>, 1> ZI(N_local_RWG);
        ZI = 0.0;
        for (int s=0 ; s<SAR_N_x_points ; ++s) {
          int index = t * SAR_N_x_points + s;
          blitz::Array<double, 1> r_src(3);
          r_src = SAR_plane_origin + SAR_x_span_offset * SAR_local_x_hat + SAR_y_span_offset * SAR_local_y_hat;
          r_src += s * Delta_x * SAR_local_x_hat + t * Delta_y * SAR_local_y_hat;
          r_SAR(index, all) = 1.0 * r_src;
          octtree.resizeSdownLevelsToZero();
          if (my_id==master) {
            if (HH || HV) blitz::cout << "\nHH and HV, r_ant = "<< r_src << blitz::endl;
            else blitz::cout << "\nVV and VH, r_ant = "<< r_src << blitz::endl;
            blitz::flush(blitz::cout);
          }
          // excitation field
          blitz::Array<std::complex<double>, 1> J_ant(3), E_0(3);
          if (HH || HV) J_ant = 100. * SAR_local_x_hat;
          else J_ant = 100. * SAR_local_y_hat;
          blitz::Array<std::complex<float>, 1> V_CFIE;
          local_target_mesh.setLocalMeshFromFile(MESH_DATA_PATH);
          local_V_CFIE_dipole (V_CFIE, J_ant, r_src, local_target_mesh, w, eps_r, mu_r, octtree.CFIE, V_FULL_PRECISION);
          local_target_mesh.resizeToZero();
          // target incoming field: reference field for RCS computation
          blitz::Array<std::complex<double>, 2> G_EJ(3, 3), G_HJ(3, 3);
          double r_dip[3], r_obs2[3];
          for (int m=0 ; m<3 ; m++) {
            r_dip[m] = r_src(m);
            r_obs2[m] = r_ref(m);
          }
          G_EJ_G_HJ (G_EJ, G_HJ, r_dip, r_obs2, eps, mu, k);
          for (int m=0 ; m<3 ; m++) E_0(m) = G_EJ(m, 0) * J_ant(0) + G_EJ(m, 1) * J_ant(1) + G_EJ(m, 2) * J_ant(2);
          // solving
          octtree.setNumberOfUpdates(0);
          if (USE_PREVIOUS_SOLUTION != 1) ZI = 0.0;
          if (SOLVER=="BICGSTAB") bicgstab(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
          else if (SOLVER=="GMRES") gmres(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, RESTART, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
          else if ((SOLVER=="RGMRES") || (SOLVER=="FGMRES")) fgmres(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, RESTART, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
          else {
            cout << "Bad solver choice!! Solver is BICGSTAB or (F)GMRES, and you chose " << SOLVER << endl;
            exit(1);
          }
          // field computation (O(N) version)
          blitz::Array<double, 2> r_obs(1, 3);
          blitz::Array<std::complex<double>, 2> E_obs(1, 3);
          r_obs(0, all) = r_ref;
          local_target_mesh.setLocalMeshFromFile(MESH_DATA_PATH);
          computeE_obs(E_obs, r_obs, local_target_mesh, ZI, static_cast<complex<float> >(eps_r), static_cast<complex<float> >(mu_r), w);
          local_target_mesh.resizeToZero();
          // filling of the RCS Arrays
          if (HH || HV) {
            blitz::Array<std::complex<double>, 1> E_H(3), E_V(3);
            E_H = sum(SAR_local_x_hat * E_obs(0, all));
            E_V = sum(SAR_local_y_hat * E_obs(0, all));
            RCS_HH(index) = real(sum(E_H * conj(E_H)))/real(sum(E_0 * conj(E_0)) * 4.0*M_PI);
            RCS_HV(index) = real(sum(E_V * conj(E_V)))/real(sum(E_0 * conj(E_0)) * 4.0*M_PI);
          }
          else {
            blitz::Array<std::complex<double>, 1> E_H(3), E_V(3);
            E_H = sum(SAR_local_x_hat * E_obs(0, all));
            E_V = sum(SAR_local_y_hat * E_obs(0, all));
            RCS_VV(index) = real(sum(E_V * conj(E_V)))/real(sum((E_0 * conj(E_0))) * 4.0*M_PI);
            RCS_VH(index) = real(sum(E_H * conj(E_H)))/real(sum(E_0 * conj(E_0)) * 4.0*M_PI);
          }
          if (my_id==master) {
            writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "SAR_RCS_HH_ASCII.txt", RCS_HH);
            writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "SAR_RCS_HV_ASCII.txt", RCS_HV);
            writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "SAR_RCS_VH_ASCII.txt", RCS_VH);
            writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "SAR_RCS_VV_ASCII.txt", RCS_VV);
            writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "r_SAR.txt", r_SAR);
          }
        }
      }
    }
  }
  // final writings
  if (my_id==master) {
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "SAR_RCS_HH_ASCII.txt", RCS_HH);
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "SAR_RCS_HV_ASCII.txt", RCS_HV);
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "SAR_RCS_VH_ASCII.txt", RCS_VH);
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "SAR_RCS_VV_ASCII.txt", RCS_VV);
    writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "r_SAR.txt", r_SAR);
    writeIntToASCIIFile(ITERATIVE_DATA_PATH + "numberOfMatvecs.txt", octtree.getNumberOfUpdates());
    writeIntToASCIIFile(ITERATIVE_DATA_PATH + "iter.txt", iter);
    cout << endl;
  }
}

/****************************************************************************/
/******************************* main ***************************************/
/****************************************************************************/

int main(int argc, char* argv[]) {

  MPI::Init();
  int ierror;
  const int num_procs = MPI::COMM_WORLD.Get_size();
  const int my_id = MPI::COMM_WORLD.Get_rank();
  const int master = 0;
  MPI_Status status;

  string simuDir = ".";
  if ( argc > 2 ) {
     if( string(argv[1]) == "--simudir" ) simuDir = argv[2];
  }

  // general variables
  const string SIMU_DIR = simuDir;
  const string RESULT_DATA_PATH = SIMU_DIR + "/result/";
  const string TMP = SIMU_DIR + "/tmp" + intToString(my_id);
  const string OCTTREE_DATA_PATH = TMP + "/octtree_data/";
  const string MESH_DATA_PATH = TMP + "/mesh/";
  const string V_CFIE_DATA_PATH = TMP + "/V_CFIE/";
  const string ITERATIVE_DATA_PATH = TMP + "/iterative_data/";

  // reading and broadcasting mesh data
  int C, N_RWG, V;
  blitz::Array<int, 1> RWGNumber_CFIE_OK, cubes_RWGsNumbers, cube_N_RWGs, cube_startIndex_RWGs;
  blitz::Array<int, 2> RWGNumber_edgeVertexes, RWGNumber_oppVertexes;
  blitz::Array<double, 2> cubes_centroids, vertexes_coord;
  if (my_id==0)
  {
    // reading the arrays sizes
    string filename = MESH_DATA_PATH + "C.txt";
    readIntFromASCIIFile(filename, C);
    
    filename = MESH_DATA_PATH + "N_RWG.txt";
    readIntFromASCIIFile(filename, N_RWG);
    
    filename = MESH_DATA_PATH + "V.txt";
    readIntFromASCIIFile(filename, V);

    // resizing
    vertexes_coord.resize(V, 3);
    cube_N_RWGs.resize(C);
    cube_startIndex_RWGs.resize(C);
    cubes_RWGsNumbers.resize(N_RWG);
    cubes_centroids.resize(C, 3);
    RWGNumber_edgeVertexes.resize(N_RWG, 2);
    RWGNumber_oppVertexes.resize(N_RWG, 2);
    RWGNumber_CFIE_OK.resize(N_RWG);

    // reading the arrays
    filename = MESH_DATA_PATH + "cubes_centroids.txt";
    readDoubleBlitzArray2DFromBinaryFile(filename, cubes_centroids);
    
    filename = MESH_DATA_PATH + "cubes_RWGsNumbers.txt";
    readIntBlitzArray1DFromBinaryFile(filename, cubes_RWGsNumbers);
    
    filename = MESH_DATA_PATH + "cube_N_RWGs.txt";
    readIntBlitzArray1DFromBinaryFile(filename, cube_N_RWGs);
    
    // creation of cube_startIndex_RWGs
    int startIndex = 0;
    for (int i=0; i<C; i++) {
      cube_startIndex_RWGs(i) = startIndex;
      startIndex += cube_N_RWGs(i);
    }
    
    filename = MESH_DATA_PATH + "vertexes_coord.txt";
    readDoubleBlitzArray2DFromBinaryFile(filename, vertexes_coord);

    filename = MESH_DATA_PATH + "RWGNumber_edgeVertexes.txt";
    readIntBlitzArray2DFromBinaryFile(filename, RWGNumber_edgeVertexes);

    filename = MESH_DATA_PATH + "RWGNumber_oppVertexes.txt";
    readIntBlitzArray2DFromBinaryFile(filename, RWGNumber_oppVertexes);

    filename = MESH_DATA_PATH + "RWGNumber_CFIE_OK.txt";
    readIntBlitzArray1DFromBinaryFile(filename, RWGNumber_CFIE_OK);
  }
  MPI_Bcast(&N_RWG, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&C, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (my_id!=0) cubes_centroids.resize(C, 3);
  MPI_Bcast(cubes_centroids.data(), cubes_centroids.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // Octtree creation based upon the cubes_centroids
  Octtree octtree(OCTTREE_DATA_PATH, cubes_centroids, my_id, num_procs);

  // partitioning of mesh
  blitz::Array<int, 1> oldIndexesOfCubes;
  octtree.computeIndexesOfCubesInOriginalMesh(oldIndexesOfCubes);
  int N_local_cubes = oldIndexesOfCubes.size();
  
  blitz::Array<int, 1> NumberOfCubesPerProcess;
  if (my_id==0) NumberOfCubesPerProcess.resize(num_procs);
  ierror = MPI_Gather(&N_local_cubes, 1, MPI_INT, NumberOfCubesPerProcess.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  blitz::Array<int, 1> MPI_Gatherv_scounts, MPI_Gatherv_displs; // it only matters for process 0
  if (my_id==0) {
    MPI_Gatherv_scounts.resize(num_procs);
    MPI_Gatherv_displs.resize(num_procs);
    int displacement = 0;
    for (int i=0 ; i<num_procs ; ++i) {
      MPI_Gatherv_scounts(i) = NumberOfCubesPerProcess(i);
      MPI_Gatherv_displs(i) = displacement;
      displacement += MPI_Gatherv_scounts(i);
    }
  }
  // process_OldIndexesOfCubes represents the cubes original numbers in mesh for each process
  // Only master process has complete overview of these cubes
  blitz::Array<int, 1> process_OldIndexesOfCubes;
  if (my_id==0) process_OldIndexesOfCubes.resize(C);
  ierror = MPI_Gatherv( oldIndexesOfCubes.data(), oldIndexesOfCubes.size(), MPI_INT, process_OldIndexesOfCubes.data(), MPI_Gatherv_scounts.data(), MPI_Gatherv_displs.data(), MPI_INT, 0,  MPI_COMM_WORLD);
  oldIndexesOfCubes.free();
  // now we need to scatter the RWG numbers and arrays, according to their original cubes, to each process
  // we have to begin by the end, to retain the good values for process 0!

  LocalMesh local_target_mesh;
  blitz::Array<int, 1> local_cubes_NRWG;
  if (my_id==master) {
    for (int receive_id=num_procs-1; receive_id>-1; receive_id--) {
      // first we find the RWGs for each cube of process receive_id
      const int startIndexOfCube = MPI_Gatherv_displs(receive_id);
      const int Ncubes = MPI_Gatherv_scounts(receive_id);
      std::vector<int> RWG_numbersTmp;
      blitz::Array<int, 1> local_cubes_NRWG_tmp(Ncubes);
      // we now find the RWG numbers for each cube of process receive_id
      for (int i=0; i<Ncubes; i++) {
        const int cubeNumber = process_OldIndexesOfCubes(startIndexOfCube + i);
        int N_RWG_in_cube = cube_N_RWGs(cubeNumber);
        int startIndex = cube_startIndex_RWGs(cubeNumber);
        for (int j=startIndex ; j<startIndex + N_RWG_in_cube ; ++j) RWG_numbersTmp.push_back(cubes_RWGsNumbers(j));
        local_cubes_NRWG_tmp(i) = N_RWG_in_cube;
      }
      // now creating a blitz::array from std::vector RWG_numbersTmp
      blitz::Array<int, 1> localRWGNumbers, localRWGNumber_CFIE_OK;
      blitz::Array<float, 2> localRWGNumber_trianglesCoord;
      localRWGNumbers.resize(RWG_numbersTmp.size());
      for (int i=0; i<localRWGNumbers.size(); i++) localRWGNumbers(i) = RWG_numbersTmp[i];
      RWG_numbersTmp.clear();
      // RWG_numbers and cube_NRWG contains all the info about the repartition of RWGs in cubes of process number="receive_id"
      // we've got to communicate these arrays, and also an array containing the coordinates of the RWGs.
      localRWGNumber_CFIE_OK.resize(localRWGNumbers.size());
      localRWGNumber_trianglesCoord.resize(localRWGNumbers.size(), 12);
      for (int i=0; i<localRWGNumbers.size(); i++) {
        const int RWG = localRWGNumbers(i);
        const int n0 = RWGNumber_oppVertexes(RWG, 0);
        const int n1 = RWGNumber_edgeVertexes(RWG, 0);
        const int n2 = RWGNumber_edgeVertexes(RWG, 1);
        const int n3 = RWGNumber_oppVertexes(RWG, 1);
        for (int j=0; j<3; j++) localRWGNumber_trianglesCoord(i, j) = vertexes_coord(n0, j);
        for (int j=0; j<3; j++) localRWGNumber_trianglesCoord(i, j+3) = vertexes_coord(n1, j);
        for (int j=0; j<3; j++) localRWGNumber_trianglesCoord(i, j+6) = vertexes_coord(n2, j);
        for (int j=0; j<3; j++) localRWGNumber_trianglesCoord(i, j+9) = vertexes_coord(n3, j);
        localRWGNumber_CFIE_OK(i) = RWGNumber_CFIE_OK(RWG);
      }
      if (receive_id!=master) {
        MPI_Send(local_cubes_NRWG_tmp.data(), local_cubes_NRWG_tmp.size(), MPI_INT, receive_id, receive_id, MPI_COMM_WORLD);
        MPI_Send(localRWGNumbers.data(), localRWGNumbers.size(), MPI_INT, receive_id, receive_id+1, MPI_COMM_WORLD);
        MPI_Send(localRWGNumber_CFIE_OK.data(), localRWGNumber_CFIE_OK.size(), MPI_INT, receive_id, receive_id+2, MPI_COMM_WORLD);
        cout << "send to P" << receive_id << " : " << "localRWGNumber_CFIE_OK.size() = " << localRWGNumber_CFIE_OK.size() << endl;
        flush(cout);
        MPI_Send(localRWGNumber_trianglesCoord.data(), localRWGNumber_trianglesCoord.size(), MPI_FLOAT, receive_id, receive_id+3, MPI_COMM_WORLD);
        cout << "send to P" << receive_id << " : " << "localRWGNumber_trianglesCoord.size() = " << localRWGNumber_trianglesCoord.size() << endl;
        flush(cout);
      }
      else {
        local_cubes_NRWG.resize(N_local_cubes);
        local_cubes_NRWG = local_cubes_NRWG_tmp;
        const int NRWG = sum(local_cubes_NRWG);
        local_target_mesh.N_local_RWG = NRWG;
        local_target_mesh.localRWGNumbers.resize(NRWG);
        local_target_mesh.localRWGNumber_CFIE_OK.resize(NRWG);
        local_target_mesh.localRWGNumber_trianglesCoord.resize(NRWG, 12);
        local_target_mesh.localRWGNumbers = localRWGNumbers;
        local_target_mesh.localRWGNumber_CFIE_OK = localRWGNumber_CFIE_OK;
        local_target_mesh.localRWGNumber_trianglesCoord = localRWGNumber_trianglesCoord;
      }
    }
  } // end of distribution loop
  if (my_id!=master) {
    local_cubes_NRWG.resize(N_local_cubes);
    MPI_Recv(local_cubes_NRWG.data(), local_cubes_NRWG.size(), MPI_INT, 0, my_id, MPI_COMM_WORLD, &status);
    const int NRWG = sum(local_cubes_NRWG);
    local_target_mesh.N_local_RWG = NRWG;
    local_target_mesh.localRWGNumbers.resize(NRWG);
    local_target_mesh.localRWGNumber_CFIE_OK.resize(NRWG);
    local_target_mesh.localRWGNumber_trianglesCoord.resize(NRWG, 12);
      
    MPI_Recv(local_target_mesh.localRWGNumbers.data(), NRWG, MPI_INT, 0, my_id+1, MPI_COMM_WORLD, &status);
    MPI_Recv(local_target_mesh.localRWGNumber_CFIE_OK.data(), NRWG, MPI_INT, 0, my_id+2, MPI_COMM_WORLD, &status);
    cout << "receive to P" << my_id << " : " << "N_local_RWGs = " << local_target_mesh.localRWGNumber_CFIE_OK.size() << endl;
    flush(cout);
    MPI_Recv(local_target_mesh.localRWGNumber_trianglesCoord.data(), local_target_mesh.localRWGNumber_trianglesCoord.size(), MPI_FLOAT, 0, my_id+3, MPI_COMM_WORLD, &status);
    cout << "receive to P" << my_id << " : " << "local_target_mesh.localRWGNumber_trianglesCoord.size() = " << local_target_mesh.localRWGNumber_trianglesCoord.size() << endl;
    flush(cout);
  }
     
  ierror = MPI_Barrier(MPI_COMM_WORLD);

  // now let's construct the octtree cubes local meshes!
  octtree.computeGaussLocatedArguments(local_cubes_NRWG, local_target_mesh.localRWGNumbers, local_target_mesh.localRWGNumber_CFIE_OK, local_target_mesh.localRWGNumber_trianglesCoord);
  octtree.RWGs_renumbering();
  // final moves
  local_target_mesh.reallyLocalRWGNumbers.resize(local_target_mesh.localRWGNumbers.size());
  for (int i=0 ; i<local_target_mesh.localRWGNumbers.size() ; ++i) local_target_mesh.reallyLocalRWGNumbers(i) = i;
  local_target_mesh.writeLocalMeshToFile(MESH_DATA_PATH);
  // OK, what kind of simulation do we want to run?
  octtree.constructArrays();
  int BISTATIC, MONOSTATIC_RCS, MONOSTATIC_SAR;
  readIntFromASCIIFile(TMP + "/BISTATIC.txt", BISTATIC);
  readIntFromASCIIFile(TMP + "/MONOSTATIC_RCS.txt", MONOSTATIC_RCS);
  readIntFromASCIIFile(TMP + "/MONOSTATIC_SAR.txt", MONOSTATIC_SAR);
  string SOLVER;
  readStringFromASCIIFile(ITERATIVE_DATA_PATH + "SOLVER.txt", SOLVER);
  if (my_id==0) cout << "SOLVER IS = " << SOLVER << endl;
  ierror = MPI_Barrier(MPI_COMM_WORLD);
  // bistatic computation
  if (BISTATIC==1) computeForOneExcitation(octtree, local_target_mesh, SOLVER, SIMU_DIR, TMP, OCTTREE_DATA_PATH, MESH_DATA_PATH, V_CFIE_DATA_PATH, RESULT_DATA_PATH, ITERATIVE_DATA_PATH);
  // monostatic RCS computation
  local_target_mesh.setLocalMeshFromFile(MESH_DATA_PATH);
  if (MONOSTATIC_RCS==1) computeMonostaticRCS(octtree, local_target_mesh, SOLVER, SIMU_DIR, TMP, OCTTREE_DATA_PATH, MESH_DATA_PATH, V_CFIE_DATA_PATH, RESULT_DATA_PATH, ITERATIVE_DATA_PATH);
  // monostatic SAR computation
  local_target_mesh.setLocalMeshFromFile(MESH_DATA_PATH);
  if (MONOSTATIC_SAR==1) computeMonostaticSAR(octtree, local_target_mesh, SOLVER, SIMU_DIR, TMP, OCTTREE_DATA_PATH, MESH_DATA_PATH, V_CFIE_DATA_PATH, RESULT_DATA_PATH, ITERATIVE_DATA_PATH);

  // Get peak memory usage of each rank
  long memusage_local = MemoryUsageGetPeak();
  std::cout << "MEMINFO " << argv[0] << " rank " << my_id << " mem=" << memusage_local/(1024*1024) << " MB" << std::endl;
  flush(std::cout);
  MPI::Finalize();
  return 0;
}
